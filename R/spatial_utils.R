library(sf)
library(tibble)
library(cppRouting)
library(dplyr)
library(lwgeom)
library(tidyr)
library(stringr)

# ------ Project to NAD83 Coloado Central StatePlane CRS ------
project_to_stateplane <- function(x, epsg = 26954) {
  stopifnot(inherits(x, "sf"))
  
  current_crs <- sf::st_crs(x)
  
  # If CRS missing entirely → error
  if (is.na(current_crs)) {
    stop("Input sf object has no CRS defined.")
  }
  
  # If EPSG missing OR different → transform
  if (is.na(current_crs$epsg) || current_crs$epsg != epsg) {
    x <- sf::st_transform(x, epsg)
  }
  
  return(x)
}

# ------ Assign event study cohorts ------

assign_groups <- function(data,
                          bgs,
                          stations,
                          type = c("treatment", "spillover", "continuous")){
  
  type <- match.arg(type)
  centroids <- sf::st_centroid(bgs)

  # Find nearest station index
  nearest_idx <- sf::st_nearest_feature(centroids, stations)
  
  # Compute distance (by element)
  dist_m <- as.numeric(
    sf::st_distance(
      centroids,
      stations[nearest_idx, ],
      by_element = TRUE
    )
  )
  
  # Add station info + distance
  bgs <- bgs %>%
    dplyr::mutate(
      nearest_station_id = stations$NAME[nearest_idx],
      cohort_year = stations$year[nearest_idx],
      distance_m = dist_m
    ) %>%
    dplyr::mutate(
      group = dplyr::case_when(
        distance_m <= 1000 ~ "treatment",
        distance_m > 1000 & distance_m <= 2000 ~ "exclude",
        distance_m > 2000 & distance_m <= 5000 ~ "spillover",
        distance_m > 5000 ~ "control"
      )
    )
  
  if(type == "treatment"){
    bgs <- bgs %>%
      mutate(
        cohort_year = ifelse(group == "treatment", cohort_year, 0)
      ) %>%
      filter(group != "exclude") # Exclude bgs 1-2 km from rail
  } else if(type == "spillover") {
    bgs <- bgs %>%
      mutate(
        cohort_year = ifelse(group == "spillover", cohort_year, 0)
      ) %>%
      filter(group %in% c("spillover","control")) # Exclude bgs under 2 km from rail
  } else {
    bgs <- bgs %>%
      mutate(
        log_dist = -log(distance_m/1000) # Negative log to show expected inverse relationship between dist and effects
      )
  }
  
  # Join dataset to assigned study bgs
  bgs_data <- bgs %>%
    sf::st_drop_geometry() %>%
    left_join(data, by = "GISJOIN")
  
  if(type == "continuous") {
    bgs_data <- bgs_data %>%
      mutate(ref_year = year - cohort_year) # Add in continuous interaction term
  } else {
    bgs_data <- bgs_data %>%
      filter(!is.na(Shape_Area)) # Filter excluded datasets
  }
  
  return(bgs_data)
}

# ------- Tract-to-bg crosswalk -------

make_tract_to_bg_xwalk <- function(bg_geoids) {
  # GEOID structure: SS CCC TTTTTT BBBB
  #                  2  3   6      4
  tibble(bg_geoid = bg_geoids) |>
    mutate(tract_geoid = str_sub(bg_geoid, 1, 11))  # state(2)+county(3)+tract(6)
}

apply_tract_shares_to_bg <- function(tract_shares, bg_geoids) {
  xwalk <- make_tract_to_bg_xwalk(bg_geoids)
  
  xwalk |>
    left_join(tract_shares, by = "tract_geoid") |>
    select(bg_geoid, mode, car_stratum, S)
}

# ------ Create road network edges ------

make_road_edges <- function(lines, weight_field, mode_tag,
                            exclude_fclass = NULL) {
  if (mode_tag == "cycle") {
    lines <- lines %>%
      mutate(Cycle_TT = Shape_Leng / 250)
  }
  
  edges <- lines_to_edges_any_vertex(lines, weight_field, mode_tag)
  
  if (!is.null(exclude_fclass)) {
    edges <- edges %>%
      filter(!fclass %in% exclude_fclass)
  }
  
  return(edges)
}

# ------ Create transit network edges ------

make_transit_edges <- function(lines, weight_field, mode_tag,
                               bidirectional = TRUE) {
  edges <- lines_to_edges_endpoint(lines, weight_field, mode_tag)
  
  # Transit is bidirectional by default
  if (bidirectional) {
    edges_rev <- edges %>%
      rename(from_id = to_id, to_id = from_id)
    edges <- bind_rows(edges, edges_rev)
  }
  
  return(edges)
}

# ------ Edge base functions ------

lines_to_edges_endpoint <- function(lines, weight_field, mode_tag) {
  
  edges <- lines %>%
    st_drop_geometry() %>%
    mutate(
      from_id = paste0(
        round(map_dbl(st_geometry(lines),
                      ~st_coordinates(.x)[1, 1]), 2), "_",
        round(map_dbl(st_geometry(lines),
                      ~st_coordinates(.x)[1, 2]), 2)
      ),
      to_id = paste0(
        round(map_dbl(st_geometry(lines),
                      ~st_coordinates(.x)[nrow(st_coordinates(.x)), 1]), 2), "_",
        round(map_dbl(st_geometry(lines),
                      ~st_coordinates(.x)[nrow(st_coordinates(.x)), 2]), 2)
      ),
      weight = .data[[weight_field]],
      mode   = mode_tag
    ) %>%
    # Keep all original columns including year_opera and Counterfac
    select(from_id, to_id, weight, mode, everything())
  
  return(edges)
}

lines_to_edges_any_vertex <- function(lines, weight_field, mode_tag) {
  
  # Get walk speed from the data itself
  # TravelTime = length / speed, so speed = length / TravelTime
  # Then apply to each sub-segment: sub_tt = sub_length / speed
  
  coords <- st_coordinates(lines) %>%
    as.data.frame() %>%
    rename(x = X, y = Y, line_id = L1)
  
  # Get full line length and travel time per feature
  line_attrs <- lines %>%
    st_drop_geometry() %>%
    mutate(
      line_id    = row_number(),
      full_tt    = .data[[weight_field]],
      full_len   = as.numeric(st_length(lines)),
      # Implied speed in meters per minute
      speed_mpm  = full_len / full_tt
    ) %>%
    select(line_id, full_tt, full_len, speed_mpm)
  
  edges <- coords %>%
    group_by(line_id) %>%
    mutate(
      from_id    = paste0(round(x, 2), "_", round(y, 2)),
      to_id      = paste0(round(lead(x), 2), "_", round(lead(y), 2)),
      # Sub-segment length
      sub_len    = sqrt((lead(x) - x)^2 + (lead(y) - y)^2)
    ) %>%
    filter(!is.na(to_id), !is.na(sub_len), sub_len > 0) %>%
    ungroup() %>%
    left_join(line_attrs, by = "line_id") %>%
    mutate(
      # Travel time proportional to sub-segment length
      weight = sub_len / speed_mpm,
      mode   = mode_tag
    ) %>%
    filter(
      !is.na(weight),
      !is.infinite(weight),
      !is.nan(weight),
      weight > 0,
      from_id != "NA_NA",
      to_id   != "NA_NA"
    ) %>%
    select(from_id, to_id, weight, mode)
  
  return(edges)
}

# ------ Drive edges (hwy + local roads) ------

get_drive_edges <- function(highway_roads, base_roads) {
  
  # Highway edges with congested travel time
  highway_edges <- lines_to_edges_any_vertex(
    lines        = highway_roads,
    weight_field = "TravelTime",
    mode_tag     = "drive"
  )
  
  # Base road edges for non-highway segments using free flow
  base_edges <- lines_to_edges_any_vertex(
    lines        = base_roads,
    weight_field = "TravelTime",
    mode_tag     = "drive"
  )
  
  bind_rows(highway_edges, base_edges)
}

# Scenario filtering

# Filter transit edges by scenario and year
filter_scenario <- function(edges, scenario) {
  max_year <- SCENARIOS$max_cohort_year[SCENARIOS$scenario == scenario]
  edges |> dplyr::filter(cohort_year <= max_year)
}

# ------ Assemble multimodal graph from component edge lists ------

assemble_graph <- function(road_edges,
                           pedestrian_edges = NULL,
                           transit_edges    = NULL,
                           connector_edges  = NULL) {
  
  components <- list(road_edges, pedestrian_edges,
                     transit_edges, connector_edges)
  components <- Filter(Negate(is.null), components)
  
  graph <- bind_rows(components) %>%
    filter(
      !is.na(from_id), !is.na(to_id), !is.na(weight),
      weight > 0,
      from_id != to_id,        # remove self loops
      from_id != "NA_NA", to_id != "NA_NA",
      !grepl("NA", from_id), !grepl("NA", to_id),
      !is.infinite(weight),
      !is.nan(weight)
    )
  
  return(graph)
}

extract_nodes <- function(graph) {
  bind_rows(
    graph %>% select(id = from_id),
    graph %>% select(id = to_id)
  ) %>%
    distinct(id) %>%
    filter(!is.na(id)) %>%
    mutate(
      x = as.numeric(sub("_.*", "", id)),
      y = as.numeric(sub(".*_", "", id))
    ) %>%
    filter(!is.na(x), !is.na(y))
}

snap_points_to_graph <- function(points, graph, max_dist = 2000,
                                 id_field = "GISJOIN") {
  
  nodes     <- extract_nodes(graph)
  pts_coord <- st_coordinates(points)
  ids       <- points[[id_field]]
  
  node_ids <- sapply(1:nrow(pts_coord), function(i) {
    dists   <- sqrt(
      (nodes$x - pts_coord[i, 1])^2 +
        (nodes$y - pts_coord[i, 2])^2
    )
    nearest <- which.min(dists)
    if (dists[nearest] <= max_dist) {
      nodes$id[nearest]
    } else {
      NA_character_
    }
  })
  
  # Return named vector with GISJOIN as names
  result <- setNames(node_ids, ids)
  
  snap_rate <- mean(!is.na(node_ids)) * 100
  message(sprintf("Snap rate: %.1f%% (%d/%d)",
                  snap_rate,
                  sum(!is.na(node_ids)),
                  length(node_ids)))
  
  return(result)
}

# ------ Calculate OD Travel Time Matrix ------

run_od_matrix <- function(graph, origin_nodes, dest_nodes, label = "") {
  
  valid_o <- origin_nodes[!is.na(origin_nodes)]
  valid_d <- dest_nodes[!is.na(dest_nodes)]
  
  message("Solving: ", label)
  message(sprintf("Origins: %d, Destinations: %d", 
                  length(valid_o), length(valid_d)))
  
  # Build cppRouting graph
  cpp_graph <- makegraph(
    graph %>%
      select(from = from_id, to = to_id, dist = weight) %>%
      filter(!is.na(from), !is.na(to), !is.na(dist), dist > 0),
    directed = TRUE
  )
  
  # Filter to nodes that exist in graph
  graph_nodes <- cpp_graph$dict$ref
  
  valid_o <- valid_o[valid_o %in% graph_nodes]
  valid_d <- valid_d[valid_d %in% graph_nodes]
  
  if (length(valid_o) == 0 || length(valid_d) == 0) {
    warning("No valid origins or destinations found in graph")
    return(data.frame(
      from_GISJOIN   = character(),
      to_GISJOIN     = character(),
      TotalTravelTime = numeric(),
      scenario       = character()
    ))
  }
  
  message(sprintf("Valid origins: %d, Valid destinations: %d",
                  length(valid_o), length(valid_d)))
  
  # Run OD matrix — cppRouting returns a matrix
  od_matrix <- get_distance_matrix(
    Graph = cpp_graph,
    from  = valid_o,
    to    = valid_d,
    allcores = TRUE  # use parallel processing
  )
  
  # Convert matrix to long format
  od_long <- as.data.table(od_matrix, keep.rownames = "from_node") %>%
    melt(
      id.vars      = "from_node",
      variable.name = "to_node",
      value.name   = "TotalTravelTime"
    ) %>%
    mutate(to_node = as.character(to_node)) %>%
    filter(
      !is.na(TotalTravelTime),
      !is.infinite(TotalTravelTime),
      TotalTravelTime > 0
    ) %>%
    mutate(
      from_GISJOIN = names(origin_nodes)[
        match(from_node, origin_nodes)
      ],
      to_GISJOIN = names(dest_nodes)[
        match(to_node, dest_nodes)
      ],
      scenario = label
    ) %>%
    filter(
      !is.na(from_GISJOIN),
      !is.na(to_GISJOIN)
    ) %>%
    select(from_GISJOIN, to_GISJOIN, TotalTravelTime, scenario)
  
  message(sprintf("Solved: %d OD pairs", nrow(od_long)))
  
  return(od_long)
}

# ------ Generate connector edges bridging nodes ------

make_connector_edges_from_nodes <- function(transit_edges,
                                            road_edges,
                                            board_penalty,
                                            alight_penalty,
                                            mode_tag,
                                            max_dist = 500) {
  
  # Get unique transit endpoint node IDs
  transit_nodes <- bind_rows(
    transit_edges %>% select(id = from_id),
    transit_edges %>% select(id = to_id)
  ) %>%
    distinct() %>%
    mutate(
      x = as.numeric(sub("_.*", "", id)),
      y = as.numeric(sub(".*_", "", id))
    )
  
  # Get unique road node IDs
  road_nodes <- bind_rows(
    road_edges %>% select(id = from_id),
    road_edges %>% select(id = to_id)
  ) %>%
    distinct() %>%
    mutate(
      x = as.numeric(sub("_.*", "", id)),
      y = as.numeric(sub(".*_", "", id))
    )
  
  # For each transit node find nearest road node
  nearest <- sapply(1:nrow(transit_nodes), function(i) {
    dists <- sqrt(
      (road_nodes$x - transit_nodes$x[i])^2 +
        (road_nodes$y - transit_nodes$y[i])^2
    )
    nearest_idx <- which.min(dists)
    if (dists[nearest_idx] <= max_dist) {
      road_nodes$id[nearest_idx]
    } else {
      NA_character_
    }
  })
  
  # Build connector edges
  connectors <- data.frame(
    transit_node = transit_nodes$id,
    road_node    = nearest,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(road_node))
  
  cat(sprintf(
    "Built %d connectors from %d transit nodes, %d unconnected\n",
    nrow(connectors),
    nrow(transit_nodes),
    sum(is.na(nearest))
  ))
  
  # Boarding edges: road → transit
  board_edges <- connectors %>%
    transmute(
      from_id = road_node,
      to_id   = transit_node,
      weight  = board_penalty,
      mode    = paste0(mode_tag, "_board")
    )
  
  # Alighting edges: transit → road
  alight_edges <- connectors %>%
    transmute(
      from_id = transit_node,
      to_id   = road_node,
      weight  = alight_penalty,
      mode    = paste0(mode_tag, "_alight")
    )
  
  bind_rows(board_edges, alight_edges)
}

build_transit_edges_from_gtfs <- function(gtfs_dirs,
                                          crs = 26954,
                                          peak_hour_start = 7,
                                          peak_hour_end = 11) {
  
  parse_gtfs_time <- function(t) {
    tryCatch({
      parts <- as.numeric(strsplit(as.character(t), ":")[[1]])
      if (length(parts) != 3 || any(is.na(parts))) return(NA_real_)
      as.numeric(parts[1] * 3600 + parts[2] * 60 + parts[3])
    }, error = function(e) NA_real_)
  }
  
  all_routes     <- list()
  all_trips      <- list()
  all_stop_times <- list()
  all_stops      <- list()
  all_calendar   <- list()
  
  for (dir in gtfs_dirs) {
    source_tag <- basename(dir)
    is_rail    <- grepl("CR|LR", source_tag,
                        ignore.case = TRUE)
    
    cat("Loading:", source_tag, 
        if (is_rail) "(rail)" else "(bus)",
        "\n")
    
    # In your for loop, also load routes.txt
    routes_raw <- read.csv(
      file.path(dir, "routes.txt"),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        route_id   = as.character(route_id),
        source     = source_tag
      ) %>%
      select(route_id, route_type, source)
    
    trips_raw <- read.csv(
      file.path(dir, "trips.txt"),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        trip_id  = as.character(trip_id),
        route_id = as.character(route_id),
        source   = source_tag
      )
    
    stop_times_raw <- read.csv(
      file.path(dir, "stop_times.txt"),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        trip_id  = as.character(trip_id),
        stop_id  = as.character(stop_id),
        source   = source_tag
      )
    
    stops_raw <- read.csv(
      file.path(dir, "stops.txt"),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        stop_id  = as.character(stop_id),
        source   = source_tag
      )
    
    calendar_raw <- read.csv(
      file.path(dir, "calendar.txt"),
      stringsAsFactors = FALSE
    )
    
    all_routes[[source_tag]]     <- routes_raw
    all_trips[[source_tag]]      <- trips_raw
    all_stop_times[[source_tag]] <- stop_times_raw
    all_stops[[source_tag]]      <- stops_raw
    all_calendar[[source_tag]]   <- calendar_raw
  }
  
  # Combine — deduplicate stops across directories by stop_id
  # (not by name, to preserve rail station IDs)
  
  routes <- bind_rows(all_routes) %>%
    distinct(route_id, .keep_all = TRUE) %>%
    mutate(
      mode = case_when(
        route_type == 3 ~ "bus",
        route_type == 2 ~ "commuter_rail",
        route_type %in% c(0, 1) ~ "light_rail",
        TRUE ~ "other"
      )
    )
  
  trips <- bind_rows(all_trips) %>%
    distinct(trip_id, .keep_all = TRUE)
  
  stops <- bind_rows(all_stops) %>%
    distinct(stop_id, .keep_all = TRUE)
  
  stop_times <- bind_rows(all_stop_times) %>%
    mutate(stop_id = as.character(stop_id)) %>%
    distinct(trip_id, stop_sequence, .keep_all = TRUE) %>%
    filter(stop_id %in% stops$stop_id)
  
  calendar <- bind_rows(all_calendar) %>%
    distinct(service_id, .keep_all = TRUE)
  
  cat(sprintf(
    "Combined: %d trips, %d stop_times, %d stops\n",
    nrow(trips), nrow(stop_times), nrow(stops)
  ))
  
  # Weekday services — MT covers M-Th rail, WK and DPSWK cover full week
  weekday_services <- calendar %>%
    filter(monday == 1, tuesday == 1,
           wednesday == 1, thursday == 1) %>%
    pull(service_id)
  
  cat("Weekday service IDs:", paste(weekday_services, collapse = ", "), "\n")
  
  weekday_trips <- trips %>%
    filter(service_id %in% weekday_services)
  
  cat("Weekday trips:", nrow(weekday_trips), "\n")
  
  # Parse stop times
  stop_times_parsed <- stop_times %>%
    filter(trip_id %in% weekday_trips$trip_id) %>%
    mutate(
      departure_sec = vapply(departure_time, parse_gtfs_time, numeric(1)),
      arrival_sec   = vapply(arrival_time,   parse_gtfs_time, numeric(1)),
      departure_hr  = departure_sec / 3600
    )
  
  # First departure per trip
  trip_start_times <- stop_times_parsed %>%
    group_by(trip_id) %>%
    summarise(
      first_departure_hr = min(departure_hr),
      .groups = "drop"
    )
  
  # Peak hour trips
  peak_trips <- weekday_trips %>%
    select(trip_id, route_id) %>%
    left_join(trip_start_times, by = "trip_id") %>%
    filter(
      first_departure_hr >= peak_hour_start,
      first_departure_hr <= peak_hour_end
    )
  
  cat("Peak hour trips:", nrow(peak_trips), "\n")
  
  # One representative trip per route — closest to 8am
  rep_trips <- peak_trips %>%
    left_join(
      stop_times_parsed %>%
        group_by(trip_id) %>%
        summarise(n_stops = n(), .groups = "drop"),
      by = "trip_id"
    ) %>%
    left_join(trip_start_times, by = "trip_id")
  
  cat("Representative trips:", nrow(rep_trips), "\n")
  cat("Routes covered:", paste(rep_trips$route_id, collapse = ", "), "\n")
  
  # Stop times for representative trips
  rep_stop_times <- stop_times_parsed %>%
    filter(trip_id %in% rep_trips$trip_id) %>%
    left_join(
      rep_trips %>% select(trip_id, route_id),
      by = "trip_id"
    ) %>%
    arrange(trip_id, stop_sequence)
  
  # Build segments
  segments <- rep_stop_times %>%
    group_by(trip_id, route_id) %>%
    mutate(
      next_stop_id = lead(stop_id),
      next_arr_sec = lead(arrival_sec),
      segment_min  = (next_arr_sec - departure_sec) / 60
    ) %>%
    filter(
      !is.na(next_stop_id),
      segment_min > 0.5,  # remove cross-street near-zero hops
      segment_min < 60
    ) %>%
    ungroup() %>%
    group_by(route_id, stop_id, next_stop_id) %>%
    summarise(
      segment_min = mean(segment_min, na.rm = TRUE),
      .groups     = "drop"
    ) %>%
    select(route_id, stop_id, next_stop_id, segment_min)
  
  cat("Segments built:", nrow(segments), "\n")
  cat("Routes with segments:", 
      paste(unique(segments$route_id), collapse = ", "), "\n")
  
  # Join stop coordinates
  stop_coords <- stops %>%
    select(stop_id, stop_lat, stop_lon)
  
  seg_coords <- segments %>%
    left_join(stop_coords, by = "stop_id") %>%
    left_join(
      stop_coords %>%
        rename(
          next_stop_id = stop_id,
          next_lat     = stop_lat,
          next_lon     = stop_lon
        ),
      by = "next_stop_id"
    ) %>%
    filter(
      !is.na(stop_lat), !is.na(stop_lon),
      !is.na(next_lat), !is.na(next_lon)
    )
  
  cat("Segments with full coordinates:", nrow(seg_coords), "\n")
  
  # Build linestring geometry
  seg_lines <- lapply(1:nrow(seg_coords), function(i) {
    st_linestring(rbind(
      c(seg_coords$stop_lon[i], seg_coords$stop_lat[i]),
      c(seg_coords$next_lon[i], seg_coords$next_lat[i])
    ))
  })
  
  # Build sf — WGS84 first then reproject
  eh_overlap_stop_ids <- c("34677", "34575", "34574", "34679", "34577", "34578", 
                           "35204", "35205", "35243", "35209", "33902", "35207",
                           "35244", "35211", "35212", "35245")  # RTD stop IDs
  
  segs_sf <- seg_coords %>%
    mutate(geometry = st_sfc(seg_lines, crs = 4326)) %>%
    st_as_sf() %>%
    st_transform(crs) %>%
    rename(TravelTime = segment_min) %>%
    mutate(
      Counterfac = 0,
      source     = "gtfs"
    ) %>%
    left_join(
      routes %>% select(route_id, mode),
      by = "route_id"
    ) %>%
    mutate(
      eh_overlap = route_id %in% c("E", "H") &
        (stop_id %in% eh_overlap_stop_ids |
         stop_id   %in% eh_overlap_stop_ids),
      cohort_year = case_when(
        route_id %in% c("A", "113B")  ~ 2016,
        route_id == "107R"            ~ 2017,
        eh_overlap                    ~ 2017,
        route_id == "113G"            ~ 2019,
        route_id == "117N"            ~ 2020,
        TRUE                          ~ 0
      )
    )
  
  message(sprintf(
    "Built %d segments across %d routes",
    nrow(segs_sf),
    n_distinct(segs_sf$route_id)
  ))
  
  return(segs_sf)
}