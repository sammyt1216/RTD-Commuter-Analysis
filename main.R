library(targets)
library(cppRouting)
library(sf)
library(here)
library(tidyverse)
library(data.table)

tar_make()

coefs <- tar_read(es_coefs_treatment)

att <- tar_read(att_raw)
