
library(dplyr)
library(stringr)
library(reshape)
library(plantecophys)
library(ggplot2)

# PIMA species

# create the file list in the directory
setwd("~/Desktop/PhD Pojects/dissertation/Chapter 3/ACi_data/Pima/Pima_csv_files/")

files_pima <- list.files(pattern = ".*\\.csv$")

files_pima <- setNames(files_pima, str_extract(basename(files_pima), 
                                     '.*(?=\\.csv)'))
# read the data and merge into a single file
merged_data_pima <- lapply(files_pima, read.csv) %>%
  merge_all() %>%
  arrange(id, elapsed)

# check unique id
unique(merged_data_pima$id)

# create a column for plants and convert into integers
merged_data_pima$plants <-  gsub("^(\\d+)_.*", "\\1", merged_data_pima$id)
merged_data_pima$plants <- as.integer(merged_data_pima$plants)

# create a column for temp setpoint
merged_data_pima$temp_setpoint <- round(merged_data_pima$Tleaf, digits = 0)
unique(merged_data_pima$temp_setpoint)
# it has 17 and 45, which was not actual set point
# adjust temp setpoint for values different from the actual setpoint

merged_data_pima$temp_setpoint <- 
  ifelse(merged_data_pima$temp_setpoint == 17, 16, merged_data_pima$temp_setpoint)
merged_data_pima$temp_setpoint <- 
  ifelse(merged_data_pima$temp_setpoint == 45, 46, merged_data_pima$temp_setpoint)
unique(merged_data_pima$temp_setpoint)

# only keep columns required to fit aci curve
pima_final <- merged_data_pima %>%
  select(plants, id, A, Ci, Ca, gsw, 
         CO2_s,	CO2_r,	H2O_s,	H2O_r,
         Qin, VPDleaf, Flow,	Tair, Tleaf, Pa, temp_setpoint) %>%
  arrange(plants)

# create unique plants id (add P for Pima, 1 changes to P1)
pima_final <- pima_final |> 
  mutate(plants = paste0("P" ,plants))


write.csv(pima_final, "~/Desktop/Photosynthesis_paper/pima_final.csv")



# Upland species

# create the file list in the directory
setwd("~/Desktop/PhD Pojects/dissertation/Chapter 3/ACi_data/Upland_Coker/upland_csv_files/")
files_upland <- list.files(pattern = ".*\\.csv$")

files_upland <- setNames(files_upland, str_extract(basename(files_upland), 
                                               '.*(?=\\.csv)'))
# read the data and merge into a single file
merged_data_upland <- lapply(files_upland, read.csv) %>%
  merge_all() %>%
  arrange(id, elapsed)

# check unique id
unique(merged_data_upland$id)

# create a colum for plants and convert into integers
merged_data_upland$plants <-  gsub("^(\\d+)_.*", "\\1", merged_data_upland$id)
merged_data_upland$plants <- as.integer(merged_data_upland$plants)

# create a column for temp setpoint
merged_data_upland$temp_setpoint <- round(merged_data_upland$Tleaf, digits = 0)
unique(merged_data_upland$temp_setpoint)
# it has 17, 44,  45, which was not actual set point
# adjust temp setpoint for values differnt from the actual setpoint

merged_data_upland$temp_setpoint <- 
  ifelse(merged_data_upland$temp_setpoint == 17, 16, merged_data_upland$temp_setpoint)
merged_data_upland$temp_setpoint <- 
  ifelse(merged_data_upland$temp_setpoint %in% c(44, 45), 46, merged_data_upland$temp_setpoint)
unique(merged_data_upland$temp_setpoint)

# only keep columns required to fit aci curve
upland_final <- merged_data_upland %>%
  select(plants, id, A, Ci, Ca, gsw, 
         CO2_s,	CO2_r,	H2O_s,	H2O_r,
         Qin, VPDleaf, Flow,	Tair, Pa, Tleaf, temp_setpoint) %>%
  arrange(plants)

# create unique plants id (add U for upland, 1 changes to U1)
upland_final <- upland_final|> 
  mutate(plants = paste0("U" ,plants))

write.csv(upland_final, "~/Desktop/Photosynthesis_paper/upland_final.csv")




