
library(readLicorData)


## read all unix licor files and convert into .csv format for upland cotton 
setwd("~/Desktop/Upland_Coker/")

## The following loop should run only once.

files <- list.files()

for (i in 1:length(files)) {
  a <- licorData(location = files[i])
  write.csv(a, paste0(files[i], ".csv"), row.names = FALSE)
}

## read all unix licor files and convert into .csv format for upland cotton 
setwd("~/Desktop/Pima/")

## read all unix licor files and convert into .csv format.
## The following loop should run only once.

files <- list.files()

for (i in 1:length(files)) {
  a <- licorData(location = files[i])
  write.csv(a, paste0(files[i], ".csv"), row.names = FALSE)
}








  
  
