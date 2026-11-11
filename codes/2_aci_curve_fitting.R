library(dplyr)
library(plantecophys)
library(tidyr)
library(ggplot2)
library(tibble)
library(minpack.lm)
library(nlme)
library(car)
library(emmeans)
library(ggpubr)
library(gridExtra)


#########################################################
# reading fit_curve function

source("~/Desktop/Photosynthesis_paper/Codes_revised/fit_curve_function.R")

##############################################################
# PIMA SPECIES
####################

# read tidy aci data 
pima_aci_data <- read.csv("~/Desktop/Photosynthesis_paper/Final_data/pima_final.csv", 
                          row.names = 1)

########################################################
# fit curves one by one

# Plant 1
# Aci curves
p1_16 <- fit_curve(data = pima_aci_data, id_name = "1_16")
plot(p1_16[[1]])
p1_16[[2]]

p1_22 <- fit_curve(data = pima_aci_data, id_name = "1_22")
plot(p1_22[[1]])
p1_22[[2]]

p1_28 <- fit_curve(data = pima_aci_data, id_name = "1_28")
plot(p1_28[[1]])

p1_34 <- fit_curve(data = pima_aci_data, id_name = "1_34")
plot(p1_34[[1]])

p1_40 <- fit_curve(data = pima_aci_data, id_name = "1_40")
plot(p1_40[[1]])

p1_46 <- fit_curve(data = pima_aci_data, id_name = "1_46")
plot(p1_46[[1]])

# 1_46 had data for two curves. So fitting one.
# first 96 rows gives better curves, so first i selected these rows
data_1_46 <- filter(pima_aci_data, id == "1_46")[1:96, ]
# removing points that represents the depression in the curve
data_1_46 <- filter(data_1_46, !(Ci >= 350 & Ci <= 600))

p1_46 <- fit_curve(data = data_1_46, id_name = "1_46", citrans = 350)
plot(p1_46[[1]])
p1_46[[2]]

# extracting the model parameters
coef_1 <- rbind(p1_16[[2]], p1_22[[2]], p1_28[[2]], p1_34[[2]], p1_40[[2]],
                 p1_46[[2]])

## Plant 2
p2_16 <- fit_curve(data = pima_aci_data, id_name =  "2_16")
plot(p2_16[[1]])

p2_22 <- fit_curve(data = pima_aci_data, id_name =  "2_22")
plot(p2_22[[1]])

p2_28 <- fit_curve(data = pima_aci_data, id_name =  "2_28")
plot(p2_28[[1]])

p2_34 <- fit_curve(data = pima_aci_data, id_name =  "2_34")
plot(p2_34[[1]])

p2_40 <- fit_curve(data = pima_aci_data, id_name =  "2_40", citrans = 350)
plot(p2_40[[1]])

# fit_curve(id_name =  "2_46") have no data

coef_2 <- rbind(p2_16[[2]], p2_22[[2]], p2_28[[2]], p2_34[[2]], p2_40[[2]])
# note: has higher vcmax values. P2_re was used as data source.


# Plant 3
p3_16 <- fit_curve(data = pima_aci_data, id_name =  "3_16")
plot(p3_16[[1]])

p3_22 <- fit_curve(data = pima_aci_data, id_name =  "3_22")
plot(p3_22[[1]])
p3_22[[2]]

p3_28 <- fit_curve(data = pima_aci_data, id_name =  "3_28")
plot(p3_28[[1]])

p3_34 <- fit_curve(data = pima_aci_data, id_name =  "3_34")
plot(p3_34[[1]])

p3_40 <- fit_curve(data = pima_aci_data, id_name =  "3_40")
plot(p3_40[[1]])

p3_46 <- fit_curve(data = pima_aci_data, id_name =  "3_46")
plot(p3_46[[1]]) # refitting the curves with some adjustments

data_3_46 <- filter(pima_aci_data, id == "3_46")
data_3_46 <- filter(data_3_46, !(Ci >= 350 & Ci <= 600)) 
# removing points that represents the depression in the curve

p3_46 <- fit_curve(data = data_3_46, id_name = "3_46")
plot(p3_46[[1]])

coef_3 <- rbind(p3_16[[2]], p3_22[[2]], p3_28[[2]], p3_34[[2]], p3_40[[2]], p3_46[[2]])


# Plant 4
p4_16 <- fit_curve(data = pima_aci_data, id_name =  "4_16_re", citrans = 400)
plot(p4_16[[1]])

p4_22 <- fit_curve(data = pima_aci_data, id_name =  "4_22_re")
plot(p4_22[[1]])

p4_28 <- fit_curve(data = pima_aci_data, id_name =  "4_28")
plot(p4_28[[1]])

p4_34 <- fit_curve(data = pima_aci_data, id_name =  "4_34")
plot(p4_34[[1]])

p4_40 <- fit_curve(data = pima_aci_data, id_name =  "4_40")
plot(p4_40[[1]])

p4_46 <- fit_curve(data = pima_aci_data, id_name =  "4_46")
plot(p4_46[[1]])

coef_4 <- rbind(p4_16[[2]], p4_22[[2]], p4_28[[2]], p4_34[[2]], p4_40[[2]], p4_46[[2]])
coef_4[1:2,1] <- c("4_16", "4_22")



# Plant 10
p10_16 <- fit_curve(data = pima_aci_data, id_name =  "10_16")
plot(p10_16[[1]])

p10_22 <- fit_curve(data = pima_aci_data, id_name =  "10_22")
plot(p10_22[[1]])

p10_28 <- fit_curve(data = pima_aci_data, id_name =  "10_28")
plot(p10_28[[1]])

p10_34 <- fit_curve(data = pima_aci_data, id_name =  "10_34")
plot(p10_34[[1]])

p10_40 <- fit_curve(data = pima_aci_data, id_name =  "10_40")
plot(p10_40[[1]])

p10_46 <- fit_curve(data = pima_aci_data, id_name =  "10_46")
plot(p10_46[[1]])

coef_10 <- rbind(p10_16[[2]], p10_22[[2]], p10_28[[2]], p10_34[[2]], p10_40[[2]], p10_46[[2]])

# vcmax values are also low


# Plant 11
p11_16 <- fit_curve(data = pima_aci_data, id_name =  "11_16")
plot(p11_16[[1]])

p11_22 <- fit_curve(data = pima_aci_data, id_name =  "11_22")
plot(p11_22[[1]])

p11_28 <- fit_curve(data = pima_aci_data, id_name =  "11_28")
plot(p11_28[[1]])

p11_34 <- fit_curve(data = pima_aci_data, id_name =  "11_34")
plot(p11_34[[1]])

p11_40 <- fit_curve(data = pima_aci_data, id_name =  "11_40")
plot(p11_40[[1]])

p11_46 <- fit_curve(data = pima_aci_data, id_name =  "11_46")
plot(p11_46[[1]]) 
# refitting this curve with some modifications
data_11_46 <- filter(pima_aci_data, id == "11_46")
data_11_46 <- filter(data_11_46, !(Ci >= 320 & Ci <= 500)) 
# removing points that represents the depression in the curve

p11_46 <- fit_curve(data = data_11_46, id_name = "11_46", citrans = 510)
plot(p11_46[[1]])

coef_11 <- rbind(p11_16[[2]], p11_22[[2]], p11_28[[2]], p11_34[[2]], p11_40[[2]], p11_46[[2]])



# Plant 12

p12_16 <- fit_curve(data = pima_aci_data, id_name =  "12_16", citrans = 450)
plot(p12_16[[1]])

p12_22 <- fit_curve(data = pima_aci_data, id_name =  "12_22")
plot(p12_22[[1]])

p12_28 <- fit_curve(data = pima_aci_data, id_name =  "12_28")
plot(p12_28[[1]])

p12_34 <- fit_curve(data = pima_aci_data, id_name =  "12_34")
plot(p12_34[[1]])

p12_40 <- fit_curve(data = pima_aci_data, id_name =  "12_40")
plot(p12_40[[1]])

p12_46 <- fit_curve(data = pima_aci_data, id_name =  "12_46_re", citrans = 400)
plot(p12_46[[1]])

coef_12 <- rbind(p12_16[[2]], p12_22[[2]], p12_28[[2]], p12_34[[2]], p12_40[[2]], p12_46[[2]])
coef_12[6,1] <- "12_46"


# Plant 17
p17_16 <- fit_curve(data = pima_aci_data, id_name =  "17_16", citrans = 450)
plot(p17_16[[1]])
p17_16[[2]]

p17_22 <- fit_curve(data = pima_aci_data, id_name =  "17_22")
plot(p17_22[[1]])

p17_28 <- fit_curve(data = pima_aci_data, id_name =  "17_28")
plot(p17_28[[1]])

p17_34 <- fit_curve(data = pima_aci_data, id_name =  "17_34")
plot(p17_34[[1]])

p17_40 <- fit_curve(data = pima_aci_data, id_name =  "17_40")
plot(p17_40[[1]])

p17_46 <- fit_curve(data = pima_aci_data, id_name =  "17_46")
plot(p17_46[[1]])
# don't have a good curve

coef_17 <- rbind(p17_16[[2]], p17_22[[2]], p17_28[[2]], p17_34[[2]], p17_40[[2]])

# Vcmax, Jmax too low. Probably will exclude for analysis.


# Plant 18
p18_16 <- fit_curve(data = pima_aci_data, id_name =  "18_16", citrans = 450)
plot(p18_16[[1]])

p18_22 <- fit_curve(data = pima_aci_data, id_name =  "18_22")
plot(p18_22[[1]])

p18_28 <- fit_curve(data = pima_aci_data, id_name =  "18_28", citrans = 400)
plot(p18_28[[1]])

p18_34 <- fit_curve(data = pima_aci_data, id_name =  "18_34", citrans = 300)
plot(p18_34[[1]])

p18_40 <- fit_curve(data = pima_aci_data, id_name =  "18_40", citrans = 400)
plot(p18_40[[1]])

p18_46 <- fit_curve(data = pima_aci_data, id_name =  "18_46")
plot(p18_46[[1]])
# refitting curve with some modification in data

data_18_46 <- filter(pima_aci_data, id == "18_46")
data_18_46 <- filter(data_18_46, !(Ci >= 350 & Ci <= 600))
# removing points that represents the depression in the curve

p18_46 <- fit_curve(data = data_18_46, id_name= "18_46")
plot(p18_46[[1]])

coef_18 <- rbind(p18_16[[2]], p18_22[[2]], p18_28[[2]], p18_34[[2]], p18_40[[2]], p18_46[[2]])


# Plant 19
p19_16 <- fit_curve(data = pima_aci_data, id_name =  "19_16")
plot(p19_16[[1]])

p19_22 <- fit_curve(data = pima_aci_data, id_name =  "19_22")
plot(p19_22[[1]])

p19_28 <- fit_curve(data = pima_aci_data, id_name =  "19_28")
plot(p19_28[[1]])

p19_34 <- fit_curve(data = pima_aci_data, id_name =  "19_34", citrans = 300)
plot(p19_34[[1]])

p19_40 <- fit_curve(data = pima_aci_data, id_name =  "19_40", citrans = 200)
plot(p19_40[[1]])

p19_46 <- fit_curve(data = pima_aci_data, id_name =  "19_46")
plot(p19_46[[1]])
# didn't give a good curve

coef_19 <- rbind(p19_16[[2]], p19_22[[2]], p19_28[[2]], p19_34[[2]], p19_40[[2]])

# Note: very low Vcmax and Jmax


# Plant 20
p20_16 <- fit_curve(data = pima_aci_data, id_name =  "20_16re", citrans = 500)
plot(p20_16[[1]])

p20_22 <- fit_curve(data = pima_aci_data, id_name =  "20_22", citrans = 400)
plot(p20_22[[1]])

p20_28 <- fit_curve(data = pima_aci_data, id_name =  "20_28")
# It has two curves, so selecting the best one
data_20_28 <- filter(pima_aci_data, id == "20_28")
data_20_28 <- data_20_28[97:192,]
p20_28 <- fit_curve(data = data_20_28, id_name =  "20_28")
plot(p20_28[[1]])

p20_34 <- fit_curve(data = pima_aci_data, id_name =  "20_34", citrans = 350)
plot(p20_34[[1]])

p20_40 <- fit_curve(data = pima_aci_data, id_name =  "20_40")
plot(p20_40[[1]])

p20_46 <- fit_curve(data = pima_aci_data, id_name =  "20_46")
plot(p20_46[[1]])
# refitting curve with some adjustments

data_20_46 <- filter(pima_aci_data, id == "20_46")
data_20_46 <- filter(data_20_46, !(Ci >= 340 & Ci <= 500)) 
# removing points that represents the depression in the curve
p20_46 <- fit_curve(data = data_20_46, id_name =  "20_46")
plot(p20_46[[1]])

coef_20 <- rbind(p20_16[[2]], p20_22[[2]], p20_28[[2]], p20_34[[2]], p20_40[[2]], p20_46[[2]])
coef_20[1,1] <- "20_16"


# Plant 25
p25_16 <- fit_curve(data = pima_aci_data, id_name =  "25_16")
plot(p25_16[[1]])

p25_22 <- fit_curve(data = pima_aci_data, id_name =  "25_22")
plot(p25_22[[1]])


p25_28 <- fit_curve(data = pima_aci_data, id_name =  "25_28")
plot(p25_28[[1]])

p25_34 <- fit_curve(data = pima_aci_data, id_name =  "25_34", citrans = 250)
plot(p25_34[[1]])

p25_40 <- fit_curve(data = pima_aci_data, id_name =  "25_40", citrans = 400)
plot(p25_40[[1]])

p25_46 <- fit_curve(data = pima_aci_data, id_name =  "25_46")
plot(p25_46[[1]])
# didn't give a good curve

coef_25 <- rbind(p25_16[[2]], p25_22[[2]], p25_28[[2]], p25_34[[2]], p25_40[[2]])



# Plant 26
p26_16 <- fit_curve(data = pima_aci_data, id_name =  "26_16", citrans = 400)
plot(p26_16[[1]])

p26_22 <- fit_curve(data = pima_aci_data, id_name =  "26_22")
plot(p26_22[[1]])

p26_28 <- fit_curve(data = pima_aci_data, id_name =  "26_28")
plot(p26_28[[1]])

p26_34 <- fit_curve(data = pima_aci_data, id_name =  "26_34", citrans = 300)
plot(p26_34[[1]])

p26_40 <- fit_curve(data = pima_aci_data, id_name =  "26_40", citrans = 350)
plot(p26_40[[1]])

p26_46 <- fit_curve(data = pima_aci_data, id_name =  "26_46")
plot(p26_46[[1]])
# didn't give a good curve

coef_26 <- rbind(p26_16[[2]], p26_22[[2]], p26_28[[2]], p26_34[[2]], p26_40[[2]])


# Plant 27
p27_16 <- fit_curve(data = pima_aci_data, id_name =  "27_16")
plot(p27_16[[1]])

p27_22 <- fit_curve(data = pima_aci_data, id_name =  "27_22")
plot(p27_22[[1]])

p27_28 <- fit_curve(data = pima_aci_data, id_name =  "27_28")
plot(p27_28[[1]])

p27_34 <- fit_curve(data = pima_aci_data, id_name =  "27_34")
plot(p27_34[[1]])

p27_40 <- fit_curve(data = pima_aci_data, id_name =  "27_40", citrans = 300 )
plot(p27_40[[1]])

p27_46 <- fit_curve(data = pima_aci_data, id_name =  "27_46")
plot(p27_46[[1]])
# refiitng the curve with some modifications

data_27_46 <- filter(pima_aci_data, id == "27_46")
data_27_46 <- filter(data_27_46, !(Ci >= 340 & Ci <= 600)) 
# removing points that represents the depression in the curve
p27_46 <- fit_curve(data = data_27_46, id_name =  "27_46", citrans = 400)
plot(p27_46[[1]])

coef_27 <- rbind(p27_16[[2]], p27_22[[2]], p27_28[[2]], p27_34[[2]], p27_40[[2]], p27_46[[2]])



# Plant 28
p28_16 <- fit_curve(data = pima_aci_data, id_name =  "28_16")
plot(p28_16[[1]])

p28_22 <- fit_curve(data = pima_aci_data, id_name =  "28_22")
plot(p28_22[[1]])

p28_28 <- fit_curve(data = pima_aci_data, id_name =  "28_28")
plot(p28_28[[1]])

p28_34 <- fit_curve(data = pima_aci_data, id_name =  "28_34")
plot(p28_34[[1]])

p28_40 <- fit_curve(data = pima_aci_data, id_name =  "28_40")
plot(p28_40[[1]])

p28_46 <- fit_curve(data = pima_aci_data, id_name =  "28_46")
plot(p28_46[[1]])
# refiting the curve with some adjustment

data_28_46 <- filter(pima_aci_data, id == "28_46")
data_28_46 <- filter(data_28_46, !(Ci >= 340 & Ci <= 550)) 
# removing points that represents the depression in the curve
p28_46 <- fit_curve(data = data_28_46, id_name =  "28_46", citrans = 350)
plot(p28_46[[1]])


coef_28 <- rbind(p28_16[[2]], p28_22[[2]], p28_28[[2]], p28_34[[2]], p28_40[[2]], p28_46[[2]])


# plant 33

p33_16 <- fit_curve(data = pima_aci_data, id_name =  "33_16", citrans = 400)
plot(p33_16[[1]])

p33_22 <- fit_curve(data = pima_aci_data, id_name =  "33_22")
plot(p33_22[[1]])

p33_28 <- fit_curve(data = pima_aci_data, id_name =  "33_28")
plot(p33_28[[1]])

p33_34 <- fit_curve(data = pima_aci_data, id_name =  "33_34")
plot(p33_34[[1]])

p33_40 <- fit_curve(data = pima_aci_data, id_name =  "33_40")
plot(p33_40[[1]])

p33_46 <- fit_curve(data = pima_aci_data, id_name =  "33_46")
plot(p33_46[[1]])

coef_33 <- rbind(p33_16[[2]], p33_22[[2]], p33_28[[2]], p33_34[[2]], p33_40[[2]], p33_46[[2]])


# plant 34

p34_16 <- fit_curve(data = pima_aci_data, id_name =  "34_16")
plot(p34_16[[1]])

p34_22 <- fit_curve(data = pima_aci_data, id_name =  "34_22")
plot(p34_22[[1]])

p34_28 <- fit_curve(data = pima_aci_data, id_name =  "34_28")
plot(p34_28[[1]])

p34_34 <- fit_curve(data = pima_aci_data, id_name =  "34_34")
plot(p34_34[[1]])

p34_40 <- fit_curve(data = pima_aci_data, id_name =  "34_40")
plot(p34_40[[1]])

p34_46 <- fit_curve(data = pima_aci_data, id_name =  "34_46", citrans = 500)
plot(p34_46[[1]])

coef_34 <- rbind(p34_16[[2]], p34_22[[2]], p34_28[[2]], p34_34[[2]], p34_40[[2]], p34_46[[2]])


# plant 35

p35_16 <- fit_curve(data = pima_aci_data, id_name =  "35_16", citrans = 400)
plot(p35_16[[1]])

p35_22 <- fit_curve(data = pima_aci_data, id_name =  "35_22")
plot(p35_22[[1]])

p35_28 <- fit_curve(data = pima_aci_data,, id_name =  "35_28")
plot(p35_28[[1]])

p35_34 <- fit_curve(data = pima_aci_data, id_name =  "35_34")
plot(p35_34[[1]])

p35_40 <- fit_curve(data = pima_aci_data, id_name =  "35_40", citrans = 200)
plot(p35_40[[1]])
# refitting with some adjustments

data_35_40 <- filter(pima_aci_data, id == "35_40")
data_35_40 <- filter(data_35_40, !(Ci >= 350 & Ci <= 500)) 
# removing points that represents the depression in the curve
p35_40 <- fit_curve(data = data_35_40, id_name =  "35_40")
plot(p35_40[[1]])


p35_46 <- fit_curve(data = pima_aci_data, id_name =  "35_46")
plot(p35_46[[1]])
# refitting with modifications

data_35_46 <- filter(pima_aci_data, id == "35_46")
data_35_46 <- filter(data_35_46, !(Ci >= 350 & Ci <= 450))
# removing points that represents the depression in the curve
p35_46 <- fit_curve(data = data_35_46, id_name =  "35_46")
plot(p35_46[[1]])

coef_35 <- rbind(p35_16[[2]], p35_22[[2]], p35_28[[2]], p35_34[[2]], p35_40[[2]])


# plant 36

p36_16 <- fit_curve(data = pima_aci_data, id_name =  "36_16", citrans = 400)
plot(p36_16[[1]])

p36_22 <- fit_curve(data = pima_aci_data, id_name =  "36_22")
plot(p36_22[[1]])

p36_28 <- fit_curve(data = pima_aci_data, id_name =  "36_28")
plot(p36_28[[1]])

p36_34 <- fit_curve(data = pima_aci_data, id_name =  "36_34")
plot(p36_34[[1]])

p36_40 <- fit_curve(data = pima_aci_data, id_name =  "36_40", citrans = 350)
plot(p36_40[[1]])

p36_46 <- fit_curve(data = pima_aci_data, id_name =  "36_46")
plot(p36_46[[1]])

coef_36 <- rbind(p36_16[[2]], p36_22[[2]], p36_28[[2]], p36_34[[2]], p36_40[[2]], p36_46[[2]])



# plant 37

p37_16 <- fit_curve(data = pima_aci_data, id_name =  "37_16", citrans = 400)
plot(p37_16[[1]])

p37_22 <- fit_curve(data = pima_aci_data, id_name =  "37_22")
plot(p37_22[[1]])

p37_28 <- fit_curve(data = pima_aci_data, id_name =  "37_28")
plot(p37_28[[1]])

p37_34 <- fit_curve(data = pima_aci_data, id_name =  "37_34")
plot(p37_34[[1]])

p37_40 <- fit_curve(data = pima_aci_data, id_name =  "37_40")
plot(p37_40[[1]])

p37_46 <- fit_curve(data = pima_aci_data, id_name =  "37_46")
plot(p37_46[[1]])
# refitting with some adjustments

data_37_46 <- filter(pima_aci_data, id == "37_46")
data_37_46 <- filter(data_37_46, !(Ci >= 350 & Ci <= 550)) 
# removing points that represents the depression in the curve
p37_46 <- fit_curve(data = data_37_46, id_name =  "37_46")
plot(p37_46[[1]])


coef_37 <- rbind(p37_16[[2]], p37_22[[2]], p37_28[[2]], p37_34[[2]], p37_40[[2]], p37_46[[2]])


# plant 38

p38_16 <- fit_curve(data = pima_aci_data, id_name =  "38_16", citrans = 400)
plot(p38_16[[1]])

p38_22 <- fit_curve(data = pima_aci_data, id_name =  "38_22")
plot(p38_22[[1]])

p38_28 <- fit_curve(data = pima_aci_data, id_name =  "38_28")
plot(p38_28[[1]])

p38_34 <- fit_curve(data = pima_aci_data, id_name =  "38_34")
plot(p38_34[[1]])

p38_40 <- fit_curve(data = pima_aci_data, id_name =  "38_40")
plot(p38_40[[1]])

coef_38 <- rbind(p38_16[[2]], p38_22[[2]], p38_28[[2]], p38_34[[2]], p38_40[[2]])

# Plant 39

p39_16 <- fit_curve(data = pima_aci_data, id_name =  "39_16", citrans = 400)
plot(p39_16[[1]])

p39_22 <- fit_curve(data = pima_aci_data, id_name =  "39_22")
plot(p39_22[[1]])

p39_28 <- fit_curve(data = pima_aci_data, id_name =  "39_28")
plot(p39_28[[1]])

p39_34 <- fit_curve(data = pima_aci_data, id_name =  "39_34")
plot(p39_34[[1]])

p39_40 <- fit_curve(data = pima_aci_data, id_name =  "39_40")
plot(p39_40[[1]])

p39_46 <- fit_curve(data = pima_aci_data, id_name =  "39_46")
plot(p39_46[[1]])

coef_39 <- rbind(p39_16[[2]], p39_22[[2]], p39_28[[2]], p39_34[[2]], p39_40[[2]], p39_46[[2]])


# Plant 40

p40_16 <- fit_curve(data = pima_aci_data, id_name =  "40_16")
plot(p40_16[[1]])

p40_22 <- fit_curve(data = pima_aci_data, id_name =  "40_22")
plot(p40_22[[1]])

p40_28 <- fit_curve(data = pima_aci_data, id_name =  "40_28")
plot(p40_28[[1]])

p40_34 <- fit_curve(data = pima_aci_data, id_name =  "40_34", citrans = 250)
plot(p40_34[[1]])

p40_40 <- fit_curve(data = pima_aci_data, id_name =  "40_40", citrans = 380)
plot(p40_40[[1]])

# p40_46 <- fit_curve(data = pima_aci_data, id_name =  "40_46")
# plot(p40_46[[1]])
# didn't produce a curve

coef_40 <- rbind(p40_16[[2]], p40_22[[2]], p40_28[[2]], p40_34[[2]], p40_40[[2]])


# Plant 41
p41_16 <- fit_curve(data = pima_aci_data, id_name =  "41_16", citrans = 400)
plot(p41_16[[1]])

p41_22 <- fit_curve(data = pima_aci_data, id_name =  "41_22")
plot(p41_22[[1]])

p41_28 <- fit_curve(data = pima_aci_data, id_name =  "41_28")
plot(p41_28[[1]])

p41_34 <- fit_curve(data = pima_aci_data, id_name =  "41_34")
plot(p41_34[[1]])

p41_40 <- fit_curve(data = pima_aci_data, id_name =  "41_40")
plot(p41_40[[1]])

p41_46 <- fit_curve(data = pima_aci_data, id_name =  "41_46", citrans = 450)
plot(p41_46[[1]])


coef_41 <- rbind(p41_16[[2]], p41_22[[2]], p41_28[[2]], p41_34[[2]], p41_40[[2]], p41_46[[2]])

# Plant 42
p42_16 <- fit_curve(data = pima_aci_data, id_name =  "42_16")
plot(p42_16[[1]])

p42_22 <- fit_curve(data = pima_aci_data, id_name =  "42_22")
plot(p42_22[[1]])

p42_28 <- fit_curve(data = pima_aci_data, id_name =  "42_28")
plot(p42_28[[1]])

p42_34 <- fit_curve(data = pima_aci_data, id_name =  "42_34")
plot(p42_34[[1]])

p42_40 <- fit_curve(data = pima_aci_data, id_name =  "42_40")
plot(p42_40[[1]])

p42_46 <- fit_curve(data = pima_aci_data, id_name =  "42_46")
plot(p42_46[[1]])
 # refitting the curve with some modifications

data_42_46 <- filter(pima_aci_data, id == "42_46")
data_42_46 <- filter(data_42_46, !(Ci >= 350 & Ci <= 500))
# removing points that represents the depression in the curve
p42_46 <- fit_curve(data = data_42_46, id_name =  "42_46", citrans = 200)
plot(p42_46[[1]])

coef_42 <- rbind(p42_16[[2]], p42_22[[2]], p42_28[[2]], p42_34[[2]], p42_40[[2]], p42_46[[2]])

# compile all Vcamx and Jmax values in a same dataframe
pima_aci_parameters <- rbind(coef_1, coef_2, coef_3, coef_4, coef_10, coef_11,
                             coef_12, coef_17, coef_18, coef_19, coef_20, coef_25,
                             coef_26, coef_27, coef_28, coef_33, coef_34, coef_35,
                             coef_36, coef_37, coef_38, coef_39, coef_40, coef_41,
                             coef_42)

write.csv(pima_aci_parameters, "~/Desktop/Photosynthesis_paper/aci_parameters_pima.csv",
        row.names = FALSE)


######################
# UPLAND SPECIES
####################

# read tidy aci data 
upland_aci_data <- read.csv("~/Desktop/Photosynthesis_paper/Final_data/upland_final.csv", 
                          row.names = 1)

unique(upland_aci_data$id)

###########################################################
# Plant 1

u1_16 <- fit_curve(data = upland_aci_data, id_name =  "1_16", citrans = 400)
plot(u1_16[[1]])

u1_22 <- fit_curve(data = upland_aci_data, id_name =  "1_22", citrans = 400)
plot(u1_22[[1]])

u1_28 <- fit_curve(data = upland_aci_data, id_name =  "1_28", citrans = 400)
plot(u1_28[[1]])

u1_34 <- fit_curve(data = upland_aci_data, id_name =  "1_34", citrans = 400)
plot(u1_34[[1]])

u1_40 <- fit_curve(data = upland_aci_data, id_name =  "1_40", citrans = 400)
plot(u1_40[[1]])

# u1_46 <- fit_curve(data = upland_aci_data, id_name =  "1_46")
# plot(u1_46[[1]])
# no data for 46 degree

coef_U1 <- rbind(u1_16[[2]], u1_22[[2]], u1_28[[2]], u1_34[[2]], u1_40[[2]])
# didn't give a good data (low values), probably won't use for further analysis.

# Plant 2

u2_16 <- fit_curve(data = upland_aci_data, id_name =  "2_16")
plot(u2_16[[1]])

u2_22 <- fit_curve(data = upland_aci_data, id_name =  "2_22")
plot(u2_22[[1]])

u2_28 <- fit_curve(data = upland_aci_data, id_name =  "2_28")
plot(u2_28[[1]])

u2_34 <- fit_curve(data = upland_aci_data, id_name =  "2_34", citrans = 200)
plot(u2_34[[1]])

u2_40 <- fit_curve(data = upland_aci_data, id_name =  "2_40")
plot(u2_40[[1]])

# u2_46 <- fit_curve(data = upland_aci_data, id_name =  "2_46")
# plot(u2_46[[1]])
# # no a good curve

coef_U2 <- rbind(u2_16[[2]], u2_22[[2]], u2_28[[2]], u2_34[[2]], u2_40[[2]])

# Plant 3

u3_16 <- fit_curve(data = upland_aci_data, id_name =  "3_16", citrans = 300)
plot(u3_16[[1]])

u3_22 <- fit_curve(data = upland_aci_data, id_name =  "3_22")
plot(u3_22[[1]])

u3_28 <- fit_curve(data = upland_aci_data, id_name =  "3_28")
plot(u3_28[[1]])

u3_34 <- fit_curve(data = upland_aci_data, id_name =  "3_34")
plot(u3_34[[1]])

u3_40 <- fit_curve(data = upland_aci_data, id_name =  "3_40")
plot(u3_40[[1]])

u3_46 <- fit_curve(data = upland_aci_data, id_name =  "3_46", citrans = 400)
plot(u3_46[[1]])

coef_U3 <- rbind(u3_16[[2]], u3_22[[2]], u3_28[[2]], u3_34[[2]], u3_40[[2]], u3_46[[2]])

# Plant 4
u4_16 <- fit_curve(data = upland_aci_data, id_name =  "4_16", citrans = 150)
plot(u4_16[[1]])

u4_22 <- fit_curve(data = upland_aci_data, id_name =  "4_22")
plot(u4_22[[1]])

u4_28 <- fit_curve(data = upland_aci_data, id_name =  "4_28")
plot(u4_28[[1]])

u4_34 <- fit_curve(data = upland_aci_data, id_name =  "4_34", citrans = 250)
plot(u4_34[[1]])

u4_40 <- fit_curve(data = upland_aci_data, id_name =  "4_40", citrans = 300)
plot(u4_40[[1]])

u4_46 <- fit_curve(data = upland_aci_data, id_name =  "4_46")
plot(u4_46[[1]])

# didn't give a good curve

coef_U4 <- rbind(u4_16[[2]], u4_22[[2]], u4_28[[2]], u4_34[[2]], u4_40[[2]])

# Plant 5
u5_16 <- fit_curve(data = upland_aci_data, id_name =  "5_16")
plot(u5_16[[1]])

u5_22 <- fit_curve(data = upland_aci_data, id_name =  "5_22")
plot(u5_22[[1]])

u5_28 <- fit_curve(data = upland_aci_data, id_name =  "5_28")
plot(u5_28[[1]])

u5_34 <- fit_curve(data = upland_aci_data, id_name =  "5_34")
plot(u5_34[[1]])

u5_40 <- fit_curve(data = upland_aci_data, id_name =  "5_40", citrans = 400)
plot(u5_40[[1]])

u5_46 <- fit_curve(data = upland_aci_data, id_name =  "5_46", citrans = 400)
plot(u5_46[[1]])


coef_U5 <- rbind(u5_16[[2]], u5_22[[2]], u5_28[[2]], u5_34[[2]], u5_40[[2]], u5_46[[2]])

# Plant 6

u6_16 <- fit_curve(data = upland_aci_data, id_name =  "6_16")
plot(u6_16[[1]])

u6_22 <- fit_curve(data = upland_aci_data, id_name =  "6_22")
plot(u6_22[[1]])

u6_28 <- fit_curve(data = upland_aci_data, id_name =  "6_28")
plot(u6_28[[1]])

u6_34 <- fit_curve(data = upland_aci_data, id_name =  "6_34", citrans = 200)
plot(u6_34[[1]])

u6_40 <- fit_curve(data = upland_aci_data, id_name =  "6_40", citrans = 250)
plot(u6_40[[1]])

# u6_46 <- fit_curve(data = upland_aci_data, id_name =  "6_46")
# plot(u6_46[[1]])
# # have no data

coef_U6 <- rbind(u6_16[[2]], u6_22[[2]], u6_28[[2]], u6_34[[2]], u6_40[[2]])


# Plant 7

u7_16 <- fit_curve(data = upland_aci_data, id_name =  "7_16")
plot(u7_16[[1]])

u7_22 <- fit_curve(data = upland_aci_data, id_name =  "7_22")
plot(u7_22[[1]])

u7_28 <- fit_curve(data = upland_aci_data, id_name =  "7_28")
plot(u7_28[[1]])

u7_34 <- fit_curve(data = upland_aci_data, id_name =  "7_34")
plot(u7_34[[1]])

u7_40 <- fit_curve(data = upland_aci_data, id_name =  "7_40", citrans = 250)
plot(u7_40[[1]])

u7_46 <- fit_curve(data = upland_aci_data, id_name =  "7_46")
plot(u7_46[[1]])
# refitting the curve with some modifications
data_u7_46 <- filter(upland_aci_data, id == "7_46")
data_u7_46 <- filter(data_u7_46, !(Ci >= 340 & Ci <= 450))
u7_46 <- fit_curve(data = data_u7_46, id_name =  "7_46")
plot(u7_46[[1]])

coef_U7 <- rbind(u7_16[[2]], u7_22[[2]], u7_28[[2]], u7_34[[2]], u7_40[[2]], u7_46[[2]])

# Plant 8

u8_16 <- fit_curve(data = upland_aci_data, id_name =  "8_16")
plot(u8_16[[1]])

u8_22 <- fit_curve(data = upland_aci_data, id_name =  "8_22")
plot(u8_22[[1]])

u8_28 <- fit_curve(data = upland_aci_data, id_name =  "8_28")
plot(u8_28[[1]])

u8_34 <- fit_curve(data = upland_aci_data, id_name =  "8_34")
plot(u8_34[[1]])

u8_40 <- fit_curve(data = upland_aci_data, id_name =  "8_40")
plot(u8_40[[1]])

u8_46 <- fit_curve(data = upland_aci_data, id_name =  "8_46")
plot(u8_46[[1]])

# refitting the curve with some modifications
data_u8_46 <- filter(upland_aci_data, id == "8_46")
data_u8_46 <- filter(data_u8_46, !(Ci >= 350 & Ci <= 550))
u8_46 <- fit_curve(data = data_u8_46, id_name =  "8_46", citrans = 400)
plot(u8_46[[1]])

coef_U8 <- rbind(u8_16[[2]], u8_22[[2]], u8_28[[2]], u8_34[[2]], u8_40[[2]], u8_46[[2]])

# Plant 9

u9_16 <- fit_curve(data = upland_aci_data, id_name =  "9_16", citrans = 400)
plot(u9_16[[1]])

u9_22 <- fit_curve(data = upland_aci_data, id_name =  "9_22")
plot(u9_22[[1]])

u9_28 <- fit_curve(data = upland_aci_data, id_name =  "9_28")
plot(u9_28[[1]])

u9_34 <- fit_curve(data = upland_aci_data, id_name =  "9_34")
plot(u9_34[[1]])

u9_40 <- fit_curve(data = upland_aci_data, id_name =  "9_40")
plot(u9_40[[1]])

u9_46 <- fit_curve(data = upland_aci_data, id_name =  "9_46")
plot(u9_46[[1]])

coef_U9 <- rbind(u9_16[[2]], u9_22[[2]], u9_28[[2]], u9_34[[2]], u9_40[[2]], u9_46[[2]])

# Plant 10

u10_16 <- fit_curve(data = upland_aci_data, id_name =  "10_16")
plot(u10_16[[1]])

u10_22 <- fit_curve(data = upland_aci_data, id_name =  "10_22")
plot(u10_22[[1]])

u10_28 <- fit_curve(data = upland_aci_data, id_name =  "10_28")
plot(u10_28[[1]])

u10_34 <- fit_curve(data = upland_aci_data, id_name =  "10_34")
plot(u10_34[[1]])

u10_40 <- fit_curve(data = upland_aci_data, id_name =  "10_40")
plot(u10_40[[1]])

u10_46 <- fit_curve(data = upland_aci_data, id_name =  "10_46", citrans = 320)
plot(u10_46[[1]])
u10_46[[2]]

coef_U10 <- rbind(u10_16[[2]], u10_22[[2]], u10_28[[2]], u10_34[[2]], u10_40[[2]], u10_46[[2]])
write.csv(coef_U10, "coef_U10.csv", row.names = FALSE)


# Plant 11

u11_16 <- fit_curve(data = upland_aci_data, id_name =  "11_16")
plot(u11_16[[1]])

u11_22 <- fit_curve(data = upland_aci_data, id_name =  "11_22")
plot(u11_22[[1]])

u11_28 <- fit_curve(data = upland_aci_data, id_name =  "11_28")
plot(u11_28[[1]])

u11_34 <- fit_curve(data = upland_aci_data, id_name =  "11_34")
plot(u11_34[[1]])

u11_40 <- fit_curve(data = upland_aci_data, id_name =  "11_40", citrans = 400)
plot(u11_40[[1]])

u11_46 <- fit_curve(data = upland_aci_data, id_name =  "11_46")
plot(u11_46[[1]])
# not a good curve 

coef_U11 <- rbind(u11_16[[2]], u11_22[[2]], u11_28[[2]], u11_34[[2]], u11_40[[2]])

# Plant 12


u12_16 <- fit_curve(data = upland_aci_data, id_name =  "12_16", citrans = 300)
plot(u12_16[[1]])

u12_22 <- fit_curve(data = upland_aci_data, id_name =  "12_22", citrans = 300)
plot(u12_22[[1]])

u12_28 <- fit_curve(data = upland_aci_data, id_name =  "12_28", citrans = 300)
plot(u12_28[[1]])

u12_34 <- fit_curve(data = upland_aci_data, id_name =  "12_34", citrans = 200)
plot(u12_34[[1]])

u12_40 <- fit_curve(data = upland_aci_data, id_name =  "12_40")
plot(u12_40[[1]])

u12_46 <- fit_curve(data = upland_aci_data, id_name =  "12_46")
plot(u12_46[[1]])
# not a good curve 

coef_U12 <- rbind(u12_16[[2]], u12_22[[2]], u12_28[[2]], u12_34[[2]], u12_40[[2]])


# Plant 13


u13_16 <- fit_curve(data = upland_aci_data, id_name =  "13_16")
plot(u13_16[[1]])

u13_22 <- fit_curve(data = upland_aci_data, id_name =  "13_22")
plot(u13_22[[1]])

u13_28 <- fit_curve(data = upland_aci_data, id_name =  "13_28")
plot(u13_28[[1]])

u13_34 <- fit_curve(data = upland_aci_data, id_name =  "13_34")
plot(u13_34[[1]])

u13_40 <- fit_curve(data = upland_aci_data, id_name =  "13_40", citrans = 450)
plot(u13_40[[1]])

u13_46 <- fit_curve(data = upland_aci_data, id_name =  "13_46")
plot(u13_46[[1]])
# not a good curve 

coef_U13 <- rbind(u13_16[[2]], u13_22[[2]], u13_28[[2]], u13_34[[2]], u13_40[[2]])

# Plant 14

u14_16 <- fit_curve(data = upland_aci_data, id_name =  "14_16", citrans = 450)
plot(u14_16[[1]])


u14_22 <- fit_curve(data = upland_aci_data, id_name =  "14_22", citrans = 450)
plot(u14_22[[1]])

u14_28 <- fit_curve(data = upland_aci_data, id_name =  "14_28", citrans = 150)
plot(u14_28[[1]])
u14_28[[2]]

u14_34 <- fit_curve(data = upland_aci_data, id_name =  "14_34", citrans = 130)
plot(u14_34[[1]])
u14_34[[2]]

u14_40 <- fit_curve(data = upland_aci_data, id_name =  "14_40", citrans = 150)
plot(u14_40[[1]])

u14_46 <- fit_curve(data = upland_aci_data, id_name =  "14_46", citrans = 200)
plot(u14_46[[1]])
u14_46[[2]]

# not a good curve

coef_U14 <- rbind(u14_16[[2]], u14_22[[2]], u14_28[[2]], u14_34[[2]], u14_40[[2]])



# Plant 43 

u43_16 <- fit_curve(data = upland_aci_data, id_name =  "43_16", citrans = 400)
plot(u43_16[[1]])

u43_22 <- fit_curve(data = upland_aci_data, id_name =  "43_22")
plot(u43_22[[1]])

u43_28 <- fit_curve(data = upland_aci_data, id_name =  "43_28")
plot(u43_28[[1]])

u43_34 <- fit_curve(data = upland_aci_data, id_name =  "43_34")
plot(u43_34[[1]])

u43_40 <- fit_curve(data = upland_aci_data, id_name =  "43_40", citrans = 400)
plot(u43_40[[1]])

u43_46 <- fit_curve(data = upland_aci_data, id_name =  "43_46")
plot(u43_46[[1]])
# not a good curve

coef_U43 <- rbind(u43_16[[2]], u43_22[[2]], u43_28[[2]], u43_34[[2]], u43_40[[2]])

# Plant 46

u46_16 <- fit_curve(data = upland_aci_data, id_name =  "46_16")
plot(u46_16[[1]])

u46_22 <- fit_curve(data = upland_aci_data, id_name =  "46_22")
plot(u46_22[[1]])

u46_28 <- fit_curve(data = upland_aci_data, id_name =  "46_28")
plot(u46_28[[1]])

u46_34 <- fit_curve(data = upland_aci_data, id_name =  "46_34")
plot(u46_34[[1]])

u46_40 <- fit_curve(data = upland_aci_data, id_name =  "46_40")
plot(u46_40[[1]])

u46_46 <- fit_curve(data = upland_aci_data, id_name =  "46_46")
plot(u46_46[[1]])

# refitting the curve with some modifications
data_u46_46 <- filter(upland_aci_data, id == "46_46")
data_u46_46 <- filter(data_u46_46, !(Ci >= 330 & Ci <= 400))
u46_46 <- fit_curve(data = data_u46_46, id_name =  "46_46")
plot(u46_46[[1]])

coef_U46 <- rbind(u46_16[[2]], u46_22[[2]], u46_28[[2]], u46_34[[2]], u46_40[[2]], u46_46[[2]] )

# Plant 47
u47_16 <- fit_curve(data = upland_aci_data, id_name =  "47_16")
plot(u47_16[[1]])

u47_22 <- fit_curve(data = upland_aci_data, id_name =  "47_22", citrans = 200)
plot(u47_22[[1]])

u47_28 <- fit_curve(data = upland_aci_data, id_name =  "47_28re", citrans = 400)
plot(u47_28[[1]])

u47_34 <- fit_curve(data = upland_aci_data, id_name =  "47_34", citrans = 400)
plot(u47_34[[1]])

u47_40 <- fit_curve(data = upland_aci_data, id_name =  "47_40", citrans = 320)
plot(u47_40[[1]])
u47_40[[2]]

u47_46 <- fit_curve(data = upland_aci_data, id_name =  "47_46", citrans = 400)
plot(u47_46[[1]])

coef_U47 <- rbind(u47_16[[2]], u47_22[[2]], u47_28[[2]], u47_34[[2]], u47_40[[2]], u47_46[[2]] )
coef_U47[3,1] <- "47_28"

# Plant 48

u48_16 <- fit_curve(data = upland_aci_data, id_name =  "48_16", citrans = 400)
plot(u48_16[[1]])

u48_22 <- fit_curve(data = upland_aci_data, id_name =  "48_22")
plot(u48_22[[1]])

u48_28 <- fit_curve(data = upland_aci_data, id_name =  "48_28")
plot(u48_28[[1]])

u48_34 <- fit_curve(data = upland_aci_data, id_name =  "48_34")
plot(u48_34[[1]])

u48_40 <- fit_curve(data = upland_aci_data, id_name =  "48_40")
plot(u48_40[[1]])

u48_46 <- fit_curve(data = upland_aci_data, id_name =  "48_46")
plot(u48_46[[1]])

# refitting the curve with some modifications
data_u48_46 <- filter(upland_aci_data, id == "48_46")
data_u48_46 <- filter(data_u48_46, !(Ci >= 330 & Ci <= 450))
u48_46 <- fit_curve(data = data_u48_46, id_name =  "48_46")
plot(u48_46[[1]])

coef_U48 <- rbind(u48_16[[2]], u48_22[[2]], u48_28[[2]], u48_34[[2]], u48_40[[2]], u48_46[[2]])

# Plant 50

u50_16 <- fit_curve(data = upland_aci_data, id_name =  "50_16")
plot(u50_16[[1]])

u50_22 <- fit_curve(data = upland_aci_data, id_name =  "50_22")
plot(u50_22[[1]])

u50_28 <- fit_curve(data = upland_aci_data, id_name =  "50_28")
plot(u50_28[[1]])

u50_34 <- fit_curve(data = upland_aci_data, id_name =  "50_34")
plot(u50_34[[1]])

u50_40 <- fit_curve(data = upland_aci_data, id_name =  "50_40")
plot(u50_40[[1]])

u50_46 <- fit_curve(data = upland_aci_data, id_name =  "50_46")
plot(u50_46[[1]])


coef_U50 <- rbind(u50_16[[2]], u50_22[[2]], u50_28[[2]], u50_34[[2]], u50_40[[2]], u50_46[[2]] )

# Plant 51
u51_16 <- fit_curve(data = upland_aci_data, id_name =  "51_16")
plot(u51_16[[1]])

u51_22 <- fit_curve(data = upland_aci_data, id_name =  "51_22")
plot(u51_22[[1]])

u51_28 <- fit_curve(data = upland_aci_data, id_name =  "51_28", citrans = 200)
plot(u51_28[[1]])

u51_34 <- fit_curve(data = upland_aci_data, id_name =  "51_34", citrans = 250)
plot(u51_34[[1]])

u51_40 <- fit_curve(data = upland_aci_data, id_name =  "51_40", citrans = 400)
plot(u51_40[[1]])

u51_46 <- fit_curve(data = upland_aci_data, id_name =  "51_46")
plot(u51_46[[1]])
# not a good curve, negative Anet values


coef_U51 <- rbind(u51_16[[2]], u51_22[[2]], u51_28[[2]], u51_34[[2]], u51_40[[2]] )


# Plant 52

u52_16 <- fit_curve(data = upland_aci_data, id_name =  "52_16re")
plot(u52_16[[1]])

u52_22 <- fit_curve(data = upland_aci_data, id_name =  "52_22")
plot(u52_22[[1]])

u52_28 <- fit_curve(data = upland_aci_data, id_name =  "52_28")
plot(u52_28[[1]])

u52_34 <- fit_curve(data = upland_aci_data, id_name =  "52_34")
plot(u52_34[[1]])

u52_40 <- fit_curve(data = upland_aci_data, id_name =  "52_40")
plot(u52_40[[1]])

u52_46 <- fit_curve(data = upland_aci_data, id_name =  "52_46")
plot(u52_46[[1]])
# refitting the curve with some modifications
data_u52_46 <- filter(upland_aci_data, id == "52_46")
data_u52_46 <- filter(data_u52_46, !(Ci >= 330 & Ci <= 500))
u52_46 <- fit_curve(data = data_u52_46, id_name =  "52_46")
plot(u52_46[[1]])

coef_U52 <- rbind(u52_16[[2]], u52_22[[2]], u52_28[[2]], u52_34[[2]], u52_40[[2]], u52_46[[2]])
coef_U52[[1,1]] <- "52_16"

# Plant 53
u53_16 <- fit_curve(data = upland_aci_data, id_name =  "53_16")
plot(u53_16[[1]])

u53_22 <- fit_curve(data = upland_aci_data, id_name =  "53_22")
plot(u53_22[[1]])

u53_28 <- fit_curve(data = upland_aci_data, id_name =  "53_28")
plot(u53_28[[1]])

u53_34 <- fit_curve(data = upland_aci_data, id_name =  "53_34")
plot(u53_34[[1]])

u53_40 <- fit_curve(data = upland_aci_data, id_name =  "53_40")
plot(u53_40[[1]])

u53_46 <- fit_curve(data = upland_aci_data, id_name =  "53_46", citrans = 300)
plot(u53_46[[1]])

coef_U53 <- rbind(u53_16[[2]], u53_22[[2]], u53_28[[2]], u53_34[[2]], u53_40[[2]], u53_46[[2]])

# Plant 54
u54_16 <- fit_curve(data = upland_aci_data, id_name =  "54_16")
plot(u54_16[[1]])

u54_22 <- fit_curve(data = upland_aci_data, id_name =  "54_22")
plot(u54_22[[1]])

u54_28 <- fit_curve(data = upland_aci_data, id_name =  "54_28")
plot(u54_28[[1]])

u54_34 <- fit_curve(data = upland_aci_data, id_name =  "54_34", citrans = 400)
plot(u54_34[[1]])

u54_40 <- fit_curve(data = upland_aci_data, id_name =  "54_40")
plot(u54_40[[1]])

u54_46 <- fit_curve(data = upland_aci_data, id_name =  "54_46", citrans = 400)
plot(u54_46[[1]])

coef_U54 <- rbind(u54_16[[2]], u54_22[[2]], u54_28[[2]], u54_34[[2]], u54_40[[2]], u54_46[[2]])

# compile all Vcamx and Jmax values in a same dataframe
upland_aci_parameters <- rbind(coef_U1, coef_U2, coef_U3, coef_U4, coef_U5, coef_U6,
                             coef_U7, coef_U8, coef_U9, coef_U10, coef_U11, coef_U12,
                             coef_U13, coef_U14, coef_U43, coef_U46, coef_U47, coef_U48,
                             coef_U50, coef_U51, coef_U52, coef_U53, coef_U54)


write.csv(upland_aci_parameters, "~/Desktop/Photosynthesis_paper/aci_parameters_upland.csv",
       row.names = FALSE)

