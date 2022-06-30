# PACKAGES ####
library(ctmm)
library(data.table)
library(raster)


# DATA ####
setwd("/home/alston92/proj/wrsf")

df <- read.csv("data/serval_example.csv")

l <- as.telemetry(df)

r1 <- raster("data/KZN_2017_Landuse_latlong.tif")


# ANALYSIS ####
SVF <- variogram(l)
GUESS <- ctmm.guess(l, variogram=SVF, interactive=FALSE)
FIT <- ctmm.select(l, GUESS, trace=2)
summary(FIT)

ud_w <- akde(l, FIT, weights=TRUE)
summary(ud_w)

e <- extent(min(l$longitude) - 0.2, max(l$longitude) + 0.2, 
            min(l$latitude) - 0.2, max(l$latitude) + 0.2)
r2 <- crop(r1, e)

m1 <- c(1, 1, 0, 2, 2, 1, 3, 43, 0, 44, 45, 1, 46, 48,0)
rclmat1 <- matrix(m1, ncol = 3, byrow = TRUE)
plantation <- reclassify(r2, rclmat1, include.lowest=TRUE, right = NA)
plantation@data@values <- as.logical(plantation@data@values)

m2 <- c(1, 11, 0, 12, 12, 1,13, 13, 0, 14, 14, 1, 15, 33, 0, 34, 35, 1, 36, 41, 0, 42, 43, 1, 44, 48, 0) 
rclmat2 <- matrix(m2, ncol = 3, byrow = TRUE)
urban <- reclassify(r2, rclmat2, include.lowest=TRUE, right = NA)
urban@data@values <- as.logical(urban@data@values)

# Weighted resource selection function with two habitat covariates
rsf_w <- ctmm:::rsf.fit(l, UD=ud_w, R=list(plantation=plantation,urban=urban), 
                              debias=TRUE, error=0.01)

summary(rsf_w)

# IID case (much slower because N<<n)
ctmm_iid <- ctmm.fit(l,CTMM=ctmm(isotropic=TRUE))
ud_iid <- akde(l,ctmm_iid)

rsf_iid <- ctmm:::rsf.fit(l, UD=ud_iid, R=list(plantation=plantation,urban=urban), 
                          debias=TRUE, error=0.01)
summary(rsf_iid)

# Create variables to populate data.frame
species <- "serval"
iid_pest <- summary(rsf_iid)$CI[1,2]
iid_plcl <- summary(rsf_iid)$CI[1,1]
iid_pucl <- summary(rsf_iid)$CI[1,3]
iid_uest <- summary(rsf_iid)$CI[2,2]
iid_ulcl <- summary(rsf_iid)$CI[2,1]
iid_uucl <- summary(rsf_iid)$CI[2,3]
w_pest <- summary(rsf_w)$CI[1,2]
w_plcl <- summary(rsf_w)$CI[1,1]
w_pucl <- summary(rsf_w)$CI[1,3]
w_uest <- summary(rsf_w)$CI[2,2]
w_ulcl <- summary(rsf_w)$CI[2,1]
w_uucl <- summary(rsf_w)$CI[2,3]
  
# Vector of resultsn
x <- data.frame(species,iid_pest,iid_plcl,iid_pucl,iid_uest,iid_ulcl,iid_uucl,w_pest,w_plcl,w_pucl,w_uest,w_ulcl,w_uucl)
  
# Store results in data.frame
write.table(x, 'results/empirical_results.csv', append=TRUE, row.names=FALSE, col.names=FALSE, sep=',') 

print("Done!")
