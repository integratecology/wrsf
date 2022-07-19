# PACKAGES ####
library(raster)
library(ctmm)

# ANALYSIS ####
# Read job number from command line
args = commandArgs(trailingOnly=TRUE)
sim_no <- as.numeric(args[1])
print(sim_no)

# Create raster
r1 <- raster(nrows = 1000, ncols = 1000,
             xmn = -0.05, xmx = 0.05, ymn = -0.05, ymx = 0.05,
             vals = as.factor(rep(c("A","B"), 500000)))
projection(r1) <- "+proj=longlat +datum=WGS84 +nodefs"

df <- fread("/home/alston92/proj/wrsf/data/subsample_sim.csv")

sim_sub <- as.telemetry(df)

# Fit the movement model to the simulated data
svf <- variogram(sim_sub)
guess <- ctmm.guess(sim_sub, variogram=svf, interactive=FALSE)
fit <- ctmm.select(sim_sub, guess, trace=2) #
summary(fit)
  
# Calculate the UDs ###
ud <- akde(sim_sub, fit, weights=TRUE)
summary(ud)
  
# Fit the RSFs ###
set.seed(sim_no)
rsf <- ctmm:::rsf.fit(sim_sub, UD=ud, R=list(test=r1), debias=TRUE, error=0.01, integrator="Riemann", interpolate=FALSE)
summary(rsf)
eTime <- Sys.time()
  
# Extract variables of interest ###
sim_no <- sim_no
samp_dur <- nd
wrsf_coef <- summary(rsf)$CI[1,2]
wrsf_lcl <- summary(rsf)$CI[1,1]
wrsf_ucl <- summary(rsf)$CI[1,3]
runtime <- difftime(eTime, sTime, units="mins")
area <- summary(rsf)$CI[2,2]
area_lcl <- summary(rsf)$CI[2,1]
area_ucl <- summary(rsf)$CI[2,3]
true_area <- trueRngArea
habitat1 <- sum(sim_sub$habitat[sim_sub$habitat==1])
habitat2 <- sum(sim_sub$habitat[sim_sub$habitat==2])/2
  
  #################################
  # Vector of results to return
  x <- data.frame(sim_no, samp_dur, wrsf_coef, wrsf_lcl, wrsf_ucl, area, area_lcl, area_ucl, true_area, habitat1, habitat2, runtime)
  
  # Store results in data.frame
  write.table(x, 'results/wrsf_sim_results_lo_wrsf_selection_duration.csv', append=TRUE, row.names=FALSE, col.names=FALSE, sep=',') 
  
  # Print indicators of progress
  print(nd)
  print(eTime)
  print(eTime - sTime)
  
}
