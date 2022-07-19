# PACKAGES ####
library(data.table)
library(raster)
library(ctmm)

# ANALYSIS ####
# Read job number from command line
sim_no <- 1
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
  
# Extract variables of interest ###
sim_no <- sim_no
est <- summary(rsf)$CI[1,2]
lcl <- summary(rsf)$CI[1,1]
ucl <- summary(rsf)$CI[1,3]
  
#################################
# Vector of results to return
x <- data.frame(sim_no, est, lcl, ucl)
  
# Store results in data.frame
write.table(x, 'results/wrsf_sim_results_lo_wrsf_selection_subsample_w.csv', append=TRUE, row.names=FALSE, col.names=FALSE, sep=',') 
  
# Print indicators of progress
print("Done!")

