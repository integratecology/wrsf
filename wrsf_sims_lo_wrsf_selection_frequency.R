# PACKAGES ####
library(raster)
library(ctmm)

# ANALYSIS ####
# Read job number from command line
args = commandArgs(trailingOnly=TRUE)
sim_no <- args[1]
print(sim_no)

# Define the OUF model parameters ###
# Length of a day in seconds
ds <- 86400

# Spatial variance in m^2
sig <- 200000

# True 95% range area
trueRngArea <- -2*log(0.05)*pi*sig

# Specify an OUF model for simulation
mod <- ctmm(tau=c(ds,ds-1), isotropic=TRUE, sigma=sig, mu=c(0,0))

# Simulation with varying sampling interval ####
# Sampling frequencies to quantify
samp <- c(2, 4, 8, 16, 32)

# Create raster
r1 <- raster(nrows = 1000, ncols = 1000,
             xmn = -0.05, xmx = 0.05, ymn = -0.05, ymx = 0.05,
             vals = as.factor(rep(c("A","B"), 500000)))
projection(r1) <- "+proj=longlat +datum=WGS84 +nodefs"

# Record start time to monitor how long replicates take to compute
sTime <- Sys.time()
print(sTime)

# Loop over sampling frequencies (samp)
for(i in 1:length(samp)){
  
  # Specify variables to manipulate sampling frequency while holding duration constant
  set.seed(sim_no) # unique seed for each simulation run
  nd <- 64 # number of days
  pd <- samp[i] # number of sampled points per day
  
  # Sampling schedule
  st <- 1:(nd*pd)*(ds/pd) 
    
  # Simulate from the movement model ###
  set.seed(sim_no)
  sim <- simulate(mod, t=st, complete=TRUE)

  # Create a new track with habitat selection (shift half of points in movement track right by 0.0001 degrees)
  df <- data.frame(sim)
  df <- df[,6:8]
  pts <- df[,1:2]
  sp <- sp::SpatialPoints(pts, proj4string = sp::CRS(projection(r1)))
  df$habitat <- raster::extract(r1, sp::spTransform(sp, sp::CRS(projection(r1))))
  df$count <- ave(df$habitat==df$habitat, df$habitat, FUN=cumsum)
  df2 <- df[df$habitat==2 & df$count %% 2,]
  df$longitude <- ifelse(row.names(df) %in% row.names(df2)==TRUE, df$longitude + 0.0001, df$longitude)
  
  pts2 <- df[,1:2]
  sp <- sp::SpatialPoints(pts2, proj4string = sp::CRS(projection(r1)))
  df$habitat <- raster::extract(r1, sp::spTransform(sp, sp::CRS(projection(r1))))
  df$count <- ave(df$habitat==df$habitat, df$habitat, FUN=cumsum)
  head(df, 20)
  print(paste0("There are ", max(df$count[df$habitat==1]), " locations in Habitat A and ", max(df$count[df$habitat==2]), " locations in Habitat B"))
  df$individual.local.identifier <- sim_no
  
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
  samp_freq <- pd
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
  x <- data.frame(sim_no, samp_freq, wrsf_coef, wrsf_lcl, wrsf_ucl, area, area_lcl, area_ucl, true_area, habitat1, habitat2, runtime)
  
  # Store results in data.frame
  write.table(x, 'results/wrsf_sim_results_lo_wrsf_selection_frequency.csv', append=TRUE, row.names=FALSE, col.names=FALSE, sep=',') 
  
  # Print indicators of progress
  print(nd)
  print(eTime)
  print(eTime - sTime)
  
}
