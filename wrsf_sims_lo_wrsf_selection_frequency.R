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
samp <- c(1, 4, 16, 64, 256)

# Create raster
r1 <- raster(nrows = 1000, ncols = 1000, 
             xmn = -5000, xmx = 5000, ymn = -5000, ymx = 5000,
             vals = as.factor(rep(c("A","B"), 500000)))
projection(r1) <- "+proj=aeqd +lon_0=0 +lat_0=0 +datum=WGS84"

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
  sim <- simulate(mod, t=st, complete = TRUE)

  # Create a new track with habitat selection (shift half of points in movement track right by 0.0001 degrees)
  sim_sub <- sim
  
  df <- data.frame(sim_sub)
  pts <- df[,2:3]
  sp <- sp::SpatialPoints(pts, proj4string = sp::CRS(projection(r1)))
  df$habitat <- raster::extract(r1, sp::spTransform(sp, sp::CRS(projection(r1))))
  df$count <- ave(df$habitat==df$habitat, df$habitat, FUN=cumsum)
  df <- df[df$habitat==2 & df$count %% 2,]
  sim_sub$longitude <- ifelse(row.names(sim) %in% row.names(df)==TRUE, sim_sub$longitude + 0.00009009009, sim_sub$longitude)

  # Check number of points in each habitat
  pts2 <- sim_sub[,6:7]
  sp <- sp::SpatialPoints(pts2, proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))
  sim_sub$habitat <- raster::extract(r1, sp::spTransform(sp, sp::CRS(projection(r1))))
  sim_sub$count <- ave(sim_sub$habitat==sim_sub$habitat, sim_sub$habitat, FUN=cumsum)
  head(sim_sub, 20)
  print(paste0("There are ", max(sim_sub$count[sim_sub$habitat==1]), " locations in Habitat A and ", max(sim_sub$count[sim_sub$habitat==2]), " locations in Habitat B"))
   
  # Fit the movement model to the simulated data
  svf <- variogram(sim_sub)
  guess <- ctmm.guess(sim_sub, variogram=svf, interactive=FALSE)
  fit <- ctmm.select(sim_sub, guess, trace=2) #
  summary(fit)
  
  # Calculate the UDs ###
  ud <- akde(sim_sub, fit, weights=TRUE)
  summary(ud)
  
  # Fit the RSFs ###
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
  print(pd)
  print(eTime)
  print(eTime - sTime)
  
}
