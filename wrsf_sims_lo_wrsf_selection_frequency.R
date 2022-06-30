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
mod <- ctmm(tau=c(ds,ds*2), isotropic=TRUE, sigma=sig, mu=c(0,0))

# Simulation with varying sampling interval ####

# Sampling frequencies to quantify
samp <- c(4, 8, 16, 32, 64, 128, 256, 512)

# Create an empty data.frame for saving results
name_df <- c("sim_no","samp_freq", "wrsf_coef", "wrsf_lcl", "wrsf_ucl", "runtime")
df_sims <- array(rep(NaN), dim = c(0, length(name_df)))
colnames(df_sims) <- name_df

# Create raster
r1 <- raster(nrows = 1000, ncols = 1000, 
             xmn = -0.05, xmx = 0.05, ymn = -0.05, ymx = 0.05,
             vals = as.factor(rep(1:2, 500000)))
projection(r1) <- "+proj=longlat +datum=WGS84"

# Record start time to monitor how long replicates take to compute
sTime <- Sys.time()
print(sTime)

# Loop over sampling frequencies (samp)
for(i in 1:length(samp)){
  
  # Specify variables to manipulate sampling frequency while holding duration constant
  set.seed(sim_no) # unique seed for each simulation run
  nd <- 100 # number of days
  pd <- samp[i] # number of sampled points per day
  
  # Sampling schedule
  st <- 1:(nd*pd)*(ds/pd) 
    
  # Simulate from the movement model ###
  sim <- simulate(mod, t=st, complete = TRUE)

  # Create a new track with habitat selection (shift half of points in movement track right by 0.0001 degrees)
  sim_sub <- sim
  sim_sub$longitude <- ifelse(sim_sub$longitude > 0 & floor(sim_sub$longitude*10000) %% 2 == 0, 
                            sim_sub$longitude + 0.0001, sim_sub$longitude)

  df2 <- data.frame(sim_sub)
  pts2 <- df2[,6:7]
  sp2 <- SpatialPoints(pts2, proj4string = CRS(projection(r1)))

  # Check number of points in each habitat
  df2$habitat <- raster::extract(r1, sp2)
   
  # Fit the movement model to the simulated data
  fit <- ctmm.fit(sim_sub, CTMM=mod, control=list(method="pNewton")) #
  print("Fitted movement model")
  
  # Calculate the UDs ###
  ud <- akde(sim_sub, fit, weights=TRUE)
  print("UD created")
  
  # Fit the RSFs ###
  rsf <- ctmm:::rsf.fit(sim_sub, UD=ud, R=list(test=r1), debias=TRUE, error=0.01)
  print("Fitted RSF")  

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
  habitat1 <- sum(df2$habitat[df2$habitat==1])
  habitat2 <- sum(df2$habitat[df2$habitat==2])/2
  
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
