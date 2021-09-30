# PACKAGES ####
library(ctmm)
library(raster)

# ANALYSIS ####
# Read job number from command line
args = commandArgs(trailingOnly=TRUE)
sim_no <- args[1]

# Define the OUF model parameters ###
# Length of a day in seconds
ds <- 86400

# Spatial variance in m^2
sig <- 100000

# True 95% range area
trueRngArea <- -2*log(0.05)*pi*sig

# Specify an OUF model for simulation
mod <- ctmm(tau=c(7*ds,ds/4), isotropic=TRUE, sigma=sig, mu=c(0,0))

# Simulation with varying sampling interval ####

# Sampling frequencies to quantify
samp <- c(1, 2, 4, 8, 16, 32, 64)

# Specify the desired number of replicates for each sampling interval
nRep <- rep(1, 7)

# Create an empty data.frame for saving results
name_df <- c("sim_no","samp_freq", "iid_coef", "wrsf_coef", "iid_lcl", "iid_ucl", "wrsf_lcl", "wrsf_ucl", "runtime")
df_sims <- array(rep(NaN), dim = c(0, length(name_df)))
colnames(df_sims) <- name_df

# Create raster
r1 <- raster(nrows = 1000, ncols = 1000, xmn = -4000, xmx = 4000, ymn = -4000, ymx = 4000, 
             vals = as.factor(rep(1:2,(ncell(r1)/2))))

# Record start time to monitor how long replicates take to compute
sTime <- Sys.time()
print(sTime)

# Loop over sampling frequencies (samp)
for(i in 1:length(samp)){
  
  # Specify variables to manipulate sampling frequency while holding duration constant
  nd <- 365 # number of days
  pd <- samp[i] # number of sampled points per day
  
  # Sampling schedule
  st <- 1:(nd*pd)*(ds/pd) 
  nReps <- nRep[i]
  
  # Simulate from the movement model ###
  sim <- simulate(mod, t=st)
  df <- data.frame(sim)
  sp <- df
  coordinates(sp) <- c("x","y")
  df$habitat <- as.factor(raster::extract(r1, sp))
  df$count <- ave(df$habitat==df$habitat, df$habitat, FUN=cumsum)
  df <- df[df$habitat==1 | df$count%%2,]
  
  # Fit the movement model to the simulated data
  fit <- ctmm.fit(sim, CTMM=mod, control=list(method="pNewton")) #
  fit_iid <- ctmm.fit(sim, CTMM=ctmm(isotropic=TRUE), control=list(method="pNewton")) #
  
  # Calculate the UDs ###
  ud <- akde(df, fit, weights=TRUE, trace=2)
  ud_iid <- akde(df, fit_iid)
  
  # Fit the RSFs ###
  rsf <- ctmm:::rsf.fit(df, UD=ud, R=list(test=r1), trace=TRUE, debias=TRUE)
  rsf_iid <- ctmm:::rsf.fit(df, UD=ud_iid, R=list(test=r1), trace=TRUE, debias=TRUE, error=0.025)
  
  eTime <- Sys.time()
  
  # Extract variables of interest ###
  samp_freq <- nReps
  iid_coef <- summary(rsf_iid)$CI[1,2]
  wrsf_coef <- summary(rsf)$CI[1,2]
  iid_lcl <- summary(rsf_iid)$CI[1,1]
  iid_ucl <- summary(rsf_iid)$CI[1,3]
  wrsf_lcl <- summary(rsf)$CI[1,1]
  wrsf_ucl <- summary(rsf)$CI[1,3]
  runtime <- sTime - eTime
  
  #################################
  # Vector of results to return
  x <- data.frame(sim_no, samp_freq, iid_coef, wrsf_coef, iid_lcl, iid_ucl, wrsf_lcl, wrsf_ucl, runtime)
  
  # Store results in data.frame
  write.table(x, 'sims_rsf.csv', append=TRUE, row.names=FALSE, col.names=FALSE, sep=',') 
  
  # Print indicators of progress
  print(pd)
  print(eTime)
  print(eTime - sTime)
  
}