# PACKAGES ####
library(raster)
library(ctmm)

# ANALYSIS ####
# Read job number from command line
# args = commandArgs(trailingOnly=TRUE)
sim_no <- 1 # args[1]
# print(sim_no)

# Define the OUF model parameters ###
# Length of a day in seconds
ds <- 86400

# Spatial variance in m^2
sig <- 200000

# Specify an OUF model for simulation
mod <- ctmm(tau=c(ds), isotropic=TRUE, sigma=sig, mu=c(0,0))

# Simulation with varying sampling interval ####

# Sampling frequencies to quantify
samp <- c(25, 50, 100, 200, 400, 800, 1600)

# Create an empty data.frame for saving results
name_df <- c("sim_no","samp_freq", "iid_coef", "wrsf_coef", "iid_lcl", "iid_ucl", "wrsf_lcl", "wrsf_ucl", "runtime")
df_sims <- array(rep(NaN), dim = c(0, length(name_df)))
colnames(df_sims) <- name_df

# Create raster
r1 <- raster(nrows = 1000, ncols = 1000, xmn = -10000, xmx = 10000, ymn = -10000, ymx = 10000, 
             vals = as.factor(c(rep(1,500000),rep(2,500000))))
projection(r1) <- "+proj=aeqd +lon_0=0 +lat_0=0 +datum=WGS84 +units=m"

# Record start time to monitor how long replicates take to compute
sTime <- Sys.time()
print(sTime)

# Loop over sampling frequencies (samp)
for(i in 1:length(samp)){
  
  # Specify variables to manipulate sampling frequency while holding duration constant
  nd <- samp[i] # number of days
  pd <- 24 # number of sampled points per day
  
  # Sampling schedule
  st <- 1:(nd*pd)*(ds/pd) 
  
  # Simulate from the movement model ###
  sim <- simulate(mod, t=st, complete = TRUE)
  df <- data.frame(sim)
  df$latitude <- ifelse(df$latitude < 0 & df$latitude > -0.001, 
                       0 + (0 - df$latitude), df$latitude)
  pts <- df[,4:5]
  sp <- SpatialPoints(pts, proj4string = CRS("+proj=aeqd +lon_0=0 +lat_0=0 +datum=WGS84 +units=m"))
  df$habitat <- raster::extract(r1, spTransform(sp, CRS(projection(r1))))
  df$count <- ave(df$habitat==df$habitat, df$habitat, FUN=cumsum)
  df <- df[df$habitat==2 | df$count %% 2,]
  sim_sub <- subset(sim, row.names(sim) %in% row.names(df))
  sim_sub$y <- ifelse(sim_sub$y < 0 & sim_sub$y > -111, 
                      0 + (0 - sim_sub$y), sim_sub$y)
  plot(sim_sub)
  xlim <- c(0,2 %#% "day")
  plot(variogram(sim_sub))
  
  # Fit the movement model to the simulated data
  fit <- ctmm.fit(sim_sub, CTMM=mod, control=list(method="pNewton")) #
  print("Fitted movement model")
  
  # Calculate the UDs ###
  ud <- akde(sim_sub, fit, weights=TRUE)
  print("UD created")
  
  sTime <- Sys.time() 

  # Fit the RSFs ###
  rsf_sw <- ctmm:::rsf.fit(sim_sub, UD=ud, R=list(test=r1), debias=TRUE, error=0.05)
  summary(rsf_sw)
  print("Fitted RSF")
  
  eTime <- Sys.time()
  
  # Extract variables of interest ###
  sim_no <- sim_no
  samp_dur<- nd
  wrsf_coef <- summary(rsf_sw)$CI[1,2]
  wrsf_lcl <- summary(rsf_sw)$CI[1,1]
  wrsf_ucl <- summary(rsf_sw)$CI[1,3]
  runtime <- difftime(eTime, sTime, units="mins")
  count1 <- sum(df$habitat[df$habitat==1])
  count2 <- sum(df$habitat[df$habitat==2])/2
 
  #################################
  # Vector of results to return
  x <- data.frame(sim_no, samp_dur, wrsf_coef, wrsf_lcl, wrsf_ucl, runtime)
  
  # Store results in data.frame
  write.table(x, 'results/wrsf_sim_results_hi_wrsf_selection.csv', append=TRUE, row.names=FALSE, col.names=FALSE, sep=',') 
  
  # Print indicators of progress
  print(nd)
  print(eTime)
  print(eTime - sTime)
  
}
