# Metadata for Alston, Fleming, et al. 2023 Methods in Ecology and Evolution
Repository for code underlying the manuscript [Preprint](https://doi.org/10.1101/2022.04.21.489059)

Files ending in ".R" are R scripts used in simulations and empirical examples

Files ending in ".sh" are bash scripts used for running R scripts in parallel on an HPC

File names denote their contents:

## Data

1. caracal_example.csv: movement data for a caracal
2. KZN_2017_Landuse_latlong.tif: remotely sensed land cover data
3. mongoose_example.csv: movement data for a mongoose
4. serval_example.csv: movement data for a serval

The movement data (1, 3, and 4 above) can only be used to replicate the results presented in Alston, Fleming, et al. Methods in Ecology and Evolution. Any other use of this data requires the written permission of Colleen Downs and Tharmalingam Ramesh (the caracal and serval data) or Colleen Downs and Jarryd Streicher (the mongoose data).

## Scripts
1. caracal_example.R: R script for caracal empirical example (Fig. 4B).
2. caracal_subsample_example.R: R script for caracal empirical example of the effect of subsampling data (Fig. 5).
3. mongoose_example.R: R script for mongoose empirical example (Fig. 4A).
4. run_caracal_example.sh: bash script for running the caracal empirical example on an HPC (Script #1).
5. run_caracal_subsample_example.sh: bash script for running the caracal empirical example of the effect of subsampling data (Script #2).
6. run_mongoose_example.sh: bash script for running the mongoose empirical example on an HPC (Script #3).
7. run_serval_example.sh: bash script for running the serval empirical example on an HPC (Script #16).
8. run_wrsf_sims_hi_iid.sh: bash script for running the IID, clustered-habitat simulation example on an HPC (Script #17).
9. run_wrsf_sims_hi_wrsf.sh: bash script for running the weighted, clustered-habitat simulation example on an HPC (Script #18).
10. run_wrsf_sims_lo_iid.sh: bash script for running the IID, unclustered-habitat simulation example on an HPC (Script #19).
11. run_wrsf_sims_lo_wrsf.sh: bash script for running the weighted, unclustered-habitat simulation example on an HPC (Script #20).
12. run_wrsf_sims_lo_wrsf_selection_duration.sh: bash script for running the simulation example on the effect of altering sampling duration on an HPC (Script #21).
13. run_wrsf_sims_lo_wrsf_selection_frequency.sh: bash script for running the simulation example on the effect of altering sampling frequency on an HPC (Script #22).
14. run_wrsf_sims_lo_wrsf_selection_subsample.sh: bash script for running the simulation example of the effect of subsampling data on an HPC (Script #23).
15. run_wrsf_sims_lo_wrsf_selection_subsample_w.sh: bash script for running the weighted RSF to compare with the example of the effect of subsampling data on an HPC (Script #24).
16. serval_example.R: R script for serval empirical example (Fig. 4C).
17. wrsf_sims_hi_iid.R: R script for running the IID, clustered-habitat simulation example on an HPC (Fig. 3B).
18. wrsf_sims_hi_wrsf.R: R script for running the weighted, clustered-habitat simulation example on an HPC (Fig. 3A).
19. wrsf_sims_lo_iid.R: R script for running the IID, unclustered-habitat simulation example (Fig. 3D).
20. wrsf_sims_lo_wrsf.R: R script for running the weighted, unclustered-habitat simulation example (Fig. 3C).
21. wrsf_sims_lo_wrsf_selection_duration.R: R script for running the simulation example on the effect of altering sampling duration in a weighted resource selection function (Fig. 1B).
22. wrsf_sims_lo_wrsf_selection_frequency.R: R script for running the simulation example on the effect of altering sampling frequency in a weighted resource selection function (Fig. 1A).
23. wrsf_sims_lo_wrsf_selection_subsample.R: R script for running the simulation example of the effect of subsampling data on resource selection functions (Fig. 2).
24. wrsf_sims_lo_wrsf_selection_subsample_w.R: R script for running the weighted RSF to compare with the example of the effect of subsampling data on resource selection functions (Fig. 2).
