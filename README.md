# Midfrontal-theta-phase-coordinates-cognitive-control

These scripts were used for the analyses described in the paper 'Midfrontal theta phase coordinates behaviorally relevant brain computations during cognitive control'

The midfrontal_theta_phase.m script performs Generalized Eigen Decomposition to increase midfrontal theta SNR and performs theta phase-resolved power-RT correlations.

The SINFIT.R script fits sine waves (using non-linear least square fitting) to the phase-resolved power-RT correlations and extracts the phase of maximum power-RT correlation for each subject and condition.

The cross_frequency_coupling.m script uses that phase to perform GED enhancing the signal around that phase. The resulting spatially filtered data is then inspected for systematic changes in power at higher frequencies by showing theta-phase-frequency power maps.
This operation is carried using a cross-validation approach so that 90% of the signal is used to compute the spatial filter which is then applied to the 10% remaining data.
This operation is carried-out 10 times and results are finally averaged over the 10 iterations.
Permutations are applied to test for significance.
