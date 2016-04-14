load("H:/Research/Completed Projects/2014 Bayesian kernel machine regression (Biostatistics)/Code/simulation/run_on_cluster_h1_h2_h3_scenarios/data/xrfdat.RData")

VarRealistic <- round(cov(xrfdat), 2)
dimnames(VarRealistic) <- NULL

dput(VarRealistic)

