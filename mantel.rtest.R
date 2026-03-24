library(ade4)
library(geodist)

# =======================================
# Mantel test for isolation-by-distance
# =======================================
coords <- read.table("./examples/GPS.csv", sep = ",", header = T, row.names = 1)
coords.m <- as.matrix(coords)
geo.dist <- as.dist(geodist(coords.m, measure = "geodesic")/1000, diag = TRUE) # Converting unit into km
geo.dist.m <- as.matrix(geo.dist)
rownames(geo.dist.m) <- rownames(coords.m)
colnames(geo.dist.m) <- rownames(coords.m)
geo.dist <- as.dist(geo.dist.m, diag = TRUE)

Fst <- read.table("./examples/Fst.csv", sep = ",", header = T, row.names = 1)
Fst.m <- as.matrix(Fst)
gen.dist <- as.dist(Fst.m/(1-Fst.m), diag = TRUE) # Linearization

IBD.output <- mantel.rtest(geo.dist, gen.dist, nrepet = 9999, alter = c("greater")) # Note that the row and column ordering of geo.dist and gen.dist must be identical
IBD.output

# =========================================
# Mantel test for isolation-by-environment
# =========================================

clim <- read.table("./examples/clim_data/pairwise_preci04_30s_clipped.csv", sep = ",", header = T, row.names = 1)
clim.m <- as.matrix(clim)
clim.dist <- as.dist(clim.m, diag = TRUE)

Fst <- read.table("./examples/Fst.csv", sep = ",", header = T, row.names = 1)
Fst.m <- as.matrix(Fst)
gen.dist <- as.dist(Fst.m/(1-Fst.m), diag = TRUE) # Linearization

IBE.output <- mantel.rtest(clim.dist, gen.dist, nrepet = 9999, alter = c("greater"))
IBE.output

# =========================================
# Mantel test for isolation-by-resistance
# =========================================
resist <- read.table("./examples/resis_surface/Results/alboobliquatus_resis_dist.csv", sep = ",", header = T, row.names = 1)
resis.dist <- as.dist(as.matrix(resist), diag = TRUE)

gen <- read.table("./examples/Fst.csv", sep = ",", header = T, row.names = 1)
Fst <- as.matrix(gen)
gen.dist <- as.dist(Fst/(1-Fst), diag = TRUE) # Linearization

IBR.output <- mantel.rtest(resis.dist, gen.dist, nrepet = 9999, alter = c("greater"))
IBR.output
