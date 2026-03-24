#devtools::install_github("nspope/corMLPE")
library(corMLPE)
library(nlme)
library(geodist)
library(tidyverse)
library(performance)
library(ggplot2)

# Import coordinates and calculate geographic distance among populations
coords <- read.csv("./examples/GPS.csv", header = TRUE, row.names = 1)
coords.m <- as.matrix(coords)
geo.dist <- as.dist(geodist(coords.m, measure = "geodesic")/1000, diag = TRUE) #converting unit into km
geo.dist.m <- as.matrix(geo.dist)
rownames(geo.dist.m) <- rownames(coords.m)
colnames(geo.dist.m) <- rownames(coords.m)

# Import pairwise Fst among populations
Fst <- read.table("./examples/Fst.csv", sep = ",", header = T, row.names = 1)
Fst.m <- as.matrix(Fst)
gen.dist <- as.dist(Fst.m/(1-Fst.m), diag = TRUE)
gen.dist.m <- as.matrix(gen.dist)
rownames(gen.dist.m) <- rownames(Fst.m)
colnames(gen.dist.m) <- rownames(Fst.m)

# Import pairwise resistance distance among populations
resis <- read.table("./examples/resis_surface/Results/alboobliquatus_resis_dist.csv", sep = ",", header = T, row.names = 1)
resis.m <- as.matrix(resis)
resis.dist <- as.dist(resis.m, diag = TRUE)
resis.dist.m <- as.matrix(resis.dist)
rownames(resis.dist.m) <- rownames(resis.m)
colnames(resis.dist.m) <- rownames(resis.m)

# Import all environmental distance matrices
env.dist.files <- list.files(path = "./examples/clim_data/",
                             pattern = "^pairwise_.*_30s_clipped\\.csv$",
                             full.names = TRUE)

# Function to read and name matrix
read_env_dist <- function(file){
  fname <- basename(file)
  varname <- gsub("pairwise_|_30s_clipped\\.csv", "", fname)
  mat <- as.matrix(read.csv(file, row.names = 1, check.names = FALSE))
  list(name = paste0(varname, ".dist"), matrix = mat)
}

env.dist.list <- lapply(env.dist.files, read_env_dist)

# Convert distance matrices to long format automatically
mat_to_long <- function(mat, varname){
  mat[upper.tri(mat, diag = TRUE)] <- NA # keep lower triangle only
  df <- as.data.frame(as.table(mat))
  df <- na.omit(df)
  colnames(df) <- c("Pop1", "Pop2", varname)
  return(df)
}

Fst.long <- mat_to_long(gen.dist.m, "Fst.trans")
geo.dist.long <- mat_to_long(geo.dist.m, "geo.dist")
resis.dist.long <- mat_to_long(resis.dist.m, "resis.dist")
env.dist.long <- lapply(env.dist.list, function(x){
  mat_to_long(x$matrix, x$name)})

# Merge Fst, geo.dist, resis.dist, and env.dist
df <- reduce(c(list(Fst.long, geo.dist.long, resis.dist.long),
               env.dist.long),
             full_join,
             by = c("Pop1", "Pop2"))

# Scale the predictor variables
pred_cols <- setdiff(names(df), c("Pop1","Pop2","Fst.trans"))
df[pred_cols] <- scale(df[pred_cols])

# Fit MLPE models
Null <- lme(Fst.trans ~ 1, random = ~ 1 | Pop1,
            correlation = corMLPE(form = ~ Pop1 + Pop2),
            method = "ML", data = df) # the Null model

IBD <- lme(Fst.trans ~ geo.dist, random = ~ 1 | Pop1,
           correlation = corMLPE(form = ~ Pop1 + Pop2),
           method = "ML", data = df) # the IBD model

IBR <- lme(Fst.trans ~ resis.dist, random = ~ 1 | Pop1,
           correlation = corMLPE(form = ~ Pop1 + Pop2),
           method = "ML", data = df) # the IBR model

head(df, 0)
clim.vars <- names(df)[6:14] # Grab clim dist
fixed_formula <- reformulate(clim.vars, response = "Fst.trans")

IBE.full <- lme(fixed = fixed_formula, random = ~ 1 | Pop1,
                correlation = corMLPE(form = ~ Pop1 + Pop2),
                method = "ML", data = df)

# Check for collinearity for the environmental factors in the IBE
car::vif(IBE.full)

# Drop the variables showing collinearity
head(df, 0)
reduced_formula.1 <- reformulate(clim.vars[-7:-8], response = "Fst.trans")
reduced_formula.1

# Select the best-fitting environmental factors and update the IBE model accordingly
IBE.preci04 <- lme(Fst.trans ~ preci04.dist, random = ~ 1 | Pop1,
                  correlation = corMLPE(form = ~ Pop1 + Pop2),
                  method = "ML", data = df)

IBE.preci05 <- lme(Fst.trans ~ preci05.dist, random = ~ 1 | Pop1,
                   correlation = corMLPE(form = ~ Pop1 + Pop2),
                   method = "ML", data = df)

IBE.preci06 <- lme(Fst.trans ~ preci06.dist, random = ~ 1 | Pop1,
                   correlation = corMLPE(form = ~ Pop1 + Pop2),
                   method = "ML", data = df)

IBE.preci07 <- lme(Fst.trans ~ preci07.dist, random = ~ 1 | Pop1,
                   correlation = corMLPE(form = ~ Pop1 + Pop2),
                   method = "ML", data = df)

IBE.srad06 <- lme(Fst.trans ~ srad06.dist, random = ~ 1 | Pop1,
                   correlation = corMLPE(form = ~ Pop1 + Pop2),
                   method = "ML", data = df)

IBE.srad08 <- lme(Fst.trans ~ srad08.dist, random = ~ 1 | Pop1,
                  correlation = corMLPE(form = ~ Pop1 + Pop2),
                  method = "ML", data = df)

IBE.reduced.1 <- lme(fixed = reduced_formula.1, random = ~ 1 | Pop1,
                   correlation = corMLPE(form = ~ Pop1 + Pop2),
                   method = "ML", data = df)


AIC(IBE.reduced.1, IBE.preci04, IBE.preci05, IBE.preci06,
    IBE.preci07, IBE.srad06, IBE.srad08)

drop1(IBE.reduced.1)

IBE.reduced.2 <- lme(Fst.trans ~ preci06.dist + srad06.dist + srad08.dist + tmin07.dist,
                     random = ~ 1 | Pop1,
                     correlation = corMLPE(form = ~ Pop1 + Pop2),
                     method = "ML", data = df)

IBE.reduced.3 <- lme(Fst.trans ~ preci06.dist + srad06.dist + srad08.dist,
                     random = ~ 1 | Pop1,
                     correlation = corMLPE(form = ~ Pop1 + Pop2),
                     method = "ML", data = df)

AIC(IBE.reduced.2, IBE.reduced.3, IBE.preci06)

car::vif(IBE.reduced.3) # Check collinearity again

IBE.final <- IBE.reduced.3 # The best-fitting IBE model

# Check the normality of residuals
par(mfrow = c(1,2))
hist(residuals(IBD))
qqnorm(residuals(IBD))
qqline(residuals(IBD))

par(mfrow = c(1,2))
hist(residuals(IBR))
qqnorm(residuals(IBR))
qqline(residuals(IBR))

par(mfrow = c(1,2))
hist(residuals(IBE.final))
qqnorm(residuals(IBE.final))
qqline(residuals(IBE.final))

# Summary of the models
summary(IBD)
summary(IBR)
summary(IBE.final)

# Model comparison and selection
AIC(Null, IBD, IBR, IBE.final)
r2(Null)
r2(IBD)
r2(IBR)
r2(IBE.final)

# Plot the selected model
df$IBD_predic <- predict(IBD, level = 0) # level=0 gives marginal predictions
ggplot(df, aes(geo.dist, Fst.trans)) +
  geom_point(size = 2, alpha = 0.3) +
  geom_line(aes(y = IBD_predic), color = "black") +
  theme_minimal(base_size = 14) +
  labs(x = "Geographic Distance (scaled)", y ="Fst/(1-Fst)",
       title = "Isolation by Distance (IBD) Test")
