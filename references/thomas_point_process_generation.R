library(spatstat)

# 1️⃣ Define window
wx <- 250
wy <- 125
W <- owin(c(0, wx), c(0, wy))

# 2️⃣ Set parameters for stronger clustering
kappa <- 0.002   # fewer cluster centers
scale <- 5       # very small cluster radius
mu <- 50         # many offspring per parent (dense clusters)

# 3️⃣ Simulate Thomas process
set.seed(123)
pp_clustered <- rThomas(kappa = kappa, scale = scale, mu = mu, win = W)

# 4️⃣ Plot it
plot(pp_clustered, main = "Highly Clustered Thomas Process", cols = "darkblue", pch = 16)

# 5️⃣ Optionally, overlay cluster density
dens <- density(pp_clustered, sigma = 10)
plot(dens, main = "Estimated Intensity (Density)")
plot(pp_clustered, add = TRUE, pch = 16, cex = 0.6)
