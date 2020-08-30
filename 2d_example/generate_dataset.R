set.seed(123)
getwd()
if(!require(COMIX)){
    devtools::install_github("MaStatLab/COMIX")
    library(COMIX)
}

pch = "."

# Draw data from generative model (unperturbed):
# ----------------------------------------------
# Three clusters:
K = 3
# Three samples:
J = 3
# Dimension:
p = 2
# Number of observations per sample and cluster:
njk = matrix(0, J, K)
njk[,1] = c(1420, 1415, 1000)
njk[,2] = c(1165, 1165, 1000)
njk[,3] = c(1670, 1665, 1500)
# Total number of observations:
n = sum(njk)
# Compute weights:
w = njk
for (j in 1:J) {
  w[j,] = w[j,] / sum(w[j,])
}

# Grand locations:
xiGM0 = matrix(c(1, 10, 2,
                 -1, .5, 5),
               nrow = p, ncol = K, byrow=TRUE)

# Covariance of cluster locations around grand locations:
E = array(0, dim = c(p, p, K)) 
E[ , , 1] = matrix(c(0.4, 0.04, 0.04, 0.4), nrow = 2, ncol = 2)
E[ , , 2] = matrix(c(0.3, 0.03, 0.03, 0.6), nrow = 2, ncol = 2)
E[ , , 3] = matrix(c(0.6, 0.05, 0.05, 0.4), nrow = 2, ncol = 2)

# Draw cluster locations (per sample):
xi0 = array(0, dim=c(J, p, K))
for (k in 1:K) {
  xi0[ , , k] = MASS::mvrnorm(J, xiGM0[, k], E[ , , k])
}

# Set cluster scale parameter (same for all samples):
Omega0 = array(0, dim = c(p, p, K))
Omega0[ , , 1] = matrix(c(2, 0.5, 0.5, 1), nrow = 2, ncol = 2)
Omega0[ , , 2] = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
Omega0[ , , 3] = matrix(c(2, 0.5, 0.5, 1), nrow = 2, ncol = 2)

# Set cluster skew parameter (same for all samples):
alpha0 = matrix(c(-6,-7, 1,
                  8, 5, 4),
                nrow = p, ncol = K, byrow = TRUE)

# Compute transformed parameters:
psi0 = matrix(0, nrow = p, ncol = K)
Sigma0 = array(0, dim = c(p, p, K))

# Compute means:
meanvec = array(0, c(J, p, K))
GMvec = matrix(0, p, K)
for (k in 1:K) {
  del.om = COMIX::transform_params(Omega0[ , , k], alpha0[, k])
  psi0[ , k] = del.om$psi
  Sigma0[ , , k] = del.om$Sigma
  for (j in 1:J) {
    meanvec[j, , k] = xi0[j, , k] + del.om$omega %*% del.om$delta * sqrt(2 / pi)
  }
  GMvec[ , k] = xiGM0[ , k] + del.om$omega %*% del.om$delta * sqrt(2 / pi)
}

# Draw observations from multivariate skew-normal latent representation:
Y = NULL
# Sample index:
C = NULL
# Cluster assignment:
Z = NULL
epsilon0 = matrix(0, nrow = p, ncol = n)
cntr = 1
for (k in 1:K) {
  for (j in 1:J) {
    C = c(C, rep(j, njk[j, k]))
    Z = c(Z, rep(k, njk[j, k]))
    z0 = truncnorm::rtruncnorm(n, a = 0, b = Inf, mean = 0, sd = 1)
    if (njk[j, k]>0) {
      for (i in 1:njk[j, k]) {
        et = MASS::mvrnorm(1, mu = c(0, 0), Sigma = Sigma0[ , , k])
        Y =  rbind(Y, xi0[j, , k] + psi0[ , k] * z0[cntr] + et)
        epsilon0[,cntr] = et
        cntr = cntr + 1
      }
    }
  }
}
n = nrow(Y)


mcol = function(col, alpha = 1) {
  return( rgb(red = (col == 1), green = 0.7 * (col == 2), blue = (col == 3), alpha) )
}

# Plot unperturbed data, top row of Figure 3:
pdf("Figure_3.pdf", width = 6, height = 4)
gg = 0.7
par(mfrow = c(2, J))
par(oma = c(0, 0, 1.5, 1.5))
par(mar = c(0, 0, 0, 0))
par(bg = "white")
xlim = range(Y[ , 1])
ylim = range(Y[ , 2])
col = c(2, 3, 4)
for (j in 1:J) {
  plot(
    Y[C == j, 1],
    Y[C == j, 2],
    pch = pch,
    xlab = "",
    ylab = "",
    col = rainbow(K)[Z[C == j]],
    cex = 1.5,
    xlim = xlim,
    ylim = ylim,
    xaxt = "n",
    yaxt = "n",
    asp = 1,
    panel.first = {
      points(0, 0, pch = 16, cex = 1e4, col = rgb(.92, .92, .92));
      grid(col = "white", lty=1)
      },
    )
  mtext(paste0("Sample #", j), side = 3, outer = TRUE, at = grconvertX(0.5, 'npc', 'nic'))
  
  mtext('"Ideal" data', side = 4, outer = FALSE, line = 0.5, cex = 0.8)
  
  abline(v = GMvec[1, ], lty = 2, col = col)
  abline(h = GMvec[2, ], lty = 2, col = col)
  for (k in 1:K) {
    if (njk[j, k] > 0) {
      points(GMvec[1, k], GMvec[2, k], pch = 20, col = col[k], cex = 1.7)
      points(GMvec[1, k], GMvec[2, k], pch = 13, cex = 1.6)
      
      points(meanvec[j, 1, k], meanvec[j, 2, k], pch = 20, col = col[k], cex = 1.5)
      points(meanvec[j, 1, k], meanvec[j, 2, k], cex = 1.4)
    }
  }
}

#--------------------------
# # Generate distorted data: 
Ynoisy = Y
Cnoisy = C
Znoisy = Z

for (j in 1:J) {
  # Distort first cluster in sample j:
  msk = (Z == 1 & C == j)
  temp = Y[msk ,]
  min_x = min(temp[ , 1])
  max_x = max(temp[ , 1])
  min_y = min(temp[ , 2])
  max_y = max(temp[ , 2])
  
  temp[ ,1] = temp[ ,1] - min_x
  temp[ ,2] = temp[ ,2] - min_y
  
  if (j==1) a = 0.02
  if (j==2) a = 0.02
  if (j==3) a = 0.025
  
  temp[ ,2] = temp[ ,2] * ( - a * temp[ , 1] ^ 2 + a * (max_x - min_x + 4) * temp[ , 1])
  temp[ ,1] = temp[ ,1] + min_x
  temp[ ,2] = temp[ ,2] + min_y
  
  Ynoisy[msk, ] = temp
  
  # Distort second cluster in sample j:
  msk = (Z == 2 & C == j)
  temp = Y[msk ,]

  min_x = min(temp[ , 1])
  max_x = max(temp[ , 1])
  min_y = min(temp[ , 2])
  max_y = max(temp[ , 2])
  
  temp[ ,1] = temp[ ,1] - min_x
  temp[ ,2] = temp[ ,2] - min_y
  
  a = 0.06

  temp[ , 2] = temp[ , 2] * ( - a * temp[ , 1] ^ 2 + a * (max_x - min_x + 2) * temp[ , 1])
  temp[ , 1] = temp[ , 1] + min_x
  temp[ , 2] = temp[ , 2] + min_y
  
  Ynoisy[msk, ] = temp
  
  # Distort third cluster in sample j:
  msk = (Z==3 & C==j)
  temp = Y[msk ,]

  min_x = min(temp[ , 1])
  max_x = max(temp[ , 1])
  min_y = min(temp[ , 2])
  max_y = max(temp[ , 2])
  
  temp[ ,1] = temp[ ,1] - min_x
  temp[ ,2] = temp[ ,2] - min_y
  
  if (j==1) a = 0.025
  if (j==2 || j==3) a = 0.04
  
  temp[ ,2] = temp[ ,2] * ( - a * temp[ , 1] ^ 2 + a * (max_x - min_x + 0.8) * temp[ , 1])
  temp[ ,1] = temp[ ,1] + min_x 
  temp[ ,2] = temp[ ,2] + min_y
  Ynoisy[msk, ] = temp
}

# Plot distorted data, bottom row of Figure 3:
par(mar = c(0, 0, 0, 0))
for (j in 1:J) {
  plot(
    Ynoisy[Cnoisy == j, 1],
    Ynoisy[Cnoisy == j, 2],
    pch = pch,
    xlab = "",
    ylab = "",
    col = rainbow(K)[Znoisy[Cnoisy == j]],
    cex = 1.5,
    xlim = xlim,
    ylim = ylim,
    xaxt = "n",
    yaxt = "n",
    asp = 1,
    panel.first = {
      points(0, 0, pch = 16, cex = 1e4, col = rgb(.92, .92, .92));
      grid(col = "white", lty=1)
      },
  )
  mtext('Distorted data', side = 4, outer = FALSE, line = 0.5, cex = 0.8)
  
  abline(v = GMvec[1, ], lty = 2, col = col)
  abline(h = GMvec[2, ], lty = 2, col = col)
}
dev.off()

saveRDS(
   list(
     Y = Y,
     C = C,
     Ynoisy = Ynoisy,
     Cnoisy = Cnoisy
     ),
   "./raw_data.rds"
  )
