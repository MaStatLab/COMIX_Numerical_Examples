# Set to TRUE to perform model fit:
do_estimate = TRUE
# Set to TRUE to plot calibrated data:
plot_compare = TRUE

getwd()

if(!require(COMIX)){
    devtools::install_github("MaStatLab/COMIX")
    library(COMIX)
}

raw_data = readRDS(file="./raw_data.rds")
Y = raw_data$Y
C = raw_data$C
J = length(unique(C))

# Range of tuning parameter zeta:
zeta_seq = c(0.2, 1)
# Fit model for each zeta:
if (do_estimate) {
  for (zeta in zeta_seq) {
    print(paste0("zeta = ", zeta))
    #ns = 500
    ns = 50
    #pmc = list(npart = 5, nburn = 3500, nsave = ns, nskip = 1, ndisplay = 100)
    pmc = list(npart = 5, nburn = 0, nsave = ns, nskip = 1, ndisplay = 100)
    # Model parameters:
    prior = list(K = 150, merge_par = 0.1, zeta = zeta, merge_step = TRUE)
    seed = as.numeric(zeta * 100)
    set.seed(seed)
    try({
      # Fit model:
      res = COMIX::comix(Y, C, prior = prior, pmc = pmc)
      # Handle label switching:
      resRelab = COMIX::relabelChain(res)
      # Calibrate:
      cal = COMIX::calibrate(resRelab)
      K_num = numeric(J)
      cntr = 1
      p = NCOL(res$data$Y)
      C = res$data$C
      sum_sv = COMIX::summarizeChain(resRelab)
      for (j in 1:J) {
        K_num[j] = length(unique(sum_sv$t[C == j]))
      }
            
      rr = 
        list(
          J = J,
          zeta = zeta,
          K_num = K_num,
          sum_chain = sum_sv,
          Y_cal = cal$Y_cal,
          seed = seed,
          prior = prior,
          pmc = pmc
          )
      
      saveRDS(rr, paste0("./fitted_base/cal_data_zeta_", zeta, ".rds"))
    })
  }
} else {
  lz = length(zeta_seq)
  K_num = matrix(NA_real_, nrow = lz, ncol = J)
  idx = 1
  for (zeta in zeta_seq) {
    try({
      rr = readRDS(paste0("./fitted_base/cal_data_zeta_", zeta, ".rds"))
      K_num[idx,] = rr$K_num
    })
    idx = idx + 1
  }
}

#================

plot_all <- function(rr, Y, C, plot_means) {
  gg = 0.4
  ust = unique(rr$sum_chain$t)
  lust = length(ust)
  par(bg = "white")
  xlim = range(Y[ , 1])
  ylim = range(Y[ , 2])
  for (j in 1:J) {
    plot(
      Y[C == j, 1],
      Y[C == j, 2],
      pch = ".",
      xlab = "",
      ylab = "",
      col = rgb(gg, gg, gg),
      cex = 1.5,
      xlim = xlim,
      ylim = ylim,
      xaxt = "n",
      yaxt = "n",
      asp = 1,
      panel.first = {
        points(0, 0, pch = 16, cex = 1e4, col = rgb(.92,.92,.92));
        grid(col = "white", lty = 1)
      }
    )
    mtext(paste0("Sample #", j), side = 3, outer = TRUE, at = grconvertX(0.5, 'npc', 'nic'))
    points(
      rr$Y_cal[C == j, 1],
      rr$Y_cal[C == j, 2],
      pch = ".",
      col = rainbow(length(unique(ord_lbl)))[ord_lbl[C == j]]
    )
    abline(
      v = rr$sum_chain$meanvec0[1, ust],
      lty = 2,
      col = rainbow(length(unique(ord_lbl)))[unique(ord_lbl)]
    )
    abline(
      h = rr$sum_chain$meanvec0[2, ust],
      lty = 2,
      col = rainbow(length(unique(ord_lbl)))[unique(ord_lbl)]
    )
    
    if (plot_means) {
      points(
        rr$sum_chain$meanvec0[1, ust],
        rr$sum_chain$meanvec0[2, ust],
        pch = 20,
        col = rainbow(length(unique(ord_lbl)))[unique(ord_lbl)],
        cex = 1.7
      )
      points(
        rr$sum_chain$meanvec0[1, ust],
        rr$sum_chain$meanvec0[2, ust],
        pch = 13,
        col = rgb(0, 0, 0),
        cex = 1.5
      )
    }
  }
  mtext(
    substitute(paste(zeta, " = ", zt), list(zt = zeta_seq[idx])),
    side = 4,
    outer = FALSE,
    line=0.5,
    cex = 0.8
  )
  
  return(invisible(NULL))
}

#================

# Figure 4:
pdf("Figure_4.pdf", width = 6, height = 4)
winner_idx = c(1,2)
if (plot_compare) {
  par(mfrow = c(length(winner_idx), J))
  par(oma=c(0, 0, 1.5, 1.5))
  par(mar = c(0, 0, 0, 0))
  
  for (idx in winner_idx) {
    rr = readRDS(paste0("./fitted_base/cal_data_zeta_",
                        zeta_seq[idx],
                        ".rds"))
    Y_cal = rr$Ycal
    ust = unique(rr$sum_chain$t)
    ord_lbl = rr$sum_chain$t
    
    for (lbl in ust) {
      ord_lbl[rr$sum_chain$t == lbl] = which(lbl == ust)
    }
    
    # Relabel to get color coding consistency:
    if (rr$zeta == 0.2) {
      temp = ord_lbl
      temp[ord_lbl == 2] = 3
      temp[ord_lbl == 3] = 2
      ord_lbl = temp
    }

    plot_all(rr, Y, C, T)
  }
}
dev.off()
