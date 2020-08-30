set.seed(123)
make_tsne = TRUE
plot_tsne = TRUE
plot_marginal_densities = TRUE

if(!require(COMIX)){
    devtools::install_github("MaStatLab/COMIX")
    library(COMIX)
}

tps = c("batch_control_fit")

zeta_seq = NULL
K_num = NULL
seed = NULL
dt = readRDS(file = "./batch_control_data.rds")
Y = dt$Y
C = dt$C
for (tp in tps) {
  fn = dir(paste0("./", tp))
  fn = fn[grepl("zeta", fn)]
  for (f in fn) {
    try({
    rr = readRDS(paste0("./", tp,"/", f))
    zeta_seq = c(zeta_seq, rr$zeta)
    K_num = rbind(K_num, rr$K_num)
    seed = c(seed, rr$seed)
    })
  }
}

oz = order(zeta_seq)
zeta_seq = zeta_seq[oz]
K_num = K_num[oz, ]
# Table S1:
capture.output(
  {
    tbl_s1 = cbind(zeta_seq, K_num)
    colnames(tbl_s1) = c("zeta", paste0("#", 1:10))
    knitr::kable(tbl_s1)
  },
  file = "./TableS1.txt"
)

winner_idx = c(2, 10)
# Compute tsne and save:
if (make_tsne) {
  library(Rtsne)
  tt = list()
  for (idx in winner_idx) {
    print(paste0("zeta = ", zeta_seq[idx]))
    rr =
      readRDS(
        paste0("./", tp, "/cal_flow_zeta_", zeta_seq[idx], ".rds")
        )
    tt[[idx]] = list(tsne = Rtsne(rr$Y_cal), zeta = zeta_seq[idx], idx = idx)
  }
  saveRDS(tt, "./tsne_favorite.rds")
  tto = Rtsne(Y)
  saveRDS(tto, "./tsne_raw.rds")
} else {
  tt = readRDS("./tsne_favorite.rds")
  tto = readRDS("./tsne_raw.rds")
}

# Plot some tsne and save:
nj = 10000
if (plot_tsne) {
  n_row = length(winner_idx) + 1
  J_seq = 1:10
  n_col = length(J_seq)
  
  jpeg(paste0("./Figure_6.jpg"), width = 200 * n_col, height = 200 * n_row)
  par(mfrow = c(n_row, n_col))
  par(mar = c(0.2, 0.2, 0.2, 0.2))
  par(oma = c(0, 0, 0, 3))
  
  tt = readRDS("./tsne_favorite.rds")
  for (j in J_seq) {
    plot(tto$Y[((j - 1) * nj + 1):(j * nj), c(2, 1)],
         xlim = rev(range(tto$Y[((j - 1) * nj + 1):(j * nj), 2])),
         xaxt = "n", yaxt = "n",
         asp = 1, pch = ".", xlab = "", ylab = "")
  }
  mtext("Raw Data", side = 4, outer = FALSE, at = grconvertX(0.5,'npc','nic'), line = 1.5, cex = 1.5)
  
  for (idx in winner_idx) {
    rr = readRDS(paste0("./",tp,"/cal_flow_zeta_",
                        zeta_seq[idx],
                        ".rds"))
    ust = unique(rr$sum_chain$t)
    post_K = length(ust)
    col_range = rainbow(post_K)
    
    ord_lbl = numeric(post_K)
    for (z in ust) {
      ord_lbl[rr$sum_chain$t == z] = which(z == ust)
    }
    for (j in J_seq) {
      plot(
        tt[[idx]]$tsne$Y[((j - 1) * nj + 1):(j * nj), ],
        asp = 1,
        pch = ".", xlab = "", ylab = "",
           xaxt = "n", yaxt = "n",
           col = col_range[ord_lbl[((j - 1) * nj + 1):(j * nj)]])
    }
    mtext(substitute(paste(zeta, " = ", zt), list(zt = zeta_seq[idx])),
          side = 4, outer = FALSE, at = grconvertX(0.5, 'npc', 'nic'),
          line = 1.5, cex = 1.5)
  }
  
  dev.off()
}

file_names = c("Figure_7.pdf", "Figure_S1.pdf", "Figure_S2.pdf")
winner_idx = list(c(2, 10), 1:5, 6:10)
if (plot_marginal_densities) {
  for (b in 1:3) {
    GP = list()
    i = 1
    pdf(
      file_names[b],
      width = 2 * 6,
      height = 2 * (length(winner_idx[[b]]) + 1)
      )
    dt = Y
    title = "Marginal density plots of raw data:"
    nrow = nrow(dt) * ncol(dt)
    markers = colnames(dt)
    dat = data.frame(value = c(dt),
                     marker = rep(markers, each = nrow(dt)),
                     sample = rep(paste("Sample", C), ncol(dt)))
    library(ggplot2)
    library(ggridges)
    
    gp =
      ggplot(dat, aes(x = value, y = sample, fill = sample)) +
      geom_density_ridges(scale = 3) + 
      scale_y_discrete(expand=c(0.01, 0)) +
      scale_x_continuous(expand=c(0.01, 0)) +
      facet_grid(. ~ marker, scales = "free_x") +
      theme_bw() +
      theme(legend.position = "none") +
      labs(y = "", x = "", title = title)
    
    options(repr.plot.width = 10, repr.plot.height = 10)
    GP[[i]] = gp
    i = i + 1
    
    for (idx in winner_idx[[b]]) {
      rr = readRDS(paste0("./",tp,"/cal_flow_zeta_",
                          zeta_seq[idx],
                          ".rds"))
      dt = rr$Y_cal
      if (which(idx == winner_idx[[b]]) == 1) {
      title = "Marginal density plots of calibrated data:"
      } else {
        title = ""
      }
      
      nrow = nrow(dt) * ncol(dt)
      markers = colnames(dt)
      dat = data.frame(value = c(dt),
                       marker = rep(markers, each = nrow(dt)),
                       sample = rep(paste("Sample", C), ncol(dt)))
      library(ggplot2)
      library(ggridges)
      
      gp =
        ggplot(dat, aes(x = value, y = sample, fill = sample)) +
        geom_density_ridges(scale = 3) + 
        scale_y_discrete(expand=c(0.01, 0)) +
        scale_x_continuous(expand=c(0.01, 0)) +
        facet_grid(. ~ marker, scales = "free_x") +
        theme_bw() +
        theme(legend.position = "none") +
        labs(y = substitute(paste(zeta, " = ", zt), list(zt = zeta_seq[idx])), x = "", title = title)

      options(repr.plot.width = 10, repr.plot.height = 10)
      GP[[i]] = gp
      i = i + 1
    }
    library(gridExtra)
    do.call("grid.arrange", c(GP, ncol = 1))
    dev.off()
  }
}
