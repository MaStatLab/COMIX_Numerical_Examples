#!/usr/bin/Rscript
argv = commandArgs(TRUE)
sim = as.numeric(argv)
# :::::::::::::::::::::::::::::::::

zeta = sim / 100
seed = sim

if(!require(COMIX)){
    devtools::install_github("MaStatLab/COMIX")
    library(COMIX)
}

batch_control = readRDS("batch_control_data.rds")
Y = scale(batch_control$Y)
C = batch_control$C
remove(batch_control)
J = length(unique(C))

ns = 500
pmc = list(npart = 50, nburn = 3500, nsave = ns, nskip = 1, ndisplay = 50)

prior = list(K = 50, merge_par = 0.1, zeta = zeta, merge_step = T)

set.seed(seed)
res = comix(Y, C, prior = prior, pmc = pmc)
print("did res")
resRelab = relabelChain(res)
print("did reLab")
cal = calibrateNoDist(resRelab)
print("Did cal")

K_num = numeric(J)
cntr = 1
p = NCOL(res$data$Y)
C = res$data$C
sum_sv = summarizeChain(resRelab)
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
    prior = prior,
    pmc = pmc,
    seed = seed
    )

saveRDS(rr, paste0("./batch_control_fit/cal_flow_zeta_",zeta,".rds"))
