argv = commandArgs(TRUE)
sim = as.numeric(argv)
# :::::::::::::::::::::::::::::::::

zeta = sim/100

getwd()

if(!require(COMIX)){
    devtools::install_github("MaStatLab/COMIX")
    library(COMIX)
}


batch_control = readRDS("batch_control_data.rds")
Y = scale(batch_control$Y)
C = batch_control$C
remove(batch_control)
J = length(unique(C))

ns = 1000
pmc = list(npart = 50, nburn = 5500, nsave = ns, nskip = 1, ndisplay = 100)
    
prior = list(K = 30, merge_par = 0.1, zeta = zeta, merge_step = TRUE)

set.seed(sim)
res = COMIX::comix(Y, C, prior = prior, pmc = pmc)
print("did res")
resRelab = relabelChain(res)
print("did reLab")
cal = calibrateNoDist(resRelab)
print("Did cal")

K_num = numeric(J)
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
	  Y_cal = cal$Y_cal
    )

saveRDS(rr, paste0("./batch_control_fit/cal_flow_zeta_",zeta,".rds"))
