make_data = TRUE

lab_seq = c(1, 2, 3)
nj = 10000

set.seed(123)
if (make_data) {
  J = length(lab_seq)
  prefix = "Exp "
  batch_control_name = "_Patient #1_bc.csv"
  data_path = "./csv_xform/"

  raw_labels =
    paste0(
      c("Dy163", "Gd155", "Nd148", "Yb172", "Gd160", "Dy164", "Nd142", "Nd145", "Nd143", "Nd146",
        "Gd156", "Yb171", "Yb174", "Gd158", "Lu175", "Dy161", "Eu153", "Er166", "Er168"),
      "Di"
    )

  labels = c("BTLA", "CD27", "CD28", "CD38", "CD39", "CD45RO", "CD57", "CD95", "CD127", "CD160",
             "CD244", "GrzB", "HLA.DR", "PD.1", "Perforin", "T.bet", "TIGIT", "TIM.3", "Ki.67")
  
  p = 19
  Y = NULL
  C = NULL
  for (lab in lab_seq) {
    temp = read.csv(file.path(paste0(data_path, prefix, lab, batch_control_name)))
    temp = temp[ , raw_labels]
    colnames(temp) = labels
    idx = sample(1:nrow(temp), nj, replace = F)
    Y_add = as.matrix(temp[idx, ])
    C = c(C, rep(which(lab == lab_seq), nj))
    Y = rbind(Y, Y_add)
  }
  data = list(Y = scale(Y), C = C)
  saveRDS(data, "batch_control_data.rds")
} else {
  data = readRDS("batch_control_data.rds") 
  C = data$C
  Y = data$Y
}
