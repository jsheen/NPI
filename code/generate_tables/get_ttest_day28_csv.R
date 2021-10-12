# ------------------------------------------------------------------------------
# @author: Justin Sheen
#
# @description: script to create csv of ttest results (one generation after intervention)
#               when intervention is enacted on day 28 of the simulation. 
#
# ------------------------------------------------------------------------------
# Get all parameters -----------------------------------------------------------
param_sets <- c()
param_sets_dex <- 1
rts <- c(1.5, 2)
overdispersions <- c(0.7, 0.4, 0.1)
clusters <- c(1000)
effects <- c(0.2, 0.4)
eits <- c(0.005)
nsamples <- c(100, 1000)
for (i in clusters) {
  for (j in rts) {
    for (k in eits) {
      for (l in overdispersions) {
        for (m in effects) {
          for (n in nsamples) {
            param_sets[[param_sets_dex]] <- list(i, j, k, l , m, n)
            param_sets_dex <- param_sets_dex + 1
          }
        }
      }
    }
  }
}
rts <- c(1.5, 2)
overdispersions <- c(0.7, 0.4, 0.1)
clusters <- c(100)
effects <- c(0.2, 0.4)
eits <- c(0.02)
nsamples <- c(10, 50, 100)
for (i in clusters) {
  for (j in rts) {
    for (k in eits) {
      for (l in overdispersions) {
        for (m in effects) {
          for (n in nsamples) {
            param_sets[[param_sets_dex]] <- list(i, j, k, l , m, n)
            param_sets_dex <- param_sets_dex + 1
          }
        }
      }
    }
  }
}
final_df <- list()
final_df_dex <- 1
# Iterate through all parameters ----------------------------------------------
for (param_set in param_sets) {
  N_cluster <- param_set[[1]][1]
  tgt_R0 <- param_set[[2]][1]
  expected_It_N  <- param_set[[3]][1]
  k_overdispersion <- param_set[[4]][1]
  effect <- param_set[[5]][1]
  nsample <- param_set[[6]][1]
  It_It1con_It1trt <- read.csv(paste0("~/NPI/code_output/res/", tgt_R0, "_", N_cluster, "_", k_overdispersion, "_", effect, "_", expected_It_N, "_day28.csv"), header=F)
  inter_day <- unname(unlist(It_It1con_It1trt[nrow(It_It1con_It1trt),][4]))
  It_It1con_It1trt <- It_It1con_It1trt[1:(nrow(It_It1con_It1trt) - 1),]
  sufficient_ncluster_res <- read.csv(paste0("~/NPI/code_output/res/res_", tgt_R0, "_", N_cluster, "_", k_overdispersion, "_", effect, "_", expected_It_N, "_", nsample, "_ttest_day28.csv"), header=F)
  sufficient_ncluster <- unname(unlist(sufficient_ncluster_res))
  if (sufficient_ncluster == 1000 | sufficient_ncluster == -1) {
    sufficient_ncluster <- ">=1000"
  }
  Rt0s <- c()
  Rt1s <- c()
  Its <- c()
  It1s_con <- c()
  It1s_trt <- c()
  for (sim_num in 1:nrow(It_It1con_It1trt)) {
    It <- It_It1con_It1trt[sim_num,3]
    It1_con <- It_It1con_It1trt[sim_num,7]
    It1_trt <- It_It1con_It1trt[sim_num,19]
    Its <- c(Its, It)
    It1s_con <- c(It1s_con, It1_con)
    It1s_trt <- c(It1s_trt, It1_trt)
    Rt0s <- c(Rt0s, (It1_con / It))
    Rt1s <- c(Rt1s, (It1_trt / It))
  }
  mean_Rt0 <- mean(Rt0s)
  var_Rt0 <- var(Rt0s)
  mean_Rt1 <- mean(Rt1s)
  var_Rt1 <- var(Rt1s)
  mean_It <- mean(Its)
  var_It <- var(Its)
  mean_It1_trt <- mean(It1s_trt)
  var_It1_trt <- var(It1s_trt)
  mean_It1_con <- mean(It1s_con)
  var_It1_con <- var(It1s_con)
  new_row <- data.frame(matrix(nrow=1, ncol=18))
  colnames(new_row) <- c("E_It", "N", "R0", "k", "effect", "inter_day", "mean_Rt0", 
                         "var_Rt0", "mean_Rt1", "var_Rt1", 
                         "mean_It", "var_It", "mean_It1_con", "var_It1_con", "mean_It1_trt", "var_It1_trt",
                         "nsample",
                         "sufficient_ncluster")
  new_row$E_It <- expected_It_N
  new_row$R0 <- tgt_R0
  new_row$effect <- effect
  new_row$k <- k_overdispersion
  new_row$N <- N_cluster
  new_row$nsample <- nsample
  new_row$inter_day <- inter_day
  new_row$mean_Rt0 <- mean_Rt0
  new_row$var_Rt0 <- var_Rt0
  new_row$mean_Rt1 <- mean_Rt1
  new_row$var_Rt1 <- var_Rt1
  new_row$sufficient_ncluster <- sufficient_ncluster
  new_row$mean_It <- mean_It / N_cluster
  new_row$var_It <- var_It / ((N_cluster)^2)
  new_row$mean_It1_trt <- mean_It1_trt / N_cluster
  new_row$var_It1_trt <- var_It1_trt / ((N_cluster)^2)
  new_row$mean_It1_con <- mean_It1_con / N_cluster
  new_row$var_It1_con <- var_It1_con / ((N_cluster)^2)
  final_df[[final_df_dex]] <- new_row
  final_df_dex <- final_df_dex + 1
}
final_df <- do.call(rbind, final_df)
write.csv(final_df, "~/NPI/code_output/tables/ttest_day28.csv")
