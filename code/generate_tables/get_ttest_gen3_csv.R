# ------------------------------------------------------------------------------
# @author: Justin Sheen
#
# @description: script to create csv of ttest results (three generations after intervention)
#
# ------------------------------------------------------------------------------
# Get all parameters -----------------------------------------------------------
param_sets <- c()
param_sets_dex <- 1
rts <- c(1.5, 2)
overdispersions <- c(0.7, 0.4, 0.1)
clusters <- c(1000, 10000)
effects <- c(0.2, 0.4)
eits <- c(0.005)
nsamples <- c(100, 1000)
for (i in clusters) {
  for (j in rts) {
    for (k in eits) {
      for (l in overdispersions) {
        for (m in effects) {
          for (n in nsamples) {
            if (!(i == 10000 & l == 0.1 & j == 1.5)) {
              if (i == 10000 & l == 0.1 & j == 2) {
                param_sets[[param_sets_dex]] <- list(i, j, 0.0045, l , m, n)
                param_sets_dex <- param_sets_dex + 1
              } else {
                param_sets[[param_sets_dex]] <- list(i, j, k, l , m, n)
                param_sets_dex <- param_sets_dex + 1
              }
            }
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
  It_It1con_It1trt <- read.csv(paste0("~/NPI/code_output/res/", tgt_R0, "_", N_cluster, "_", k_overdispersion, "_", effect, "_", expected_It_N, ".csv"), header=F)
  adjustment <- unname(unlist(It_It1con_It1trt[nrow(It_It1con_It1trt),][1] / 3000))
  inter_day <- unname(unlist(It_It1con_It1trt[nrow(It_It1con_It1trt),][4]))
  It_It1con_It1trt <- It_It1con_It1trt[1:(nrow(It_It1con_It1trt) - 1),]
  sufficient_ncluster_res <- read.csv(paste0("~/NPI/code_output/res/res_", tgt_R0, "_", N_cluster, "_", k_overdispersion, "_", effect, "_", expected_It_N, "_", nsample, "_ttest_gen3.csv"), header=F)
  sufficient_ncluster <- unname(unlist(sufficient_ncluster_res))
  if(sufficient_ncluster == 1000 | sufficient_ncluster == -1) {
    sufficient_ncluster <- ">=1000"
    adjusted_sufficient_ncluster <- ">=1000"
  } else {
    adjusted_sufficient_ncluster <- sufficient_ncluster * adjustment
  }
  Rt0s <- c()
  Rt1s <- c()
  Its <- c()
  It1s_con <- c()
  It1s_trt <- c()
  for (sim_num in 1:nrow(It_It1con_It1trt)) {
    It <- It_It1con_It1trt[sim_num,3]
    It1_con <- It_It1con_It1trt[sim_num,15]
    It1_trt <- It_It1con_It1trt[sim_num,27]
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
  lower_bound <- SScalc.LB(Rt1=mean_Rt1, Rt0=mean_Rt0, 
                           k1=k_overdispersion, k0=k_overdispersion, # Need to have read in this function from Lee's script
                           m=nsample, n=N_cluster, EIt=(mean_It / N_cluster),
                           pow=.8, alpha=.05, N=NULL)
  new_row <- data.frame(matrix(nrow=1, ncol=20))
  colnames(new_row) <- c("E_It", "N", "R0", "k", "effect", "inter_day", "mean_Rt0", 
                         "var_Rt0", "mean_Rt1", "var_Rt1", 
                         "mean_It", "var_It", "mean_It1_con", "var_It1_con", "mean_It1_trt", "var_It1_trt",
                         "nsample", "lower_bound",
                         "sufficient_ncluster", "adjusted_sufficient_ncluster")
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
  new_row$lower_bound <- lower_bound
  new_row$sufficient_ncluster <- sufficient_ncluster
  new_row$adjusted_sufficient_ncluster <- adjusted_sufficient_ncluster
  new_row$mean_It <- mean_It
  new_row$var_It <- var_It
  new_row$mean_It1_trt <- mean_It1_trt
  new_row$var_It1_trt <- var_It1_trt
  new_row$mean_It1_con <- mean_It1_con
  new_row$var_It1_con <- var_It1_con
  final_df[[final_df_dex]] <- new_row
  final_df_dex <- final_df_dex + 1
}
final_df <- do.call(rbind, final_df)
write.csv(final_df, "~/NPI/code_output/tables/ttest_gen3.csv")
