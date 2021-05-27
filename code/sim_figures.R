require(reshape2)
require(ggplot2)

lb <- read.csv("~/NPI/code_output/tables/ttest.csv")

# Sample size vs. overdispersion figure ----------------------------------------
param_sets <- list(c(1.5, 1000, 0.4, 100),
                   c(1.5, 100, 0.4, 10),
                   c(2, 1000, 0.4, 100),
                   c(2, 10000, 0.4, 100),
                   c(2, 100, 0.4, 10),
                   c(1.5, 1000, 0.2, 100),
                   c(1.5, 100, 0.2, 10),
                   c(2, 1000, 0.2, 100),
                   c(2, 10000, 0.2, 100),
                   c(2, 100, 0.2, 10),
                   
                   c(1.5, 1000, 0.4, 1000),
                   c(1.5, 1000, 0.2, 1000),
                   c(2, 1000, 0.4, 1000),
                   c(2, 1000, 0.2, 1000),
                   c(2, 10000, 0.4, 1000),
                   c(2, 10000, 0.2, 1000),
                   
                   c(1.5, 100, 0.4, 50),
                   c(1.5, 100, 0.2, 50),
                   c(1.5, 100, 0.4, 100),
                   c(1.5, 100, 0.2, 100),
                   c(2, 100, 0.4, 50),
                   c(2, 100, 0.2, 50),
                   c(2, 100, 0.4, 100),
                   c(2, 100, 0.2, 100)
                   )

toplot_ls <- list()
toplot_ls_dex <- 1
for (param_set in param_sets) {
  new_row <- data.frame(matrix(ncol=4, nrow=1))
  colnames(new_row) <- c("one", "four", "seven", "Effect")
  one <- lb[which(lb$R0 == param_set[1] & lb$N == param_set[2] & lb$effect == param_set[3] & lb$nsample == param_set[4] & lb$k == 0.1),]$sufficient_ncluster
  new_row$one <- as.numeric(ifelse(one == ">=1000", 1000, as.numeric(one)))
  four <- lb[which(lb$R0 == param_set[1] & lb$N == param_set[2] & lb$effect == param_set[3] & lb$nsample == param_set[4] & lb$k == 0.4),]$sufficient_ncluster
  new_row$four <- as.numeric(ifelse(four == ">=1000", 1000, as.numeric(four)))
  seven <- lb[which(lb$R0 == param_set[1] & lb$N == param_set[2] & lb$effect == param_set[3] & lb$nsample == param_set[4] & lb$k == 0.7),]$sufficient_ncluster
  new_row$seven <- as.numeric(ifelse(seven == ">=1000", 1000, as.numeric(seven)))
  new_row$Effect <- param_set[3]
  toplot_ls[[toplot_ls_dex]] <- new_row
  toplot_ls_dex <- toplot_ls_dex + 1
}
toplot <- do.call(rbind, toplot_ls)
toplot <- toplot[-which(toplot$one == 1000 & toplot$four == 1000 & toplot$seven == 1000),]
#toplot <- toplot[-which(toplot$one == 1000 | toplot$four == 1000 | toplot$seven == 1000),]
toplot <- data.frame(t(toplot))
toplot$x <- c("0.1", "0.4", "0.7", "Effect")
toplot_m <- melt(toplot, id.vars="x")
toplot_m$Effect <- NA
next_eff <- toplot_m$value[4]
for (i in 1:nrow(toplot_m)) {
  if (i %% 4 == 0 & i != nrow(toplot_m)) {
    next_eff <- toplot_m$value[4 + i]
  } else {
    toplot_m$Effect[i] <- next_eff
  }
}
toplot_m <- toplot_m[which(toplot_m$x != 'Effect'),]
toplot_m$Effect <- as.factor(toplot_m$Effect)
#toplot_m$value <- ifelse(toplot_m$value == 1000, NA, toplot_m$value)
ggplot() +
  geom_line(data = toplot_m, aes(x = x, y = value, group = variable, colour=Effect), linetype='dashed')+theme_bw()+
  geom_point(data = toplot_m, aes(x = x, y = value, group = variable, colour=Effect), size=2)+
  labs(title = "Required Sample Size vs. Overdispersion") + labs(y = "Number of Clusters Per Arm") + 
  labs(x = "Overdispersion (k)")+
  scale_x_discrete(expand = c(0.05, 0.05))


# Community size ---------------------------------------------------------------
param_sets <- list(c(1.5, 0.4, 0.2),
                   c(1.5, 0.7, 0.2),
                   c(2, 0.1, 0.2),
                   c(2, 0.4, 0.2),
                   c(2, 0.7, 0.2),
                   c(1.5, 0.4, 0.4),
                   c(1.5, 0.7, 0.4),
                   c(2, 0.1, 0.4),
                   c(2, 0.4, 0.4),
                   c(2, 0.7, 0.4))
toplot_ls <- list()
toplot_ls_dex <- 1
for (param_set in param_sets) {
  new_row <- data.frame(matrix(ncol=3, nrow=1))
  colnames(new_row) <- c("one", "two", "three")
  one <- lb[which(lb$R0 == param_set[1] & lb$N == 100 & lb$k == param_set[2] & lb$effect == param_set[3] & lb$nsample == 100),]$sufficient_ncluster
  new_row$one <- as.numeric(ifelse(one == ">=1000", 1000, as.numeric(one)))
  two <- lb[which(lb$R0 == param_set[1] & lb$N == 1000 & lb$k == param_set[2] & lb$effect == param_set[3] & lb$nsample == 100),]$sufficient_ncluster
  new_row$two <- as.numeric(ifelse(two == ">=1000", 1000, as.numeric(two)))
  three <- lb[which(lb$R0 == param_set[1] & lb$N == 10000 & lb$k == param_set[2] & lb$effect == param_set[3] & lb$nsample == 100),]$sufficient_ncluster
  new_row$three <- as.numeric(ifelse(three == ">=1000", 1000, as.numeric(three)))
  toplot_ls[[toplot_ls_dex]] <- new_row
  toplot_ls_dex <- toplot_ls_dex + 1
}
toplot <- do.call(rbind, toplot_ls)
toplot <- toplot[-which(toplot$one == 1000 | toplot$two == 1000 | toplot$three == 1000),]
toplot <- data.frame(t(toplot))
toplot$x <- c("100", "1000", "10000")
toplot_m <- melt(toplot, id.vars="x")
ggplot() +
  geom_line(data = toplot_m, aes(x = x, y = value, group = variable), col='dodgerblue2')+theme_bw()+
  labs(title = "Required Sample Size vs. Community Size") + labs(y = "Number of Clusters Per Arm") + 
  labs(x = "Community Size (N)")+
  scale_x_discrete(expand = c(0.05, 0.05))


# Number of clusters needed if we wait three generations -----------------------
lb2 <- read.csv("~/NPI/code_output/tables/lb_gen2.csv")
lb3 <- read.csv("~/NPI/code_output/tables/lb_gen3.csv")
param_sets <- list()
param_sets_dex <- 1
R0s <- c(1.5, 2)
Ns <- c(1000)
ks <- c(0.1, 0.4, 0.7)
effects <- c(0.2, 0.4)
nsamples <- c(100, 1000)
for (R0 in R0s) {
  for (N in Ns) {
    for (k in ks) {
      for (effect in effects) {
        for (nsample in nsamples) {
          if (!(k == 0.1 & N == 10000 & R0 == 1.5)) {
            param_sets[[param_sets_dex]] <- c(R0, N, k, effect, nsample)
            param_sets_dex <- param_sets_dex + 1
          }
        }
      }
    }
  }
}
Ns <- c(100)
nsamples <- c(10, 50, 100)
for (R0 in R0s) {
  for (N in Ns) {
    for (k in ks) {
      for (effect in effects) {
        for (nsample in nsamples) {
          if (!(k == 0.1 & N == 10000 & R0 == 1.5)) {
            param_sets[[param_sets_dex]] <- c(R0, N, k, effect, nsample)
            param_sets_dex <- param_sets_dex + 1
          }
        }
      }
    }
  }
}
toplot_ls <- list()
toplot_ls_dex <- 1
for (param_set in param_sets) {
  new_row <- data.frame(matrix(ncol=3, nrow=1))
  colnames(new_row) <- c("one", "two", "three")
  one <- lb[which(lb$R0 == param_set[1] & lb$N == param_set[2] & lb$k == param_set[3] & lb$effect == param_set[4] & lb$nsample == param_set[5]),]$sufficient_ncluster
  new_row$one <- as.numeric(ifelse(one == ">=1000", 1000, as.numeric(one)))
  two <- lb2[which(lb2$R0 == param_set[1] & lb2$N == param_set[2] & lb2$k == param_set[3] & lb2$effect == param_set[4] & lb2$nsample == param_set[5]),]$sufficient_ncluster
  new_row$two <- as.numeric(ifelse(two == ">=1000", 1000, as.numeric(two)))
  three <- lb3[which(lb3$R0 == param_set[1] & lb3$N == param_set[2] & lb3$k == param_set[3] & lb3$effect == param_set[4] & lb3$nsample == param_set[5]),]$sufficient_ncluster
  new_row$three <- as.numeric(ifelse(three == ">=1000", 1000, as.numeric(three)))
  toplot_ls[[toplot_ls_dex]] <- new_row
  toplot_ls_dex <- toplot_ls_dex + 1
}
toplot <- do.call(rbind, toplot_ls)
toplot <- toplot[-which(toplot$one == 1000 | toplot$two == 1000 | toplot$three == 1000),]
toplotsm <- toplot[-which(toplot$one > 225 | toplot$two > 225 | toplot$three > 225),]
toplot <- data.frame(t(toplot))
toplot$x <- factor(c("First Gen.", "Second Gen.", "Third Gen."), levels=c("First Gen.", "Second Gen.", "Third Gen."))
toplot_m <- melt(toplot, id.vars="x")
ggplot() +
  geom_line(data = toplot_m, aes(x = x, y = value, group = variable), col='dodgerblue2')+theme_bw()+
  labs(title = "C) Required Sample Size vs. Sampling Time (n=10000)") + labs(y = "Number of Clusters Per Arm") + 
  labs(x = "Sampling Time After Intervention")+
  scale_x_discrete(expand = c(0.05, 0.05))

toplot_t <- data.frame(t(toplot))
toplot_t <- toplot_t[1:(nrow(toplot_t) - 1),]
as.numeric(toplot_t$two[which(toplot_t$two < toplot_t$three)]) - as.numeric(toplot_t$three[which(toplot_t$two < toplot_t$three)])
which((toplot$two < toplot$three))
param_sets[[29]]
param_sets[[30]]
param_sets[[34]]
param_sets[[45]]
param_sets[[47]]
param_sets[[48]]
param_sets[[50]]
param_sets[[52]]
toplotsm <- data.frame(t(toplotsm))
toplotsm$x <- factor(c("First Gen.", "Second Gen.", "Third Gen."), levels=c("First Gen.", "Second Gen.", "Third Gen."))
toplotsm_m <- melt(toplotsm, id.vars="x")
ggplot() +
  geom_line(data = toplotsm_m, aes(x = x, y = value, group = variable), col='dodgerblue2')+theme_bw()+
  labs(title = "Required Sample Size vs. Community Size") + 
  #labs(y = "Number of Clusters Per Arm") + 
  #labs(x = "Community Size (N)")+
  scale_x_discrete(expand = c(0.05, 0.05)) + ylim(0, 225) +
  theme(title = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())



# Sampling more from less or sampling less from more ---------------------------
param_sets <- list(c(1.5, 0.4, 0.2, 10000),
                   c(1.5, 0.7, 0.2, 10000),
                   c(2, 0.1, 0.2, 10000),
                   c(2, 0.4, 0.2, 10000),
                   c(2, 0.7, 0.2, 10000),
                   c(1.5, 0.4, 0.4, 10000),
                   c(1.5, 0.7, 0.4, 10000),
                   c(2, 0.1, 0.4, 10000),
                   c(2, 0.4, 0.4, 10000),
                   c(2, 0.7, 0.4, 10000),
                   
                   c(1.5, 0.1, 0.2, 1000),
                   c(1.5, 0.4, 0.2, 1000),
                   c(1.5, 0.7, 0.2, 1000),
                   c(2, 0.1, 0.2, 1000),
                   c(2, 0.4, 0.2, 1000),
                   c(2, 0.7, 0.2, 1000),
                   c(1.5, 0.1, 0.4, 1000),
                   c(1.5, 0.4, 0.4, 1000),
                   c(1.5, 0.7, 0.4, 1000),
                   c(2, 0.1, 0.4, 1000),
                   c(2, 0.4, 0.4, 1000),
                   c(2, 0.7, 0.4, 1000))
toplot_ls <- list()
toplot_ls_dex <- 1
for (param_set in param_sets) {
  new_row <- data.frame(matrix(ncol=3, nrow=1))
  colnames(new_row) <- c("one", "two", "Effect")
  one <- lb[which(lb$R0 == param_set[1] & lb$N == param_set[4] & lb$k == param_set[2] & lb$effect == param_set[3] & lb$nsample == 100),]$sufficient_ncluster
  new_row$one <- as.numeric(ifelse(one == ">=1000", 1000, as.numeric(one) * 100))
  two <- lb[which(lb$R0 == param_set[1] & lb$N == param_set[4] & lb$k == param_set[2] & lb$effect == param_set[3] & lb$nsample == 1000),]$sufficient_ncluster
  new_row$two <- as.numeric(ifelse(two == ">=1000", 1000, as.numeric(two) * 1000))
  new_row$Effect <- param_set[3]
  toplot_ls[[toplot_ls_dex]] <- new_row
  toplot_ls_dex <- toplot_ls_dex + 1
}
toplot <- do.call(rbind, toplot_ls)
toplot <- toplot[-which(toplot$one == 1000 | toplot$two == 1000),]
toplot <- data.frame(t(toplot))
toplot$x <- factor(c("100", "1000", "Effect"), levels=c("100", "1000", "Effect"))
toplot_m <- melt(toplot, id.vars="x")
toplot_m$Effect <- NA
next_eff <- toplot_m$value[3]
for (i in 1:nrow(toplot_m)) {
  if (i %% 3 == 0 & i != nrow(toplot_m)) {
    next_eff <- toplot_m$value[3 + i]
  } else {
    toplot_m$Effect[i] <- next_eff
  }
}
toplot_m <- toplot_m[which(toplot_m$x != 'Effect'),]
toplot_m$Effect <- as.factor(toplot_m$Effect)
#toplot_m$value <- ifelse(toplot_m$value == 1000, NA, toplot_m$value)
ggplot() +
  geom_line(data = toplot_m, aes(x = x, y = value, group = variable, colour=Effect), linetype='dashed')+theme_bw()+
  geom_point(data = toplot_m, aes(x = x, y = value, group = variable, colour=Effect), size=2)+
  labs(title = "Total Required Sample Size vs. Percent Sampled") + labs(y = "Total Required Sample Size") + 
  labs(x = "Number Sampled From Each Cluster")+
  scale_x_discrete(expand = c(0.05, 0.05))

# Showing effect of R0 ---------------------------------------------------------
param_sets <- list(c(100, 0.4, 0.2, 10000),
                   c(1000, 0.4, 0.2, 10000),
                   c(100, 0.7, 0.2, 10000),
                   c(1000, 0.7, 0.2, 10000),
                   c(100, 0.4, 0.4, 10000),
                   c(1000, 0.4, 0.4, 10000),
                   c(100, 0.7, 0.4, 10000),
                   c(1000, 0.7, 0.4, 10000),
                   
                   c(100, 0.1, 0.2, 1000),
                   c(1000, 0.1, 0.2, 1000),
                   c(100, 0.4, 0.2, 1000),
                   c(1000, 0.4, 0.2, 1000),
                   c(100, 0.7, 0.2, 1000),
                   c(1000, 0.7, 0.2, 1000),
                   c(100, 0.1, 0.4, 1000),
                   c(1000, 0.1, 0.4, 1000),
                   c(100, 0.4, 0.4, 1000),
                   c(1000, 0.4, 0.4, 1000),
                   c(100, 0.7, 0.4, 1000),
                   c(1000, 0.7, 0.4, 1000)
                   )
toplot_ls <- list()
toplot_ls_dex <- 1
for (param_set in param_sets) {
  new_row <- data.frame(matrix(ncol=3, nrow=1))
  colnames(new_row) <- c("one", "two", "Effect")
  one <- lb[which(lb$R0 == "1.5" & lb$N == param_set[4] & lb$k == param_set[2] & lb$effect == param_set[3] & lb$nsample == param_set[1]),]$sufficient_ncluster
  new_row$one <- as.numeric(ifelse(one == ">=1000", 1000, as.numeric(one)))
  two <- lb[which(lb$R0 == "2" & lb$N == param_set[4] & lb$k == param_set[2] & lb$effect == param_set[3] & lb$nsample == param_set[1]),]$sufficient_ncluster
  new_row$two <- as.numeric(ifelse(two == ">=1000", 1000, as.numeric(two)))
  new_row$Effect <- param_set[3]
  toplot_ls[[toplot_ls_dex]] <- new_row
  toplot_ls_dex <- toplot_ls_dex + 1
}
toplot <- do.call(rbind, toplot_ls)
toplot <- toplot[-which(toplot$one == 1000 | toplot$two == 1000),]
toplot <- data.frame(t(toplot))
toplot$x <- factor(c("1.5", "2.0", "Effect"), levels=c("1.5", "2.0", "Effect"))
toplot_m <- melt(toplot, id.vars="x")
toplot_m$Effect <- NA
next_eff <- toplot_m$value[3]
for (i in 1:nrow(toplot_m)) {
  if (i %% 3 == 0 & i != nrow(toplot_m)) {
    next_eff <- toplot_m$value[3 + i]
  } else {
    toplot_m$Effect[i] <- next_eff
  }
}
toplot_m <- toplot_m[which(toplot_m$x != 'Effect'),]
toplot_m$Effect <- as.factor(toplot_m$Effect)
#toplot_m$value <- ifelse(toplot_m$value == 1000, NA, toplot_m$value)
ggplot() +
  geom_line(data = toplot_m, aes(x = x, y = value, group = variable, colour=Effect), linetype='dashed')+theme_bw()+
  geom_point(data = toplot_m, aes(x = x, y = value, group = variable, colour=Effect), size=2)+
  labs(title = expression(paste("Required Sample Size vs. ", R[0]))) + labs(y = "Number of Clusters Per Arm") + 
  labs(x = expression(R[0]))+
  scale_x_discrete(expand = c(0.05, 0.05))

