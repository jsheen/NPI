# -----------------------------------------------------------------------------
# @author: Justin Sheen
#
# @description: script used to create the figure showing effect of lag on sample size
# -----------------------------------------------------------------------------

require(reshape2)
require(ggplot2)
require(gridExtra)
lb <- read.csv("~/NPI/code_output/tables/ttest_rerun.csv", stringsAsFactors = F)
lb2 <- read.csv("~/NPI/code_output/tables/ttest_gen2.csv", stringsAsFactors = F)
lb3 <- read.csv("~/NPI/code_output/tables/ttest_gen3.csv", stringsAsFactors = F)
param_sets <- list()
param_sets_dex <- 1
R0s <- c(1.5, 2)
Ns <- c(1000, 10000)
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
  new_row <- data.frame(matrix(ncol=7, nrow=1))
  colnames(new_row) <- c("one", "two", "three", "effect", "R0", "N", "k")
  one <- lb[which(lb$R0 == param_set[1] & lb$N == param_set[2] & lb$k == param_set[3] & lb$effect == param_set[4] & lb$nsample == param_set[5]),]
  new_row$one <- as.numeric(ifelse(one$sufficient_ncluster == ">=1000" | one$sufficient_ncluster == "-1", 1000, one$sufficient_ncluster))
  two <- lb2[which(lb2$R0 == param_set[1] & lb2$N == param_set[2] & lb2$k == param_set[3] & lb2$effect == param_set[4] & lb2$nsample == param_set[5]),]
  new_row$two <- as.numeric(ifelse(two$sufficient_ncluster == ">=1000" | two$sufficient_ncluster == "-1", 1000, two$sufficient_ncluster))
  three <- lb3[which(lb3$R0 == param_set[1] & lb3$N == param_set[2] & lb3$k == param_set[3] & lb3$effect == param_set[4] & lb3$nsample == param_set[5]),]
  new_row$three <- as.numeric(ifelse(three$sufficient_ncluster == ">=1000" | three$sufficient_ncluster == "-1", 1000, three$sufficient_ncluster))
  new_row$effect <- param_set[4]
  new_row$R0 <- param_set[1]
  new_row$N <- param_set[2]
  new_row$k <- param_set[3]
  toplot_ls[[toplot_ls_dex]] <- new_row
  toplot_ls_dex <- toplot_ls_dex + 1
}
toplot <- do.call(rbind, toplot_ls)
toplot$one <- ifelse(toplot$one == 1000, NA, toplot$one)
toplot$two <- ifelse(toplot$two == 1000, NA, toplot$two)
toplot$three <- ifelse(toplot$three == 1000, NA, toplot$three)
Ns <- c(100, 1000, 10000)
dummy <- data.frame(matrix(NA, nrow=4, ncol=5))
colnames(dummy) <- c('x', '20%', '40%')
dummy$x <- 1:4
dummy$y1 <- 1:4
dummy$y2 <- 2:5
dummy$y3 <- 3:6
dummy$y4 <- 4:7
dummy_plot <- ggplot() +
  geom_line(dummy, mapping=aes(x=x, y=y1, color='20%, 1.5', linetype='20%, 1.5'))+
  geom_line(dummy, mapping=aes(x=x, y=y3, color='20%, 2.0', linetype='20%, 2.0'))+
  geom_line(dummy, mapping=aes(x=x, y=y2, color='40%, 1.5', linetype='40%, 1.5'))+
  geom_line(dummy, mapping=aes(x=x, y=y4, color='40%, 2.0', linetype='40%, 2.0'))+
  scale_color_manual(name="Effect Size, "~R[0]^0,
                     labels=c("20%, 1.5", "20%, 2.0","40%, 1.5", "40%, 2.0"),
                     values=c('black', 'black', 'dodgerblue', 'dodgerblue'))+
  scale_linetype_manual(name="Effect Size, "~R[0]^0,
                        labels=c("20%, 1.5", "20%, 2.0","40%, 1.5", "40%, 2.0"),
                        values=c('dashed', 'solid', 'dashed', 'solid'))+
  theme(legend.position='bottom',
        text = element_text(size=18))
#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend <- g_legend(dummy_plot)
ps <- list()
ps_dex <- 1
for (N in Ns) {
  print(N)
  helper_function <- function(toplot, k, N) {
    subplot <- toplot[which(toplot$N == N & toplot$k == k),]
    subplot_20_1.5 <- subplot[which(subplot$effect == "0.2" & subplot$R0 == "1.5"),]
    subplot_20_2 <- subplot[which(subplot$effect == "0.2" & subplot$R0 == "2"),]
    subplot_40_1.5 <- subplot[which(subplot$effect == "0.4" & subplot$R0 == "1.5"),]
    subplot_40_2 <- subplot[which(subplot$effect == "0.4" & subplot$R0 == "2"),]
    subplot_20_1.5 <- data.frame(t(subplot_20_1.5))[1:3,]
    subplot_20_1.5$x <- c('First Gen.', 'Second Gen.', 'Third Gen.')
    subplot_20_1.5 <- melt(subplot_20_1.5, 'x')
    subplot_20_2 <- data.frame(t(subplot_20_2))[1:3,]
    subplot_20_2$x <- c('First Gen.', 'Second Gen.', 'Third Gen.')
    subplot_20_2 <- melt(subplot_20_2, 'x')
    subplot_40_1.5 <- data.frame(t(subplot_40_1.5))[1:3,]
    subplot_40_1.5$x <- c('First Gen.', 'Second Gen.', 'Third Gen.')
    subplot_40_1.5 <- melt(subplot_40_1.5, 'x')
    subplot_40_2 <- data.frame(t(subplot_40_2))[1:3,]
    subplot_40_2$x <- c('First Gen.', 'Second Gen.', 'Third Gen.')
    subplot_40_2 <- melt(subplot_40_2, 'x')
    subplot_20_2$title <- paste0('N = ', N, ", k = ", k)
    if (N == 10000 & k == 0.1) {
      subplot_final <- ggplot() + 
        geom_line(subplot_20_2, mapping=aes(x=x,y=value,group=variable,color='20'), linetype='solid')+
        geom_line(subplot_40_2, mapping=aes(x=x,y=value,group=variable,color='40'), linetype='solid')+
        ylab('Number of Clusters Per Arm')+
        xlab('')+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.position='none',
              text = element_text(size=12))+
        facet_grid(. ~ title)+
        scale_color_manual(values=c('black', 'dodgerblue'))
    } else if (N == 10000 & k == 0.4) {
      subplot_final <- ggplot() + 
        geom_line(subplot_20_1.5, mapping=aes(x=x,y=value,group=variable,color='20'), linetype='dashed')+
        geom_line(subplot_20_2, mapping=aes(x=x,y=value,group=variable,color='20'), linetype='solid')+
        geom_line(subplot_40_1.5, mapping=aes(x=x,y=value,group=variable,color='40'), linetype='dashed')+
        geom_line(subplot_40_2, mapping=aes(x=x,y=value,group=variable,color='40'), linetype='solid')+
        ylab('Number of Clusters Per Arm')+
        xlab('Generation of Sampling')+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.position='none',
              text = element_text(size=12))+
        facet_grid(. ~ title)+
        scale_color_manual(values=c('black', 'dodgerblue'))
    } else {
      subplot_final <- ggplot() + 
        geom_line(subplot_20_1.5, mapping=aes(x=x,y=value,group=variable,color='20'), linetype='dashed')+
        geom_line(subplot_20_2, mapping=aes(x=x,y=value,group=variable,color='20'), linetype='solid')+
        geom_line(subplot_40_1.5, mapping=aes(x=x,y=value,group=variable,color='40'), linetype='dashed')+
        geom_line(subplot_40_2, mapping=aes(x=x,y=value,group=variable,color='40'), linetype='solid')+
        ylab('Number of Clusters Per Arm')+
        xlab('')+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.position='none',
              text = element_text(size=12))+
        facet_grid(. ~ title)+
        scale_color_manual(values=c('black', 'dodgerblue'))
    }

    return(subplot_final)
  }
  subplot_1 <- helper_function(toplot, 0.1, N)
  subplot_4 <- helper_function(toplot, 0.4, N)
  subplot_7 <- helper_function(toplot, 0.7, N)
  if (N == 100) {
    p <- grid.arrange(arrangeGrob(subplot_1 + ylab('') + xlab(''),
                                  subplot_4 + ylab('') + xlab(''),
                                  subplot_7 + ylab('') + xlab(''),
                                  nrow=1),
                      nrow=1,heights=c(10))
  } else if (N == 1000) {
    p <- grid.arrange(arrangeGrob(subplot_1 + xlab(''),
                                  subplot_4 + ylab('') + xlab(''),
                                  subplot_7 + ylab('') + xlab(''),
                                  nrow=1),
                      nrow=1,heights=c(10))
  } else {
    p <- grid.arrange(arrangeGrob(subplot_1 + ylab(''),
                                  subplot_4 + ylab(''),
                                  subplot_7 + ylab(''),
                                  nrow=1),
                      nrow=1, heights=c(10))
  }
  ps[[ps_dex]] <- p
  ps_dex <- ps_dex + 1
}
ggsave(filename="~/NPI/code_output/figs/Fig5A.png", plot=ps[[1]],
       width=12, height=3, units="in", dpi=600)
ggsave(filename="~/NPI/code_output/figs/Fig5B.png", plot=ps[[2]],
       width=12, height=3, units="in", dpi=600)
ggsave(filename="~/NPI/code_output/figs/Fig5C.png", plot=ps[[3]],
       width=12, height=3, units="in", dpi=600)



