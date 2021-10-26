# -----------------------------------------------------------------------------
# @description: script used to create the figure showing effect of lag on sample size
# -----------------------------------------------------------------------------

require(reshape2)
require(ggplot2)
require(gridExtra)
require(grid)
lb <- read.csv("~/NPI/code_output/tables/ttest.csv", stringsAsFactors = F)
lb2 <- read.csv("~/NPI/code_output/tables/ttest_gen2.csv", stringsAsFactors = F)
lb3 <- read.csv("~/NPI/code_output/tables/ttest_gen3.csv", stringsAsFactors = F)
param_sets <- list()
param_sets_dex <- 1
R0s <- c(1.5, 2)
ns <- c(1000, 10000)
ks <- c(0.1, 0.4, 0.7)
effects <- c(0.2, 0.4)
nsamples <- c(100, 1000)
for (R0 in R0s) {
  for (n in ns) {
    for (k in ks) {
      for (effect in effects) {
        for (nsample in nsamples) {
          if (!(k == 0.1 & n == 10000 & R0 == 1.5)) {
            param_sets[[param_sets_dex]] <- c(R0, n, k, effect, nsample)
            param_sets_dex <- param_sets_dex + 1
          }
        }
      }
    }
  }
}
ns <- c(100)
nsamples <- c(10, 50, 100)
for (R0 in R0s) {
  for (n in ns) {
    for (k in ks) {
      for (effect in effects) {
        for (nsample in nsamples) {
          if (!(k == 0.1 & n == 10000 & R0 == 1.5)) {
            param_sets[[param_sets_dex]] <- c(R0, n, k, effect, nsample)
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
  colnames(new_row) <- c("one", "two", "three", "effect", "R0", "n", "k")
  one <- lb[which(lb$R0 == param_set[1] & lb$N == param_set[2] & lb$k == param_set[3] & lb$effect == param_set[4] & lb$nsample == param_set[5]),]
  new_row$one <- as.numeric(ifelse(one$sufficient_ncluster == ">=1000" | one$sufficient_ncluster == "-1", 1000, one$sufficient_ncluster))
  two <- lb2[which(lb2$R0 == param_set[1] & lb2$N == param_set[2] & lb2$k == param_set[3] & lb2$effect == param_set[4] & lb2$nsample == param_set[5]),]
  new_row$two <- as.numeric(ifelse(two$sufficient_ncluster == ">=1000" | two$sufficient_ncluster == "-1", 1000, two$sufficient_ncluster))
  three <- lb3[which(lb3$R0 == param_set[1] & lb3$N == param_set[2] & lb3$k == param_set[3] & lb3$effect == param_set[4] & lb3$nsample == param_set[5]),]
  new_row$three <- as.numeric(ifelse(three$sufficient_ncluster == ">=1000" | three$sufficient_ncluster == "-1", 1000, three$sufficient_ncluster))
  new_row$effect <- param_set[4]
  new_row$R0 <- param_set[1]
  new_row$n <- param_set[2]
  new_row$k <- param_set[3]
  toplot_ls[[toplot_ls_dex]] <- new_row
  toplot_ls_dex <- toplot_ls_dex + 1
}
toplot <- do.call(rbind, toplot_ls)
toplot$one <- ifelse(toplot$one == 1000, NA, toplot$one)
toplot$two <- ifelse(toplot$two == 1000, NA, toplot$two)
toplot$three <- ifelse(toplot$three == 1000, NA, toplot$three)
ns <- c(100, 1000, 10000)
ps <- list()
ps_dex <- 1
for (n in ns) {
  print(n)
  helper_function <- function(toplot, k, n) {
    subplot <- toplot[which(toplot$n == n & toplot$k == k),]
    subplot_20_1.5 <- subplot[which(subplot$effect == "0.2" & subplot$R0 == "1.5"),]
    subplot_20_2 <- subplot[which(subplot$effect == "0.2" & subplot$R0 == "2"),]
    subplot_40_1.5 <- subplot[which(subplot$effect == "0.4" & subplot$R0 == "1.5"),]
    subplot_40_2 <- subplot[which(subplot$effect == "0.4" & subplot$R0 == "2"),]
    subplot_20_1.5 <- data.frame(t(subplot_20_1.5))[1:3,]
    subplot_20_1.5$x <- c('1st', '2nd', '3rd')
    subplot_20_1.5 <- melt(subplot_20_1.5, 'x')
    subplot_20_2 <- data.frame(t(subplot_20_2))[1:3,]
    subplot_20_2$x <- c('1st', '2nd', '3rd')
    subplot_20_2 <- melt(subplot_20_2, 'x')
    subplot_40_1.5 <- data.frame(t(subplot_40_1.5))[1:3,]
    subplot_40_1.5$x <- c('1st', '2nd', '3rd')
    subplot_40_1.5 <- melt(subplot_40_1.5, 'x')
    subplot_40_2 <- data.frame(t(subplot_40_2))[1:3,]
    subplot_40_2$x <- c('1st', '2nd', '3rd')
    subplot_40_2 <- melt(subplot_40_2, 'x')
    subplot_20_2$title <- paste0('n = ', n, ", k = ", k)
    if (n == 10000 & k == 0.1) {
      subplot_final <- ggplot() + 
        geom_line(size=1.5, subplot_20_2, mapping=aes(x=x,y=value,group=variable,color='20'), linetype='solid')+
        geom_line(size=1.5, subplot_40_2, mapping=aes(x=x,y=value,group=variable,color='40'), linetype='solid')+
        ylab('')+
        xlab('')+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.position='none',
              text = element_text(size=30),
              axis.title.x = element_text(size=23),
              axis.title.y = element_text(size=21))+
        facet_grid(. ~ title)+
        scale_color_manual(values=c('steelblue1', 'steelblue4'))
    } else if (n == 1000 & k == 0.1) {
      subplot_final <- ggplot() + 
        geom_line(size=1.5, subplot_20_1.5, mapping=aes(x=x,y=value,group=variable,color='20'), linetype='dashed')+
        geom_line(size=1.5, subplot_20_2, mapping=aes(x=x,y=value,group=variable,color='20'), linetype='solid')+
        geom_line(size=1.5, subplot_40_1.5, mapping=aes(x=x,y=value,group=variable,color='40'), linetype='dashed')+
        geom_line(size=1.5, subplot_40_2, mapping=aes(x=x,y=value,group=variable,color='40'), linetype='solid')+
        ylab('Number of Clusters Per Arm')+
        xlab('')+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.position='none',
              text = element_text(size=30),
              axis.title.x = element_text(size=23),
              axis.title.y = element_text(size=21))+
        facet_grid(. ~ title)+
        scale_color_manual(values=c('steelblue1', 'steelblue4'))
    } else if (n == 10000 & k == 0.4) {
      subplot_final <- ggplot() + 
        geom_line(size=1.5, subplot_20_1.5, mapping=aes(x=x,y=value,group=variable,color='20'), linetype='dashed')+
        geom_line(size=1.5, subplot_20_2, mapping=aes(x=x,y=value,group=variable,color='20'), linetype='solid')+
        geom_line(size=1.5, subplot_40_1.5, mapping=aes(x=x,y=value,group=variable,color='40'), linetype='dashed')+
        geom_line(size=1.5, subplot_40_2, mapping=aes(x=x,y=value,group=variable,color='40'), linetype='solid')+
        ylab('')+
        xlab('Generation of Sampling')+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.position='none',
              text = element_text(size=30),
              axis.title.x = element_text(size=23),
              axis.title.y = element_text(size=21))+
        facet_grid(. ~ title)+
        scale_color_manual(values=c('steelblue1', 'steelblue4'))
    } else {
      subplot_final <- ggplot() + 
        geom_line(size=1.5, subplot_20_1.5, mapping=aes(x=x,y=value,group=variable,color='20'), linetype='dashed')+
        geom_line(size=1.5, subplot_20_2, mapping=aes(x=x,y=value,group=variable,color='20'), linetype='solid')+
        geom_line(size=1.5, subplot_40_1.5, mapping=aes(x=x,y=value,group=variable,color='40'), linetype='dashed')+
        geom_line(size=1.5, subplot_40_2, mapping=aes(x=x,y=value,group=variable,color='40'), linetype='solid')+
        ylab('')+
        xlab('')+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.position='none',
              text = element_text(size=30),
              axis.title.x = element_text(size=23),
              axis.title.y = element_text(size=21))+
        facet_grid(. ~ title)+
        scale_color_manual(values=c('steelblue1', 'steelblue4'))
    }

    return(subplot_final)
  }
  subplot_1 <- helper_function(toplot, 0.1, n)
  subplot_4 <- helper_function(toplot, 0.4, n)
  subplot_7 <- helper_function(toplot, 0.7, n)
  if (n == 100) {
    p <- grid.arrange(arrangeGrob(subplot_1,
                                  subplot_4,
                                  subplot_7,
                                  nrow=1),
                      nrow=1,heights=c(10))
  } else if (n == 1000) {
    p <- grid.arrange(arrangeGrob(subplot_1,
                                  subplot_4,
                                  subplot_7,
                                  nrow=1),
                      nrow=1,heights=c(10))
  } else {
    p <- grid.arrange(arrangeGrob(subplot_1,
                                  subplot_4,
                                  subplot_7,
                                  nrow=1),
                      nrow=1, heights=c(10))
  }
  ps[[ps_dex]] <- p
  ps_dex <- ps_dex + 1
}
ps[[4]] <- ggplot() + theme_void()
ggsave(filename="~/NPI/code_output/figs/Sheen_Fig5.tiff", marrangeGrob(grobs = ps, nrow=4, ncol=1, top=NULL),
       width=13, height=16, units='in', dpi=600)
ggsave(filename="~/NPI/code_output/figs/Sheen_Fig5.png", marrangeGrob(grobs = ps, nrow=4, ncol=1, top=NULL),
       width=13, height=16, units='in', dpi=72)

### Create legend
legend_plot <- ggplot(lb_ttest) + 
  geom_line(size=1.5, mapping=aes(x=factor(paste0(perc_sample*100,"%"), levels=c("1%","10%","50%","100%")),
                                  y=indivs_per_arm/1000, color=effect, linetype=factor(R0), group=paste0(R0,effect,k))) + 
  facet_wrap(facets=vars(factor(n, levels=c(100,1000,10000),
                                labels=c("n = 100","n = 1,000","n = 10,000"))), 
             nrow=3, ncol=1, scales="free") +
  scale_x_discrete("Percent of Individuals \nSampled Per Cluster",
                   breaks=c("1%","10%","50%","100%")) +
  scale_color_manual("Effect Size", values = c("0.2" = "steelblue1", "0.4" = "steelblue4")) +
  #scale_color_continuous("Effect Size",
  #                       breaks=c(.2,.4),
  #                       labels=c("20%","40%")) +
  scale_linetype_discrete(expression(R[0]^C),
                          breaks=c(1.5,2),
                          labels=c("1.5","2.0")) +
  theme_bw() + theme(text=element_text(size=30),
                     axis.title.y = element_text(size=30),
                     axis.title.x = element_text(size=30),
                     legend.title = element_text(size = 30), 
                     legend.text = element_text(size = 25),
                     title=element_text(size=20),
                     plot.margin=unit(c(0,0,0,0),"cm"),
                     legend.position="bottom")+
  labs(y="Number of Individuals Per Arm / 1000",
       title="C) Required Individuals by \nPercent Sampled") +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(legend.key.size = grid::unit(5, "lines"))
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(legend_plot)
#legend_grid <- grid.draw(mylegend)
ggsave(filename="~/NPI/code_output/figs/Sheen_Fig5_legend.tiff", mylegend,
       width=12, height=1, units='in', dpi=600)
ggsave(filename="~/NPI/code_output/figs/Sheen_Fig5_legend.png", mylegend,
       width=12, height=1, units='in', dpi=72)

