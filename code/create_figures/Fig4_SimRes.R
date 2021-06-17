### Requires tidyverse package and NPI_SS_Formulae.R commands
require(tidyverse)
require(gridExtra)
source(file="~/NPI/code/generate_results/NPI_SS_Formulae.R")

### Load results of simulation sample size calculations
lb_ttest <- read.csv(file="~/NPI/code_output/tables/ttest_rerun.csv", stringsAsFactors = F)

lb_ttest <- as_tibble(lb_ttest) %>% mutate(sufficient_ncluster=ifelse(sufficient_ncluster==">=1000",NA,as.numeric(sufficient_ncluster)),
                                           adjusted_sufficient_ncluster=ifelse(adjusted_sufficient_ncluster==">=1000",NA,
                                                                               as.numeric(adjusted_sufficient_ncluster)),
                                           perc_sample = nsample/N) %>%
  mutate(indivs_per_arm = nsample*sufficient_ncluster)

### Fig. 4: Simulation Results:
plot.4A <- ggplot(lb_ttest) + 
  geom_line(mapping=aes(x=k, y=sufficient_ncluster, color=effect, linetype=factor(R0), group=paste0(nsample,R0,effect,E_It))) + 
  facet_wrap(facets=vars(factor(N, levels=c(100,1000,10000),
                                labels=c("N = 100","N = 1,000","N = 10,000"))), 
             nrow=3, ncol=1, scales="free_x") +
  scale_x_continuous(name="k\n", breaks=c(0.1,0.4,0.7), minor_breaks=NULL) +
  scale_color_continuous("Effect Size",
                         breaks=c(.2,.4),
                         labels=c("20%","40%")) +
  scale_linetype_discrete(expression(R[0]^0),
                          breaks=c(1.5,2),
                          labels=c("1.5","2.0")) +
  theme_bw() + theme(text=element_text(size=30),
                   axis.title.y = element_text(size=30),
                   axis.title.x = element_text(size=30),
                   title=element_text(size=20),
                   plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
                   legend.position="none")+
  labs(y="Number of Clusters Per Arm",
       title="A) Required Clusters by k\n")

plot.4B <- ggplot(lb_ttest) + 
  geom_line(mapping=aes(x=factor(paste0(perc_sample*100,"%"), levels=c("1%","10%","50%","100%")),
                        y=sufficient_ncluster, color=effect, linetype=factor(R0), group=paste0(R0,effect,k))) + 
  facet_wrap(facets=vars(factor(N, levels=c(100,1000,10000),
                                labels=c("N = 100","N = 1,000","N = 10,000"))), 
             nrow=3, ncol=1, scales="free_x") +
  scale_x_discrete("Percent of Individuals \nSampled Per Cluster",
                   breaks=c("1%","10%","50%","100%")) +
  scale_color_continuous("Effect Size",
                         breaks=c(.2,.4),
                         labels=c("20%","40%")) +
  scale_linetype_discrete(expression(R[0]^0),
                          breaks=c(1.5,2),
                          labels=c("1.5","2.0")) +
  theme_bw() + theme(text=element_text(size=30),
                     axis.title.y = element_text(size=30),
                     axis.title.x = element_text(size=30),
                     title=element_text(size=20),
                     plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
                     legend.position="none")+
  labs(y="Number of Clusters Per Arm",
       title="B) Required Clusters by \nPercent Sampled")

plot.4C <- ggplot(lb_ttest) + 
  geom_line(mapping=aes(x=factor(paste0(perc_sample*100,"%"), levels=c("1%","10%","50%","100%")),
                        y=indivs_per_arm/1000, color=effect, linetype=factor(R0), group=paste0(R0,effect,k))) + 
  facet_wrap(facets=vars(factor(N, levels=c(100,1000,10000),
                                labels=c("N = 100","N = 1,000","N = 10,000"))), 
             nrow=3, ncol=1, scales="free") +
  scale_x_discrete("Percent of Individuals \nSampled Per Cluster",
                   breaks=c("1%","10%","50%","100%")) +
  scale_color_continuous("Effect Size",
                         breaks=c(.2,.4),
                         labels=c("20%","40%")) +
  scale_linetype_discrete(expression(R[0]^0),
                          breaks=c(1.5,2),
                          labels=c("1.5","2.0")) +
  theme_bw() + theme(text=element_text(size=30),
                     axis.title.y = element_text(size=30),
                     axis.title.x = element_text(size=30),
                     title=element_text(size=20),
                     plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
                     legend.position="none")+
  labs(y="Number of Clusters Per Arm",
       title="C) Required Individuals by \nPercent Sampled")

l <- list(plot.4A, ggplot() + theme_void(), plot.4B, ggplot() + theme_void(), plot.4C)
ggsave(filename="~/NPI/code_output/figs/Fig4.tiff", marrangeGrob(grobs = l, nrow=2, ncol=3, top=NULL),
       width=16.5, height=24, units='in', dpi=600)
