### Requires tidyverse package and NPI_SS_Formulae.R commands
require(tidyverse)
require(gridExtra)
source(file="~/NPI/code/generate_results/NPI_SS_Formulae.R")

### Figure 2: Full Sampling Formula Results
R_full <- expand_grid(Rt0=seq(1,4,by=.01),k=seq(.1,1,by=.3)) %>%
  cbind(tibble(effect=0.4,popsize=1000,EIt=0.005))

R_full$full <- apply(R_full, MARGIN=1,
                     FUN=function(vals) SScalc.full(Rt1=vals["Rt0"]*(1-vals["effect"]), Rt0=vals["Rt0"], 
                                                    k1=vals["k"], k0=vals["k"], 
                                                    n=vals["popsize"],
                                                    EIt=.005,
                                                    pow=.8, alpha=.05, N=NULL))

plot_r_k_full <- ggplot(R_full, mapping=aes(x=Rt0, y=full, group=k, color=k)) + geom_line() +
  labs(title=expression("A) Required Clusters vs. "*R[t]^0*" and "*k),
       x=expression(R[t]^0), y="Number of Clusters Per Arm",
       color="k") +
  scale_color_continuous(breaks=unique(R_full$k)) +
  theme_bw() + theme(text=element_text(size=12),
                     legend.position="bottom") +
  coord_cartesian(xlim=c(1,4.1), ylim=c(0,160), expand=FALSE, clip="on")

ggsave(filename="~/NPI/code_output/figs/Full_R_k.png", plot=plot_r_k_full,
       width=6, height=6, units="in", dpi=600)

R_full_EIt <- expand_grid(EIt=seq(.001,.05, by=.001),popsize=seq(100,1000,by=300)) %>%
  cbind(tibble(effect=0.4,Rt0=1.5,k=0.4))

R_full_EIt$full <- apply(R_full_EIt, MARGIN=1,
                         FUN=function(vals) SScalc.full(Rt1=vals["Rt0"]*(1-vals["effect"]), Rt0=vals["Rt0"], 
                                                        k1=vals["k"], k0=vals["k"], 
                                                        n=vals["popsize"],
                                                        EIt=vals["EIt"],
                                                        pow=.8, alpha=.05, N=NULL))

plot_EIt_full <- ggplot(R_full_EIt, mapping=aes(x=EIt, y=full, group=popsize, color=popsize)) + geom_line() +
  labs(title=expression("B) Required Clusters vs. E["*I[t]*"] and n"),
       x=expression("E["*I[t]*"]"), y="Number of Clusters Per Arm",
       color="Cluster Size") +
  scale_color_continuous(breaks=unique(R_full_EIt$popsize)) +
  coord_cartesian(ylim=c(0,300), xlim=c(0,0.051), expand=FALSE, clip="on") +
  theme_bw() + theme(text=element_text(size=12),
                     legend.position="bottom")

ggsave(filename="~/NPI/code_output/figs/Full_EIt.png", plot=plot_EIt_full,
       width=6, height=6, units="in", dpi=600)


### Figure 3: Approximations Accounting for Sampling:
Baseline.vals <- tibble(Rt0=1.5, effect=0.4, k=0.4, popsize=10000,
                        sampsize=100, EIt=0.005)

##### Fig 3A: R_t and k
R_sample <- expand_grid(Rt0=seq(1,4,by=.01),k=seq(.1,1,by=.3)) %>%
  cbind(Baseline.vals %>% dplyr::select(!c(Rt0,k)))
R_sample$samp <- apply(R_sample, MARGIN=1,
                     FUN=function(vals) SScalc.samp(Rt1=vals["Rt0"]*(1-vals["effect"]), Rt0=vals["Rt0"],
                                                    k1=vals["k"], k0=vals["k"],
                                                    n=vals["popsize"], m=vals["sampsize"],
                                                    EIt=.005,
                                                    pow=.8, alpha=.05, N=NULL))

plot_r_k <- ggplot(R_sample, mapping=aes(x=Rt0, y=samp, group=k, color=k)) + geom_line() +
  labs(title=expression("A) Required Clusters vs. "*R[t]^0*" and "*k),
       x=expression(R[t]^0), y="Number of Clusters Per Arm",
       color="k") +
  scale_color_continuous(breaks=unique(R_sample$k)) +
  theme_bw() + theme(text=element_text(size=12),
                     legend.position="bottom") +
  coord_cartesian(ylim=c(0,200), xlim=c(1,4.1), expand=FALSE, clip="on") #+
  # geom_point(aes(x=Rt0, y=suffNo, group=k, color=k))

ggsave(filename="~/NPI/code_output/figs/R_k.png", plot=plot_r_k,
       width=6, height=6, units="in", dpi=600)


##### Fig 3B: R_t and Effect Size
R_delta_sample <- expand_grid(Rt0=seq(1,4,by=.01), effect=seq(.2,.9,by=.1)) %>%
  cbind(Baseline.vals %>% dplyr::select(!c(Rt0,effect)))

R_delta_sample$samp <- apply(R_delta_sample, MARGIN=1,
                             FUN=function(vals) SScalc.samp(Rt1=vals["Rt0"]*(1-vals["effect"]), Rt0=vals["Rt0"],
                                                            k1=vals["k"], k0=vals["k"],
                                                            n=vals["popsize"], m=vals["sampsize"],
                                                            EIt=.005,
                                                            pow=.8, alpha=.05, N=NULL))

plot_r_delta <- ggplot(R_delta_sample, mapping=aes(x=Rt0, y=samp, group=effect, color=effect)) + geom_line() +
  labs(title=expression("B) Required Clusters vs. "*R[t]^0*" and Effect Size"),
       x=expression(R[t]^0), y="Number of Clusters Per Arm",
       color="Effect Size") +
  scale_color_continuous(breaks=seq(.1,.9,by=.2), labels=paste0(seq(10,90,by=20),"%")) +
  theme_bw() + theme(text=element_text(size=12),
                     legend.position="bottom") +
  coord_cartesian(ylim=c(0,800), xlim=c(1,4.1), expand=FALSE, clip="on")

ggsave(filename="~/NPI/code_output/figs/R_Delta.png", plot=plot_r_delta,
       width=6, height=6, units="in", dpi=600)


##### Fig 3C: Sample Size and Cluster Size
V_sample <- expand_grid(sampsize=20:250, popsize=seq(1000,10000,by=1500)) %>% 
  cbind(Baseline.vals %>% dplyr::select(!c(sampsize,popsize)))

V_sample$samp <- apply(V_sample, MARGIN=1,
                       FUN=function(vals) SScalc.samp(Rt1=vals["Rt0"]*(1-vals["effect"]), Rt0=vals["Rt0"], 
                                                  k1=vals["k"], k0=vals["k"], 
                                                  n=vals["popsize"], m=vals["sampsize"],
                                                  EIt=.005,
                                                  pow=.8, alpha=.05, N=NULL))

plot_sizes <- ggplot(V_sample, mapping=aes(x=sampsize, y=samp, group=popsize, color=popsize)) + geom_line() +
  labs(title="C) Required Clusters vs. Number Sampled Per Cluster and Cluster Size",
       x="Number Sampled Per Cluster", y="Number of Clusters Per Arm",
       color="Cluster Size") +
  scale_color_continuous(breaks=unique(V_sample$popsize)[c(1,3,5,7)]) +
  coord_cartesian(ylim=c(0,400), xlim=c(0,255), expand=FALSE, clip="on")+
  theme_bw() + theme(text=element_text(size=12),
                     legend.position="bottom")

ggsave(filename="~/NPI/code_output/figs/Sizes.png", plot=plot_sizes,
       width=6, height=6, units="in", dpi=600)


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
                   legend.position="none")+#,
  #legend.position="bottom") +
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
                     legend.position="none")+#,
  #legend.position="bottom") +
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
                     legend.position="none")+#,
  labs(y="Number of Clusters Per Arm",
       title="C) Required Individuals by \nPercent Sampled")

#ggsave(filename="~/NPI/code_output/figs/Fig4A.png", plot=plot.4A,
#       width=5.5, height=12, units="in", dpi=600)
#ggsave(filename="~/NPI/code_output/figs/Fig4B.png", plot=plot.4B,
#       width=5.5, height=12, units="in", dpi=600)
#ggsave(filename="~/NPI/code_output/figs/Fig4C.png", plot=plot.4C,
#       width=5.5, height=12, units="in", dpi=600)
l <- list(plot.4A, ggplot() + theme_void(), plot.4B, ggplot() + theme_void(), plot.4C)
ggsave(filename="~/NPI/code_output/figs/Fig4.tiff", marrangeGrob(grobs = l, nrow=2, ncol=3, top=NULL),
       width=16.5, height=24, units='in', dpi=600)
