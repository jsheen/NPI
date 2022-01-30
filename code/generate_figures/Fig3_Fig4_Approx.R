### Requires tidyverse package and NPI_SS_Formulae.R commands
require(tidyverse)
require(gridExtra)
source(file="~/NPI/code/NPI_SS_Formulae.R")


### Figure 3: Full Sampling Formula Results
R_full <- expand_grid(Rt0=seq(1,4,by=.01),k=seq(.1,1,by=.3)) %>%
  cbind(tibble(effect=0.4,popsize=1000,EIt=0.005))

R_full$full <- apply(R_full, MARGIN=1,
                     FUN=function(vals) SScalc.full(Rt1=vals["Rt0"]*(1-vals["effect"]), Rt0=vals["Rt0"], 
                                                    k1=vals["k"], k0=vals["k"], 
                                                    n=vals["popsize"],
                                                    EIt=.005,
                                                    pow=.8, alpha=.05, N=NULL))
R_full$k <- factor(R_full$k)
plot_r_k_full <- ggplot(R_full, mapping=aes(x=Rt0, y=full, group=k, color=k)) + geom_line(size=1.5) +
  labs(title=expression("A) Required Clusters vs. "*R[t]^C*" and "*k),
       x=expression(R[t]^C), y="Number of Clusters Per Arm",
       color="k") +
  #scale_color_continuous(breaks=unique(R_full$k)) +
  theme_bw() + theme(text=element_text(size=25),
                     title=element_text(size=24),
                     axis.title=element_text(size=30),
                     legend.position="bottom",
                     legend.title = element_text(size = 35), 
                     legend.text = element_text(size = 30),
                     axis.text.x = element_text(size=30),
                     axis.text.y = element_text(size=30),
                     plot.margin=unit(c(0.1,2,0.1,0.1),"cm"),
                     legend.key.size = unit(1, 'cm')) +
  scale_color_manual(values = c("0.1" = "steelblue1", "0.4" = "steelblue2",
                                "0.7" = "steelblue3", "1" = "steelblue4")) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  coord_cartesian(xlim=c(1,4.1), ylim=c(0,175), expand=FALSE, clip="on")

R_full_EIt <- expand_grid(EIt=seq(.001,.05, by=.001),popsize=seq(100,1000,by=300)) %>%
  cbind(tibble(effect=0.4,Rt0=1.5,k=0.4))

R_full_EIt$full <- apply(R_full_EIt, MARGIN=1,
                         FUN=function(vals) SScalc.full(Rt1=vals["Rt0"]*(1-vals["effect"]), Rt0=vals["Rt0"], 
                                                        k1=vals["k"], k0=vals["k"], 
                                                        n=vals["popsize"],
                                                        EIt=vals["EIt"],
                                                        pow=.8, alpha=.05, N=NULL))
R_full_EIt$popsize <- factor(R_full_EIt$popsize)
plot_EIt_full <- ggplot(R_full_EIt, mapping=aes(x=EIt, y=full, group=popsize, color=popsize)) + geom_line(size=1.5) +
  labs(title=expression("B) Required Clusters vs. E["*I[t]*"] and n"),
       x=expression("E["*I[t]*"]"), y="",
       color="Cluster Size") +
  #scale_color_continuous(breaks=unique(R_full_EIt$popsize)) +
  coord_cartesian(ylim=c(0,800), xlim=c(0,0.051), expand=FALSE, clip="on") +
  theme_bw() + theme(text=element_text(size=25),
                     title=element_text(size=24),
                     axis.title=element_text(size=30),
                     legend.position="bottom",
                     legend.title = element_text(size = 28), 
                     legend.text = element_text(size = 24),
                     axis.text.x = element_text(size=30),
                     axis.text.y = element_text(size=30),
                     plot.margin=unit(c(0.1,2.2,0.1,0.1),"cm"),
                     legend.key.size = unit(1, 'cm')) +
                     scale_color_manual(values = c("100" = "steelblue1", "400" = "steelblue2",
                                                   "700" = "steelblue3", "1000" = "steelblue4")) +
  guides(colour = guide_legend(override.aes = list(size=10)))

l <- list(plot_r_k_full, plot_EIt_full)
ggsave(filename="~/NPI/code_output/figs/Sheen_Fig3.tiff", marrangeGrob(grobs = l, nrow=1, ncol=2, top=NULL),
       width=16, height=8, units='in', dpi=600)
ggsave(filename="~/NPI/code_output/figs/Sheen_Fig3.eps", marrangeGrob(grobs = l, nrow=1, ncol=2, top=NULL),
       width=16, height=8, units='in', dpi=600)
ggsave(filename="~/NPI/code_output/figs/Sheen_Fig3.png", marrangeGrob(grobs = l, nrow=1, ncol=2, top=NULL),
       width=16, height=8, units='in', dpi=72)

### Figure 4: Approximations Accounting for Sampling:
Baseline.vals <- tibble(Rt0=1.5, effect=0.4, k=0.4, popsize=10000,
                        sampsize=100, EIt=0.005)

##### Fig 4A: R_t and k
R_sample <- expand_grid(Rt0=seq(1,4,by=.01),k=seq(.1,1,by=.3))
R_sample$effect <- 0.4
R_sample$popsize <- 10000
R_sample$sampsize <- 100
R_sample$samp <- apply(R_sample, MARGIN=1,
                       FUN=function(vals) SScalc.samp(Rt1=vals["Rt0"]*(1-vals["effect"]), Rt0=vals["Rt0"],
                                                      k1=vals["k"], k0=vals["k"],
                                                      n=vals["popsize"], m=vals["sampsize"],
                                                      EIt=.005,
                                                      pow=.8, alpha=.05, N=NULL))
R_sample$k <- factor(R_sample$k)
plot_r_k <- ggplot(R_sample, mapping=aes(x=Rt0, y=samp, group=k, color=k)) + geom_line(size=1.5) +
  labs(title=expression(atop("A) Required Clusters vs.", R[t]^C*" and "*k)),
       x=expression(R[t]^C), y="Number of Clusters Per Arm",
       color="k") +
  #scale_color_continuous(breaks=unique(R_sample$k)) +
  theme_bw() + theme(text=element_text(size=25),
                     title=element_text(size=24),
                     axis.title=element_text(size=30),
                     legend.position="bottom",
                     legend.title = element_text(size = 35), 
                     legend.text = element_text(size = 30),
                     axis.text.x = element_text(size=30),
                     axis.text.y = element_text(size=30),
                     plot.margin=unit(c(0.1,3,0.1,0.1),"cm"),
                     legend.key.size = unit(1, 'cm')) +
  coord_cartesian(ylim=c(0,200), xlim=c(1,4.1), expand=FALSE, clip="on") +
  scale_color_manual(values = c("0.1" = "steelblue1", "0.4" = "steelblue2",
                                "0.7" = "steelblue3", "1" = "steelblue4")) +
  guides(colour = guide_legend(override.aes = list(size=10)))

##### Fig 4B: R_t and Effect Size
R_delta_sample <- expand_grid(Rt0=seq(1,4,by=.01), effect=seq(.2,.8,by=.2)) 
R_delta_sample$k <- 0.4
R_delta_sample$popsize <- 10000
R_delta_sample$sampsize <- 100
R_delta_sample$samp <- apply(R_delta_sample, MARGIN=1,
                             FUN=function(vals) SScalc.samp(Rt1=vals["Rt0"]*(1-vals["effect"]), Rt0=vals["Rt0"],
                                                            k1=vals["k"], k0=vals["k"],
                                                            n=vals["popsize"], m=vals["sampsize"],
                                                            EIt=.005,
                                                            pow=.8, alpha=.05, N=NULL))
R_delta_sample$effect <- factor(R_delta_sample$effect)
plot_r_delta <- ggplot(R_delta_sample, mapping=aes(x=Rt0, y=samp, group=effect, color=effect)) + geom_line(size=1.5) +
  labs(title=expression(atop("B) Required Clusters vs.", R[t]^C*" and Effect Size")),
       x=expression(R[t]^C), y="",
       color="Effect Size") +
  #scale_color_continuous(breaks=seq(.1,.9,by=.2), labels=paste0(seq(10,90,by=20),"%")) +
  theme_bw() + theme(text=element_text(size=25),
                     title=element_text(size=24),
                     axis.title=element_text(size=30),
                     legend.position="bottom",
                     legend.title = element_text(size = 30), 
                     legend.text = element_text(size = 28),
                     axis.text.x = element_text(size=30),
                     axis.text.y = element_text(size=30),
                     plot.margin=unit(c(0.1,3,0.1,0.1),"cm"),
                     legend.key.size = unit(1, 'cm')) +
  coord_cartesian(ylim=c(0,800), xlim=c(1,4.1), expand=FALSE, clip="on") +
  scale_color_manual(values = c("0.2" = "steelblue1", "0.4" = "steelblue2",
                                "0.6" = "steelblue3", "0.8" = "steelblue4")) +
  guides(colour = guide_legend(override.aes = list(size=10)))

##### Fig 4C: Sample Size and Cluster Size
V_sample <- expand_grid(sampsize=10:250, popsize=seq(1000,10000,by=3000)) 
V_sample$effect <- 0.4
V_sample$k <- 0.4
V_sample$Rt0 <- 1.5
V_sample$samp <- apply(V_sample, MARGIN=1,
                       FUN=function(vals) SScalc.samp(Rt1=vals["Rt0"]*(1-vals["effect"]), Rt0=vals["Rt0"], 
                                                      k1=vals["k"], k0=vals["k"], 
                                                      n=vals["popsize"], m=vals["sampsize"],
                                                      EIt=.005,
                                                      pow=.8, alpha=.05, N=NULL))
V_sample$popsize <- factor(V_sample$popsize)
plot_sizes <- ggplot(V_sample, mapping=aes(x=sampsize, y=samp, group=popsize, color=popsize)) + geom_line(size=1.5) +
  labs(title="C) Required Clusters vs. Number \nSampled Per Cluster and Cluster Size",
       x="Number Sampled Per Cluster", y="",
       color="Cluster Size") +
  #scale_color_continuous(breaks=unique(V_sample$popsize)[c(1,3,5,7)]) +
  coord_cartesian(ylim=c(0,800), xlim=c(0,255), expand=FALSE, clip="on")+
  theme_bw() + theme(text=element_text(size=25),
                     title=element_text(size=24),
                     axis.title=element_text(size=30),
                     legend.position="bottom",
                     legend.title = element_text(size = 28), 
                     legend.text = element_text(size = 20),
                     axis.text.x = element_text(size=30),
                     axis.text.y = element_text(size=30),
                     plot.margin=unit(c(0.1,3,0.1,0.1),"cm"),
                     legend.key.size = unit(1, 'cm')) +
  scale_color_manual(values = c("1000" = "steelblue1", "4000" = "steelblue2",
                                "7000" = "steelblue3", "10000" = "steelblue4")) +
  guides(colour = guide_legend(override.aes = list(size=10)))

l <- list(plot_r_k, plot_r_delta, plot_sizes)
ggsave(filename="~/NPI/code_output/figs/Sheen_Fig4.tiff", marrangeGrob(grobs = l, nrow=1, ncol=3, top=NULL),
       width=24, height=8, units='in', dpi=600)
ggsave(filename="~/NPI/code_output/figs/Sheen_Fig4.eps", marrangeGrob(grobs = l, nrow=1, ncol=3, top=NULL),
       width=24, height=8, units='in', dpi=600)
ggsave(filename="~/NPI/code_output/figs/Sheen_Fig4.png", marrangeGrob(grobs = l, nrow=1, ncol=3, top=NULL),
       width=24, height=8, units='in', dpi=72)

