# Import libraries and set input and output folders ----------------------------
library(ggplot2)
input_folder <- "/Users/jsheen/SW-CRT-outbreak/NPI_study/EoN/code_output/csvs_1_2/It1_It/dists/"
output_folder <- "/Users/jsheen/SW-CRT-outbreak/NPI_study/EoN/code_output/csvs_1_2/It1_It/pdfs/"

beta_vals <- c(0.02, 0.03, 0.04, 0.05)
ncomms <- c(20, 40, 60, 80, 100) / 2
effects <- c(0, 0.2, 0.4, 0.6)
days <- c(14, 21, 28)
gens <- 3

plot_ls_It <- list()
plot_ls_It1_It <- list()
plot_ls_log_It1_It <- list()
plot_ls_It_dex <- 1
plot_ls_It1_It_dex <- 1
plot_ls_log_It1_It_dex <- 1
for (day in days) {
  for (beta in beta_vals) {
    if (beta == 0.02) {
      R0 <- 1.34
    } else if (beta == 0.03) {
      R0 <- 1.93
    } else if (beta == 0.04) {
      R0 <- 2.47
    } else {
      R0 <- 3.00
    }
    for (ncomm in ncomms) {
      for (effect in effects) {
        file_name_It <- paste0(input_folder, "beta_", beta, "/effect_", effect, "/intervention_", day, "/ncomm_", ncomm, ".0/It.csv")
        file_name_It1 <- paste0(input_folder, "beta_", beta, "/effect_", effect, "/intervention_", day, "/ncomm_", ncomm, ".0/It1.csv")
        file_name_It1_It <- paste0(input_folder, "beta_", beta, "/effect_", effect, "/intervention_", day, "/ncomm_", ncomm, ".0/It1_It.csv")
        file_name_log_It1_It <- paste0(input_folder, "beta_", beta, "/effect_", effect, "/intervention_", day, "/ncomm_", ncomm, ".0/log_It1_It.csv")
        It <- data.frame(t(read.csv(file_name_It, header=F)[,1:300])) / (ncomm * 500)
        It1 <- data.frame(t(read.csv(file_name_It1, header=F)[,1:300])) / (ncomm * 500)
        It1_It <- data.frame(t(read.csv(file_name_It1_It, header=F)[,1:300]))
        log_It1_It <- data.frame(t(read.csv(file_name_log_It1_It, header=F)[,1:300]))
        plot_ls_It[[plot_ls_It_dex]] <- ggplot(It, aes(x=X1)) + geom_histogram(binwidth=0.001) + xlab("I_t / N") + ggtitle(paste0('Day: ', day, '; ',
                                                                                                                'R0: ', R0, '\n',
                                                                                                                'NComm: ', ncomm, '; ',
                                                                                                                'Effect: ', effect)) + ylab('') + xlim(0, 0.15)
        plot_ls_It_dex <- plot_ls_It_dex + 1
        plot_ls_It[[plot_ls_It_dex]] <- ggplot(It1, aes(x=X1)) + geom_histogram(binwidth=0.001) + xlab('I_(t+1) / N') + ggtitle('\n') + ylab('') + xlim(0, 0.15)
        plot_ls_It_dex <- plot_ls_It_dex + 1
        plot_ls_It[[plot_ls_It_dex]] <- ggplot(It1, aes(x=X2)) + geom_histogram(binwidth=0.001) + xlab('I_(t+2) / N') + ggtitle('\n') + ylab('') + xlim(0, 0.15)
        plot_ls_It_dex <- plot_ls_It_dex + 1
        plot_ls_It[[plot_ls_It_dex]] <- ggplot(It1, aes(x=X3)) + geom_histogram(binwidth=0.001) + xlab('I_(t+3) / N') + ggtitle('\n') + ylab('') + xlim(0, 0.15)
        plot_ls_It_dex <- plot_ls_It_dex + 1
        
        plot_ls_It1_It[[plot_ls_It1_It_dex]] <- ggplot(It1_It, aes(x=X1)) + geom_histogram(binwidth=0.1) + xlab('I_(t+1) / I_t') + ggtitle(paste0('Day: ', day, '; ',
                                                                                                                                      'R0: ', R0, '\n',
                                                                                                                                      'NComm: ', ncomm, '; ',
                                                                                                                                      'Effect: ', effect)) + ylab('') + xlim(0, 3)
        plot_ls_It1_It_dex <- plot_ls_It1_It_dex + 1
        plot_ls_It1_It[[plot_ls_It1_It_dex]] <- ggplot(It1_It, aes(x=X2)) + geom_histogram(binwidth=0.1) + xlab('I_(t+2) / I_t') + ggtitle('\n') + ylab('') + xlim(0, 3)
        plot_ls_It1_It_dex <- plot_ls_It1_It_dex + 1
        plot_ls_It1_It[[plot_ls_It1_It_dex]] <- ggplot(It1_It, aes(x=X3)) + geom_histogram(binwidth=0.1) + xlab('I_(t+3) / I_t') + ggtitle('\n') + ylab('') + xlim(0, 3)
        plot_ls_It1_It_dex <- plot_ls_It1_It_dex + 1
        
        plot_ls_log_It1_It[[plot_ls_log_It1_It_dex]] <- ggplot(log_It1_It, aes(x=X1)) + geom_histogram(binwidth=0.1) + xlab('log(I_(t+1) / I_t)') + ggtitle(paste0('Day: ', day, '; ',
                                                                                                                                                   'R0: ', R0, '\n',
                                                                                                                                                   'NComm: ', ncomm, '; ',
                                                                                                                                                   'Effect: ', effect)) + ylab('') + xlim(-2, 2)
        plot_ls_log_It1_It_dex <- plot_ls_log_It1_It_dex + 1
        plot_ls_log_It1_It[[plot_ls_log_It1_It_dex]] <- ggplot(log_It1_It, aes(x=X2)) + geom_histogram(binwidth=0.1) + xlab('log(I_(t+2) / I_t)') + ggtitle('\n') + ylab('') + xlim(-2, 2)
        plot_ls_log_It1_It_dex <- plot_ls_log_It1_It_dex + 1
        plot_ls_log_It1_It[[plot_ls_log_It1_It_dex]] <- ggplot(log_It1_It, aes(x=X3)) + geom_histogram(binwidth=0.1) + xlab('log(I_(t+3) / I_t)') + ggtitle('\n') + ylab('') + xlim(-2, 2)
        plot_ls_log_It1_It_dex <- plot_ls_log_It1_It_dex + 1
      }
    }
  }
}
pdf(paste0(output_folder, "It.pdf"), width=10, height=10)
p <- marrangeGrob(grobs=plot_ls_It, nrow=4, ncol=4)
p
dev.off()
pdf(paste0(output_folder, "It1_It.pdf"), width=12, height=10)
p <- marrangeGrob(grobs=plot_ls_It1_It, nrow=3, ncol=4)
p
dev.off()
pdf(paste0(output_folder, "log_It1_It.pdf"), width=12, height=10)
p <- marrangeGrob(grobs=plot_ls_log_It1_It, nrow=3, ncol=4)
p
dev.off()


