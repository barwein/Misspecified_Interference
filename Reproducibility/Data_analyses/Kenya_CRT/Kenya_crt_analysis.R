###
# Script that analyze the CRT in Kenya conducted by Venturo-Conerly et al., 2022
###


# Load data and libraries -------------------------------------------------


library(data.table)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(ggtext)

source("Reproducibility/Estimation_functions.R")

library(CIwithMIS)

school.df <- fread("Reproducibility/Data_analyses/Kenya_CRT/imputed_data.csv")

# Clean Data --------------------------------------------------------------


# length(unique(school.df$Participant_ID)) # 895 particpants from 2 schools and
# length(unique(school.df$Class)) # 24 classes (12 per school)
#
# table(school.df$Condition) # Study Skills = Control group
# table(school.df$Condition, school.df$School) # Study Skills = Control group

school.df.pre <- school.df[Time==0,]
school.df.post <- school.df[Time==1,]

setkeyv(school.df.pre, "Participant_ID")
setkeyv(school.df.post, "Participant_ID")

# Combine pre- and post- data
school.relevant.only <- school.df.pre[school.df.post]

school.relevant.only <- school.relevant.only[,.(Class, School, Participant_ID, Condition,
                                                GAD_Total, i.GAD_Total,
                                                PHQ_Total, i.PHQ_Total)]

# Compute difference in outcome (post-pre)
school.relevant.only[,`:=`(gad.diff = i.GAD_Total - GAD_Total,
                           phq.diff = i.PHQ_Total - PHQ_Total)]

# Define new version of treatment
school.relevant.only[,Treatment := as.numeric(Condition %in% c("Values","Growth"))]


# ATE under the null ------------------------------------------------------

# Estimate ATE using HT estimator (p(T=1)=0.5 for each unit thus weighting by 1/0.5 = 2)

n <- nrow(school.relevant.only)

gad.class.sum <- tapply(school.relevant.only$gad.diff, school.relevant.only$Class, sum)
gad.trt.classs <- tapply(school.relevant.only$Treatment, school.relevant.only$Class, function(x){all(x==1)})


# ATE and SE based on Aronow&Midelton 2013; they expressed cluster estimators using HT.

null.ATE.GAD <- n^{-1}*2*(sum(gad.class.sum[gad.trt.classs]) - sum(gad.class.sum[!gad.trt.classs]))

null.Var.ATE.GAD <- n^{-2}*(sum(gad.class.sum[gad.trt.classs]^2) - sum(gad.class.sum[!gad.trt.classs]^2))

null.gad.ci.low <- null.ATE.GAD - (1.96 * sqrt(null.Var.ATE.GAD))
null.gad.ci.high <- null.ATE.GAD + (1.96 * sqrt(null.Var.ATE.GAD))

paste0("GAD effect estimate [CI] is ", round(null.ATE.GAD,3),
       " [", round(null.gad.ci.low,3), ", ",
       round(null.gad.ci.high,3),"]")


# PHQ
phq.class.sum <- tapply(school.relevant.only$phq.diff, school.relevant.only$Class, sum)
phq.trt.classs <- tapply(school.relevant.only$Treatment, school.relevant.only$Class, function(x){all(x==1)})

null.ATE.phq <- n^{-1}*2*(sum(phq.class.sum[phq.trt.classs]) - sum(phq.class.sum[!phq.trt.classs]))

null.Var.ATE.phq <- n^{-2}*(sum(phq.class.sum[phq.trt.classs]^2) - sum(phq.class.sum[!phq.trt.classs]^2))

null.phq.ci.low <- null.ATE.phq - (1.96 * sqrt(null.Var.ATE.phq))
null.phq.ci.high <- null.ATE.phq + (1.96 * sqrt(null.Var.ATE.phq))

paste0("PHQ effect estimate [CI] is ", round(null.ATE.phq,3),
       " [", round(null.phq.ci.low,3), ", ",
       round(null.phq.ci.high,3),"]")


# SA ----------------------------------------------------------------------

# Order by Class, School
setorderv(school.relevant.only, c("Class","School"))


school.names <- unique(school.relevant.only$School) # "KDSS" "MGSS"

n.each.KDSS <- c(table(school.relevant.only[School == "KDSS", Class]))
n.each.MGSS <- c(table(school.relevant.only[School == "MGSS", Class]))

# Function that run SA
kenya_school_SA <- function(Z.obs,
                            N_each_cluster_vec,
                            phq.diff,
                            gad.diff,
                            theta_vec,
                            n,
                            R,
                            M,
                            N_clusters){
  SA.results <- lapply(theta_vec, function(theta){
                  rbindlist(lapply(seq(M), function(m){
                  # Generate Adj mat from SBM
                  S.11 <- generate_Q_matrix(N_clusters = N_clusters/2,
                                            within_prob_vec = rep(1,N_clusters/2),
                                            between_prob_vec = theta)
                  S.22 <- generate_Q_matrix(N_clusters = N_clusters/2,
                                            within_prob_vec = rep(1,N_clusters/2),
                                            between_prob_vec = theta)
                  S.12 <- matrix(0, nrow = N_clusters/2, ncol=N_clusters/2)
                  # Generate block-Q matrix
                  Q.mat <- rbind(cbind(S.11,S.12),
                                  cbind(S.12,S.22))
                  # Sample SBM
                  adj.mat <- as_adjacency_matrix(sample_sbm(n = n,
                                               pref.matrix = Q.mat,
                                                 block.sizes = N_each_cluster_vec,
                                                   directed = FALSE))
                  # Compute prob matrix
                  prob.mat <- Get_prob_matrices_list(R = R,
                                                     n = n,
                                                     Pz_function = Z_ber_clusters,
                                                     pz_func_args = list(N_clusters = N_clusters,
                                                                         N_each_cluster_vec = N_each_cluster_vec,
                                                                         p = 0.5),
                                                     A.list = list(A=adj.mat),
                                                     exposures_contrast = list(c("c11","c00")),
                                                     exposures_vec = c("c11","c00"),
                                                     threshold = rep(0,n))
                  # Compute CE (c11-c00)
                  CE.GAD <- MR_CE_estimator(Z.obs = Z.obs,
                                            Y.obs = gad.diff,
                                             A.list = list(A=adj.mat),
                                             exposures_contrast = list(c("c11","c00")),
                                             exposures_vec = c("c11","c00"),
                                             Prob_matrices_list = prob.mat,
                                             threshold = rep(0,n))
                  CE.PHQ <- MR_CE_estimator(Z.obs = Z.obs,
                                            Y.obs = phq.diff,
                                             A.list = list(A=adj.mat),
                                             exposures_contrast = list(c("c11","c00")),
                                             exposures_vec = c("c11","c00"),
                                             Prob_matrices_list = prob.mat,
                                             threshold = rep(0,n))
                  # Save results
                  data.table(iter = m,
                              theta = theta,
                             gad.ht = CE.GAD$`c11-c00`$ht_ce,
                             gad.hajek = CE.GAD$`c11-c00`$hajek_ce,
                             gad.ht.var = CE.GAD$`c11-c00`$var_ht_ce,
                             gad.hajek.var = CE.GAD$`c11-c00`$var_hajek_ce,
                             phq.ht = CE.PHQ$`c11-c00`$ht_ce,
                             phq.hajek = CE.PHQ$`c11-c00`$hajek_ce,
                             phq.ht.var = CE.PHQ$`c11-c00`$var_ht_ce,
                             phq.hajek.var = CE.PHQ$`c11-c00`$var_hajek_ce)
                  }))
  })
  return(rbindlist(SA.results))
}


# Perform SA!
theta_vec <- seq(0.001,0.015,0.001)
set.seed(59542)
SA_run <- kenya_school_SA(Z.obs = school.relevant.only$Treatment,
                       N_each_cluster_vec = c(n.each.KDSS, n.each.MGSS),
                       phq.diff = school.relevant.only$phq.diff,
                       gad.diff = school.relevant.only$gad.diff,
                       theta_vec = theta_vec,
                       n = n,
                       R = 10^3,
                       M = 500,
                       N_clusters = 24)

SA.results <- fread("Reproducibility/Data_analyses/Kenya_CRT/kenya_SA_M500.csv")

nul.results <- data.table(iter = 1, theta = 0,
                          gad.ht = null.ATE.GAD, gad.hajek = null.ATE.GAD ,
                          gad.ht.var = null.Var.ATE.GAD, gad.hajek.var = null.Var.ATE.GAD,
                          phq.ht = null.ATE.phq, phq.hajek = null.ATE.phq,
                          phq.ht.var = null.Var.ATE.phq, phq.hajek.var = null.Var.ATE.phq)

combined.results <- rbindlist(list(nul.results,SA.results))


phq.results.full <- combined.results[,.(theta, phq.hajek, phq.ht)]
gad.results.full <- combined.results[,.(theta, gad.hajek, gad.ht)]
phq.results.full[,outcome := "PHQ"]
gad.results.full[,outcome := "GAD"]
phq.results.full[theta==0,`:=`(ci.high = null.phq.ci.high, ci.low = null.phq.ci.low)]
gad.results.full[theta==0,`:=`(ci.high = null.gad.ci.high, ci.low = null.gad.ci.low)]
names(phq.results.full)[2:3] <- c("hajek","ht")
names(gad.results.full)[2:3] <- c("hajek","ht")
melted.results.full <- rbindlist(list(gad.results.full,phq.results.full))


gad.boxplot <- ggplot(melted.results.full[theta>0 & theta <= 0.01 & outcome=="GAD",],
                aes(x=factor(theta,
                             labels = sprintf("%g%%",seq(0.1,1,0.1))),
                    y=hajek,
                    group=factor(theta,
                                 labels = sprintf("%g%%",seq(0.1,1,0.1))))) +
                geom_boxplot(fill = "gray80", alpha = 0.8) +
                geom_hline(yintercept = 0, lty = "dashed", linewidth = 1) +
                geom_hline(yintercept = c(null.ATE.GAD), col = "#1D8A99", linewidth = 1.2) +
                scale_y_continuous(breaks = round(c(seq(-4,-2,1),null.ATE.GAD,c(0,2,1)),2),
                                   limits = c(-4.9,2.9)) +
                labs(x=TeX("$\\theta$"), y ="", title = "GAD") +
                theme_pubclean() +
                theme(axis.text.x = element_text(size =20, face = "bold", hjust = 0.4),
                      axis.text.y = element_markdown(size =20, face = "bold",
                                                 color = c(rep("grey30",3),"#1D8A99",rep("grey30",3))),
                                                 # color = "#009E73"),
                      axis.title.x = element_text(size = 35),
                      plot.title = element_text(size=28, face = "bold",hjust = 0.5, vjust = -5))

phq.boxplot <- ggplot(melted.results.full[theta>0 & theta <= 0.01 & outcome=="PHQ",],
                      aes(x=factor(theta,
                                   labels = sprintf("%g%%",seq(0.1,1,0.1))),
                          y=hajek,
                          group=factor(theta,
                                       labels = sprintf("%g%%",seq(0.1,1,0.1))))) +
  geom_boxplot(fill = "gray80", alpha = 0.8) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 1) +
  geom_hline(yintercept = c(null.ATE.phq), col = "#1D8A99", linewidth = 1.2) +
  scale_y_continuous(breaks = round(c(seq(-4,-1,1),null.ATE.phq,c(0,2,1)),2),
                     limits = c(-4.9,2.9)) +
  labs(x=TeX("$\\theta$"), y ="", title = "PHQ") +
  theme_pubclean() +
  theme(axis.text.x = element_text(size =20, face = "bold", hjust = 0.4),
        axis.text.y = element_markdown(size =20, face = "bold",
                                       color = c(rep("grey30",4),"#1D8A99",rep("grey30",3))),
        # color = "#009E73"),
        axis.title.x = element_text(size = 35),
        plot.title = element_text(size=28, face = "bold",hjust = 0.5, vjust = -5))

boxplot.kenya <- ggarrange(gad.boxplot,phq.boxplot,
                  nrow = 1, align = "hv")

summary.combined.results <- combined.results[,.(mean.gad.ht = mean(gad.ht),
                                                mean.gad.hajek = mean(gad.hajek),
                                                var.gad.ht = var(gad.ht),
                                                var.gad.hajek = var(gad.hajek),
                                                perc025.gad.ht = quantile(gad.ht, 0.025),
                                                perc975.gad.ht = quantile(gad.ht, 0.975),
                                                perc025.gad.hajek = quantile(gad.hajek, 0.025),
                                                perc975.gad.hajek = quantile(gad.hajek, 0.975),
                                                max.gad.hajek = max(gad.hajek),
                                                min.gad.hajek = min(gad.hajek),
                                                mean.phq.ht = mean(phq.ht),
                                                mean.phq.hajek = mean(phq.hajek),
                                                var.phq.ht = var(phq.ht),
                                                var.phq.hajek = var(phq.hajek),
                                                perc025.phq.ht = quantile(phq.ht, 0.025),
                                                perc975.phq.ht = quantile(phq.ht, 0.975),
                                                perc025.phq.hajek = quantile(phq.hajek, 0.025),
                                                perc975.phq.hajek = quantile(phq.hajek, 0.975),
                                                max.phq.hajek = max(phq.hajek),
                                                min.phq.hajek = min(phq.hajek)),
                                             by = "theta"]
summary.combined.results[1, `:=`(var.gad.ht = null.Var.ATE.GAD,
                                 var.gad.hajek = null.Var.ATE.GAD,
                                 perc025.gad.ht = null.gad.ci.low,
                                 perc975.gad.ht = null.gad.ci.high,
                                 perc025.gad.hajek = null.gad.ci.low,
                                 perc975.gad.hajek = null.gad.ci.high,
                                 max.gad.hajek = NA,
                                 min.gad.hajek = NA,
                                 var.phq.ht = null.Var.ATE.phq,
                                 var.phq.hajek = null.Var.ATE.phq,
                                 perc025.phq.ht = null.phq.ci.low,
                                 perc975.phq.ht = null.phq.ci.high,
                                 perc025.phq.hajek = null.phq.ci.low,
                                 perc975.phq.hajek = null.phq.ci.high,
                                 max.phq.hajek = NA,
                                 min.phq.hajek = NA)
                                 ]

summary.combined.results[,is.null := theta==0]

phq.results <- summary.combined.results[,.(theta, is.null, mean.phq.hajek,
                                           perc025.phq.hajek, perc975.phq.hajek,
                                           max.phq.hajek, min.phq.hajek)]
gad.results <- summary.combined.results[,.(theta, is.null, mean.gad.hajek,
                                           perc025.gad.hajek, perc975.gad.hajek,
                                           max.gad.hajek, min.gad.hajek)]
phq.results[,outcome := "PHQ"]
gad.results[,outcome := "GAD"]
names(phq.results)[3:7] <- c("mean.hajek","perc025","perc975","max","min")
names(gad.results)[3:7] <- c("mean.hajek","perc025","perc975","max","min")

melted.results <- rbindlist(list(gad.results,phq.results))


ci.plot <- ggplot(melted.results[theta <= 0.01,],
                     aes(x=factor(theta, labels = sprintf("%g%%",seq(0,1,0.1))),
                         y = mean.hajek, group = factor(theta),
                    col = is.null)) +
                geom_errorbar(aes(ymin = perc025,
                                  ymax = perc975,
                                  width = is.null/2),
                              # width = 0,
                              linewidth = 1.1, alpha = 0.75,
                              show.legend = F) +
                # scale_linetype_manual(values = c(1,9)) +
                geom_point(cex = 10, show.legend = F, aes(shape=is.null)) +
                # geom_point(aes(y=min),shape = 21, cex = 6, col = "red", alpha = 0.6) +
                # geom_point(aes(y=max),shape = 21, cex = 6, col = "red", alpha = 0.6) +
                scale_shape_manual(values = c(16,18)) +
                geom_hline(yintercept = 0, lty = 2, alpha = 0.5,linewidth = 1) +
                scale_color_manual(values = c("TRUE" = "#009E73","FALSE" = "#000000")) +
                scale_y_continuous(breaks = seq(-2.5,1.5,0.5)) +
                # scale_x_discrete(limits = sprintf("%g%%",seq(0,1,0.1))) +
                labs(x=TeX("$\\theta$"), y ="") +
                facet_wrap(~outcome, nrow = 1) +
                theme_pubclean() +
                theme(axis.text.x = element_text(size =20, face = "bold", hjust = 0.4),
                      axis.text.y = element_text(size =20, face = "bold"),
                      axis.title.x = element_text(size = 35),
                      strip.background = element_blank(),
                      strip.text = element_text(size=26, face = "bold")
                      )


