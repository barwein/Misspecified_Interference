###
# Script that generates the graphic results of the NMR bias-variance tradeoff simulations
###


# Load --------------------------------------------------------------------


library(ggplot2)
library(data.table)
library(latex2exp)
library(ggpubr)



# Graphics ----------------------------------------------------------------

MR.sim.results <- fread("Reproducibility/Simulations/results/MR_bias_var_PA_n3000_K6_eta0.25.csv")

MR.sim.results[,true_effect := rep(c(1,0.75,0.5,0.25,0.5),nrow(MR.sim.results)/5)]

mean.sim.results <- MR.sim.results[,.(m_ht = abs(100*mean((ht_ce-true_effect)/true_effect)),
                                   m_hajek = abs(100*mean((hajek_ce-true_effect)/true_effect)),
                                   se_ht = 100*sd(ht_ce/true_effect,na.rm=T),
                                   se_hajek = 100*sd(hajek_ce/true_effect,na.rm=T),
                                   rmse_ht = sqrt((100*mean((ht_ce-true_effect)/true_effect))^2 +
                                                    (100*sd(ht_ce/true_effect,na.rm=T))^2),
                                   rmse_hajek = sqrt((100*mean((hajek_ce-true_effect)/true_effect))^2 +
                                                       (100*sd(hajek_ce/true_effect,na.rm=T))^2)
                                   ),
                                by = c("ce_contrast","K","adj.mat.used","with.true.adj")]



# Main plot c11-c10 -------------------------------------------------------

# [ce_contrast == "c11-c10",]
p.hajek.bias <- ggplot(mean.sim.results[ce_contrast == "c11-c10",],
                         aes(x=factor(K),y = m_hajek,
                             col = with.true.adj, shape = with.true.adj)) +
                    labs(title = "Bias", 
                         x = "",
                         # y = TeX('$\\tau - \\hat{\\tau}$')
                         y = ""
                         ) +
                geom_point(size = 10, stroke = 2) +
                # scale_color_manual(name = "With A*", values = c("red4","green4")) +
                # scale_shape_manual(name = "With A*", values = c(1,4)) +
                scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
                scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
                # scale_y_continuous(breaks = seq(0,0.6,0.1),limits = c(0,0.63)) +
                scale_y_continuous(breaks = seq(0,120,20),limits = c(0,120)) +
                # facet_wrap(~ce_contrast, nrow = 1, scales = "free") +
                theme_pubclean() +
                guides(color = guide_legend(override.aes = list(size = 12))) +
               theme(axis.title.x = element_text(size = 28, face = "bold"),
                  title = element_text(size = 32, face="bold"),
                  axis.text = element_text(size=26, face = "bold"),
                  # legend.title = element_text(size = 14),
                  legend.title = element_blank(),
                  legend.text = element_text(face = "bold", size = 28),
                  legend.key.size = unit(1.2,"cm"),
                  plot.margin=unit(c(1,0.15,1,1),"cm")) 


                
                
# [ce_contrast == "c11-c10",]
p.hajek.sd <- ggplot(mean.sim.results[ce_contrast == "c11-c10",],
                     aes(x=factor(K),y = se_hajek,
                         col = with.true.adj, shape = with.true.adj)) +
                labs(title = "SD", 
                     # x = "",
                     x = "# networks used",
                     # y = TeX('$\\tau - \\hat{\\tau}$')
                     y = ""
                     ) +
              geom_point(size = 10, stroke = 2) +
              # scale_color_manual(name = "With A*", values = c("red4","green4")) +
              # scale_shape_manual(name = "With A*", values = c(1,4)) +
              scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
              scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
              # scale_y_continuous(breaks = seq(0,0.6,0.1),limits = c(0,0.63)) +
              scale_y_continuous(breaks = seq(0,120,20),limits = c(0,120)) +
              # facet_wrap(~ce_contrast, nrow = 1, scales = "free") +
              theme_pubclean() +
              guides(color = guide_legend(override.aes = list(size = 12))) +
              theme(axis.title.x = element_text(size = 28, face = "bold"),
                    title = element_text(size = 32, face="bold"),
                    axis.text = element_text(size=26, face = "bold"),
                    # legend.title = element_text(size = 14),
                    axis.text.y = element_blank(),
                    legend.title = element_blank(),
                    legend.text = element_text(face = "bold", size = 28),
                    legend.key.size = unit(1.2,"cm"),
                    plot.margin=unit(c(1,0.15,1,0.15),"cm")) 



# [ce_contrast == "c11-c10",]
p.hajek.rmse <- ggplot(mean.sim.results[ce_contrast == "c11-c10",],
                     aes(x=factor(K),y = rmse_hajek,
                         col = with.true.adj, shape = with.true.adj)) +
                labs(title = "RMSE", 
                     # x = "# networks used",
                     x = "",
                     # y = TeX('$\\tau - \\hat{\\tau}$')
                     y = ""
                     ) +
                geom_point(size = 10, stroke = 2) +
                # scale_color_manual(name = "With A*", values = c("red4","green4")) +
                # scale_shape_manual(name = "With A*", values = c(1,4)) +
                scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
                scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
                # scale_y_continuous(breaks = seq(0,0.6,0.1),limits = c(0,0.63)) +
                scale_y_continuous(breaks = seq(0,120,20),limits = c(0,120)) +
                # facet_wrap(~ce_contrast, nrow = 1, scales = "free") +
                theme_pubclean() +
                guides(color = guide_legend(override.aes = list(size = 12))) +
                theme(axis.title.x = element_text(size = 28, face = "bold"),
                      title = element_text(size = 32, face="bold"),
                      axis.text = element_text(size=26, face = "bold"),
                      # legend.title = element_text(size = 14),
                      axis.text.y = element_blank(),
                      legend.title = element_blank(),
                      legend.text = element_text(face = "bold", size = 28),
                      legend.key.size = unit(1.2,"cm"),
                      plot.margin=unit(c(1,1,1,0.15),"cm")) 
  



p.hajek.both <- ggarrange(p.hajek.bias,p.hajek.sd,p.hajek.rmse,
                          nrow = 1,ncol = 3,
                          widths = c(1,1,1),
                          common.legend = T, legend = "top", align = "h") 
                          # %>%
                          # annotate_figure(top = text_grob("Hajek estimator; c10-c00",
                          #                 face = "bold", size =14))


# APPENDIX plot -------------------------------------------------------


p.hajek.bias.apdx <- ggplot(mean.sim.results[ce_contrast %in% c("c01-c00","c11-c00"),],
                       aes(x=factor(K),y = m_hajek,
                           col = with.true.adj, shape = with.true.adj)) +
  labs(title = "Bias", 
       x = "",
       # y = TeX('$\\tau - \\hat{\\tau}$')
       y = ""
  ) +
  geom_point(size = 8, stroke = 1.6) +
  # scale_color_manual(name = "With A*", values = c("red4","green4")) +
  # scale_shape_manual(name = "With A*", values = c(1,4)) +
  scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
  scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
  # scale_y_continuous(breaks = seq(0,0.5,0.1),limits = c(0,0.55)) +
  scale_y_continuous(breaks = seq(0,140,20),limits = c(0,140)) +
  facet_wrap(~ce_contrast, nrow = 1, scales = "free") +
  theme_pubclean() +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme(axis.title.x = element_text(size = 26, face = "bold"),
        title = element_text(size = 28, face="bold"),
        axis.text = element_text(size=22, face = "bold"),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 26),
        legend.key.size = unit(1.1,"cm"),
        strip.text = element_text(size=26, face = "bold")) 



p.hajek.sd.apdx <- ggplot(mean.sim.results[ce_contrast %in% c("c01-c00","c11-c00"),],
                     aes(x=factor(K),y = se_hajek,
                         col = with.true.adj, shape = with.true.adj)) +
  labs(title = "SD", 
       x = "",
       # y = TeX('$\\tau - \\hat{\\tau}$')
       y = ""
  ) +
  geom_point(size = 8, stroke = 1.6) +
  # scale_color_manual(name = "With A*", values = c("red4","green4")) +
  # scale_shape_manual(name = "With A*", values = c(1,4)) +
  scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
  scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
  # scale_y_continuous(breaks = seq(0,0.5,0.1),limits = c(0,0.53)) +
  scale_y_continuous(breaks = seq(0,140,20),limits = c(0,140)) +
  facet_wrap(~ce_contrast, nrow = 1, scales = "free") +
  theme_pubclean() +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme(axis.title.x = element_text(size = 26, face = "bold"),
        title = element_text(size = 28, face="bold"),
        axis.text = element_text(size=22, face = "bold"),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 26),
        legend.key.size = unit(1.1,"cm"),
        strip.text = element_text(size=26, face = "bold")) 



p.hajek.rmse.apdx <- ggplot(mean.sim.results[ce_contrast %in% c("c01-c00","c11-c00"),],
                       aes(x=factor(K),y = rmse_hajek,
                           col = with.true.adj, shape = with.true.adj)) +
                    labs(title = "RMSE", 
                         x = "# networks used",
                         # y = TeX('$\\tau - \\hat{\\tau}$')
                         y = ""
                    ) +
                    geom_point(size = 8, stroke = 1.6) +
                    # scale_color_manual(name = "With A*", values = c("red4","green4")) +
                    # scale_shape_manual(name = "With A*", values = c(1,4)) +
                    scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
                    scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
                    # scale_y_continuous(breaks = seq(0,0.5,0.1),limits = c(0,0.53)) +
                    scale_y_continuous(breaks = seq(0,140,20),limits = c(0,140)) +
                    facet_wrap(~ce_contrast, nrow = 1, scales = "free") +
                    theme_pubclean() +
                    guides(color = guide_legend(override.aes = list(size = 10))) +
                    theme(axis.title.x = element_text(size = 26, face = "bold"),
                          title = element_text(size = 28, face="bold"),
                          axis.text = element_text(size=22, face = "bold"),
                          # legend.title = element_text(size = 14),
                          legend.title = element_blank(),
                          legend.text = element_text(face = "bold", size = 26),
                          legend.key.size = unit(1.1,"cm"),
                          strip.text = element_text(size=26, face = "bold")) 
                  
                  



p.hajek.both.apdx <- ggarrange(p.hajek.bias.apdx,p.hajek.sd.apdx,p.hajek.rmse.apdx,
                          nrow = 3,ncol = 1,
                          # widths = c(1,1,2),
                          common.legend = T, legend = "top")



# SE SD comparison --------------------------------------------------------

se_sd_comparison <- MR.sim.results[with.true.adj==TRUE ,
                                   .(hajek.se.sd = sqrt(var_hajek_ce)/sd(hajek_ce),
                                       ht.se.sd = sqrt(var_ht_ce)/sd(ht_ce)),
                                   by = c("adj.mat.used","ce_contrast","K")]

se_sd_comparison <- melt.data.table(se_sd_comparison,
                      id.vars = c("adj.mat.used","ce_contrast","K"),
                      measure.vars = c("hajek.se.sd","ht.se.sd"),
                      variable.name = "esti",
                      value.name = "se.sd")

se_sd_comparison_mean <- MR.sim.results[with.true.adj==TRUE ,
                                   .(hajek.se.sd = mean(sqrt(var_hajek_ce)/sd(hajek_ce),na.rm=T),
                                       ht.se.sd = mean(sqrt(var_ht_ce)/sd(ht_ce),na.rm=T)),
                                   by = c("adj.mat.used","ce_contrast","K")]

se_sd_comparison_mean <- melt.data.table(se_sd_comparison_mean,
                      id.vars = c("adj.mat.used","ce_contrast","K"),
                      measure.vars = c("hajek.se.sd","ht.se.sd"),
                      variable.name = "esti",
                      value.name = "se.sd")

# [ce_contrast=="c01-c00",]
se.to.sd.plot <- ggplot(se_sd_comparison_mean[ce_contrast=="c01-c00",],
               aes(x=factor(K), y = se.sd, col=esti, shape = esti)) +
            geom_point(size=10,alpha= 0.7) +
            geom_hline(yintercept = 1, lty = "dashed", linewidth = 1.1) +
            scale_color_manual(values = c("ht.se.sd" = "#990000","hajek.se.sd" = "#0065A9")
                               ,labels = c("Hajek","HT")
                               ) +
            scale_shape_manual(values = c("ht.se.sd" = 16,"hajek.se.sd" = 18)
                               ,labels = c("Hajek","HT")
                               ) +
            labs(x="# networks used",
                 y = "SE/SD") +
          guides(col = guide_legend(override.aes = list(size = 12))) +
          theme_pubclean() +
          theme(axis.text.x = element_text(size =20, face = "bold"),
                axis.text.y = element_text(size =20, face = "bold"),
                axis.title.x = element_text(size = 26, face="bold", vjust = 0.2),
                axis.title.y = element_text(size = 22, face="bold"),
                legend.position = "top",
                legend.title = element_blank(),
                legend.text = element_text(size = 26, face = "bold"),
                legend.key.size = unit(1.2,"cm"),
          )

