## Make figures for ltr and lpsac
##
## First edit: 20201120
## Last edit: 20201120
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

require(devtools)
# install_version("ggplot2", version = "3.3.0", repos = "http://cran.us.r-project.org")
require(ggplot2)
# require(ggpubr)
require(cowplot)
require(data.table)

## 2. Define or source functions used in this script ---------------------------

theme0 <- function(...) theme(legend.position = "none",
                              panel.background = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.ticks = element_blank(),
                              axis.text.x = element_blank(),
                              axis.text.y = element_blank(),
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank(),
                              axis.ticks.length = unit(0, "null"),
                              panel.border = element_blank(), ...)

## 3. Load and explore data ----------------------------------------------------

sir("clean")

repl_diff <- read.csv("clean/hsac_collector_raw_replication_diff.csv")
gdiv_comp <- read.csv("clean/hsac_collector_raw_gdiv_comp.csv")
specpool_gamma <- read.csv("clean/specpool_gamma.csv") 
site_pred <- read.csv("clean/hsac_collector_raw_site_pred.csv")
pred_perc <- read.csv("clean/hsac_collector_pred_perc.csv")            
max_perc <- read.csv("clean/hsac_collector_max_perc.csv")             
pred_tsp <- read.csv("clean/hsac_collector_pred_tsp.csv")            

## 4. Make graphs for hsac model evaluation ------------------------------------

##...

## 5. Make graphs for hsac with nr. of tree species ----------------------------

head(pred_tsp)

levels(pred_tsp$div_metric) <- c("alpha diversity", 
                                 "true beta diversity", 
                                 "gamma diversity")

p1 <- ggplot(pred_tsp, aes(x = nr_tsp, y = X50.))
p2 <- geom_line(size = 1, linetype = "dashed", colour = "grey")
p3 <- geom_point(size = 3)
p4 <- geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), width = 0.2)
p5 <- facet_grid(div_metric ~ ., scales = "free_y")

png("figures/hsac_tsp.png", 6000/4, 10000/4, "px", res = 600/4)
p1 + p2 + p3 + p4 + p5 +
  ylab("number/fraction of epiphytic lichen species") + 
  xlab("number of tree species") +
  theme_classic(40) 
dev.off()

## 5. All tree species group percentages combined ------------------------------

## All percentages combined:

## Chose diversity index here:
# div <- "adiv"
# div <- "bdiv"
div <- "gdiv"

head(pred_perc)
head(max_perc)

## Style for both:
s1 <- scale_color_manual(breaks = c("dec", "pine", "spruce"), 
                         values = c("#ffc425", "#d11141", "black"))
s2 <- scale_fill_manual(breaks = c("dec", "pine", "spruce"),
                        values = c("#ffc425", "#d11141", "black"))

## Ylab for ggplot & ggtern:
if(div == "adiv"){ylab <- "alpha diversity"} 
if(div == "bdiv"){ylab <- "true beta diversity"} 
if(div == "gdiv"){ylab <- "gamma diversity"} 

## Predictions:
q1 <- ggplot(droplevels(pred_perc[pred_perc$div_metric == div, ]), 
             aes(x = perc_pred*100, y = X50., 
                 fill = L1, color = L1, lty = L1))
q2 <- geom_line(size = 2)
q3 <- geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = .2, colour = NA)

Q <- q1 + q2 + q3 +
     ylab(ylab) + xlab("percentage of trees") +
     s1 + s2 +
     theme_classic(40) +                  
     theme(legend.position = c(0.4, 0.10), 
           legend.title = element_blank(),
           legend.key.size = unit(3, 'lines'),
           legend.direction = "horizontal")

## Maxima:
if(div == "bdiv"){
  max_perc$b_max[max_perc$b_max < 0 | max_perc$b_max >1] <- NA
  p1 <- ggplot(max_perc, aes(b_max*100, color = L1, fill = L1, lty = L1))
}
if(div == "gdiv"){
  max_perc$g_max[max_perc$g_max < 0 | max_perc$g_max >1] <- NA
  p1_g_max <- ggplot(max_perc, aes(g_max*100, color = L1, fill = L1, lty = L1))
}
if(div %in% c("gdiv", "bdiv")){
  p2 <- geom_density(alpha = .2)
  P <- p1 + p2 + xlim(0, 100) + s1 + s2 + theme0()
  P <- plot_grid(P, Q, align = "v", nrow = 2, rel_heights = c(1/5, 4/5))
}
  
## Ternary plot:

head(site_pred)

y_all <- rbind(as.matrix(plot_estimates[[1]]),
               as.matrix(plot_estimates[[2]]),
               as.matrix(plot_estimates[[3]]))

d_all <- data.frame("r" = apply(y_all, 2, mean),
                    "sd" = apply(y_all, 2, sd),
                    "dec" = rep(plot_estimates$dec, 2),
                    "spruce" = rep(plot_estimates$spruce, 2),
                    "pine" = rep(plot_estimates$pine, 2),
                    "diversity" = c(rep("alpha", length(plot_estimates$dec)),
                                    rep("beta", length(plot_estimates$dec))))

require(ggtern) ## ggtern breaks ggplot2, load  ggplot2 version 3.2.1 to solve problem
R <- ggtern(droplevels(d_all[d_all$diversity == div, ]), 
            aes(x = dec, y = spruce, z = pine)) +
  # stat_density_tern(geom = 'polygon',
  #                   n         = 500,
  #                   aes(colour  = ..level.., alpha = ..level..)) +
  # geom_interpolate_tern(aes(value = richness, colour = ..level..),
  #                       bins = 50,
  #                       alpha = 1) +
  geom_mask() +
  geom_point(aes(colour = r), size = 10) +
  scale_colour_gradient(low = "yellow", high = "blue") + 
  # scale_fill_gradient(low = "white", high = "yellow")  +
  # guides(colour = "", fill = "none", alpha = "none") +
  xlab("") + ylab("") + ggtern::zlab("") +
  labs(colour = ylab, size = 10) +
  theme_classic(40) + 
  theme(legend.position = "bottom",#"c(0.5, -0.05)", 
        legend.key.size = unit(3, 'lines'),
        legend.direction = "horizontal") +
  Tarrowlab("% spruce") + 
  Larrowlab("% deciduous") + 
  Rarrowlab("% pine") +
  theme_showarrows()

R

# # ## Annotate and export combined plots:
# # 
# # PQ <- annotate_figure(PQ, 
# #                       fig.lab = " (a)",
# #                       fig.lab.pos = "top.left", 
# #                       fig.lab.size = 35)
# # R <- annotate_figure(R,
# #                      fig.lab = " (b)", 
# #                      fig.lab.pos = "top.left",
# #                      fig.lab.size = 35)
# 
# png("figures/figure_3_new.png", 16000/4, 8000/4, "px", res = 600/4)
# plot_grid(PQ, R, align = "h", axis = "t", ncol = 2, rel_widths = c(0.52, 0.48))
# dev.off()

png(paste0("figures/figure_3a_new", div, ".png"), 8000/4, 8000/4, "px", res = 600/4)
PQ
dev.off()

png(paste0("figures/figure_3b_new", div, ".png"), 8000/4, 8000/4, "px", res = 600/4)
R
dev.off()

## -------------------------------END-------------------------------------------
