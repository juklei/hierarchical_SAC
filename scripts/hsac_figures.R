## Make figures for ltr and lpsac
##
## First edit: 20201120
## Last edit: 20201120
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

require(ggtern)
require(ggplot2)
require(ggtern)
require(cowplot)
require(data.table)
require(scales)

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

dir("clean")

repl_cv <- read.csv("clean/hsac_collector_raw_replication_cv.csv")
site_pred <- read.csv("clean/hsac_collector_raw_site_pred.csv")
pred_perc <- read.csv("clean/hsac_collector_pred_perc.csv")            
max_perc <- read.csv("clean/hsac_collector_max_perc.csv")             
pred_tsp <- read.csv("clean/hsac_collector_pred_tsp.csv") 
specpool_gamma <- read.csv("clean/specpool_gamma.csv")

## 4. Make graphs for hsac model evaluation ------------------------------------

## Compare varition in estimated diversity metrics depending on sample order:

head(repl_cv)

g1 <- ggplot(repl_cv, aes(ntree, X50., color = div_metric))
g2a <- geom_smooth(aes(ntree, X50.), method = "loess", se = FALSE)
g2b <- geom_smooth(aes(ntree, X2.5.), method = "loess", se = FALSE, linetype = "dashed")
g2c <- geom_smooth(aes(ntree, X97.5.), method = "loess", se = FALSE, linetype = "dashed")
g3 <- geom_point(size = 3, alpha= 0.4)
g4 <- geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), 
                    width = 0.5, alpha = 0.4)

png("figures/hsac_repl_cv.png", 8000/4, 6000/4, "px", res = 600/4)
ggplot2:::print.ggplot(g1 + g2a + g2b + g2c + g3 + g4 +
                         ylab("coefficient of variation") + 
                         xlab("number of sampled trees per site") +
                         scale_color_manual(breaks = c("adiv", "bdiv", "gdiv"),
                                            labels = c("alpha", "true beta", "gamma"),
                                            values = c("orange", "blue", "black")) +
                         theme_classic(40) +
                         theme(legend.position = c(0.6, .35), 
                               legend.title = element_blank(),
                               legend.key.size = unit(3, 'lines'),
                               legend.background =  element_blank(),
                               legend.direction = "horizontal"))
dev.off()

## Compare gamma estimates with those of vegan::specpool depending on number of
## tree species:

head(specpool_gamma)
head(site_pred)

sg_mean <- specpool_gamma[, c(1:3, 5, 8, 10)]
sg_mean <- reshape2::melt(sg_mean, 
                          id.vars = c("plot", "Species", "n"), 
                          value.name = "mean_sp")
sg_se <- specpool_gamma[, c(1, 4, 6, 9)]
sg_se <- reshape2::melt(sg_se, id.vars = "plot")
sg_mean$lower_sp <- sg_mean$mean - sg_se$value*1.96 ## SE to 95% CI
sg_mean$upper_sp <- sg_mean$mean + sg_se$value*1.96 ## SE to 95% CI

gcomp <- cbind(sg_mean, site_pred[site_pred$div_metric == "gdiv", 
                                  c("mean", "X2.5.", "X97.5.")])

h1 <- ggplot(gcomp, aes(mean, mean_sp))
h2 <- geom_point(size = 3, alpha = 0.4)#, position = position_dodge(width = 1))
h3 <- geom_errorbarh(aes(xmin = X2.5., xmax = X97.5.), height = 1, alpha = 0.4)
h4 <- geom_errorbar(aes(ymin = lower_sp, ymax = upper_sp), width = 0.5, alpha = 0.4)
h5 <- facet_grid(variable ~ ., scales = "free_y")
h6 <- geom_abline(intercept = 0, slope = 1, linetype = "dashed")

png("figures/gamma_comp.png", 6000/4, 10000/4, "px", res = 600/4)
ggplot2:::print.ggplot(h1 + h6 + h2 + h3 + h4 + h5 +
                       ylab("gamma diversity vegan::specpool") + 
                       xlab("gamma diversity hsac model") +
                       theme_classic(40))
dev.off()

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
ggplot2:::print.ggplot(p1 + p2 + p3 + p4 + p5 +
                       ylab("number/fraction of epiphytic lichen species") + 
                       xlab("number of tree species") +
                       theme_classic(40))
dev.off()

## 6. All tree species group percentages combined ------------------------------

## All percentages combined:

## Chose diversity index here:
# div <- "adiv"
div <- "bdiv"
# div <- "gdiv"

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
# if(div == "bdiv"){Q <- Q + coord_trans(y = "log10")}

## Maxima:
if(div == "bdiv"){
  max_perc$b_max[max_perc$b_max < 0 | max_perc$b_max >1] <- NA
  p1 <- ggplot(max_perc, aes(b_max*100, color = L1, fill = L1, lty = L1)) 
}
if(div == "gdiv"){
  max_perc$g_max[max_perc$g_max < 0 | max_perc$g_max >1] <- NA
  p1 <- ggplot(max_perc, aes(g_max*100, color = L1, fill = L1, lty = L1))
}
if(div %in% c("gdiv", "bdiv")){
  p2 <- geom_density(alpha = .2)
  P <- p1 + p2 + xlim(0, 100) + s1 + s2 + theme0()
  Q <- plot_grid(P, Q, align = "v", nrow = 2, rel_heights = c(1/5, 4/5))
}

png(paste0("figures/hsac_perc_quad_", div, ".png"), 
    8000/4, ifelse(div == "adiv", 6000/4, 8000/4), "px", res = 600/4)
ggplot2:::print.ggplot(Q)
dev.off()
  
## Ternary plot:

head(site_pred)

require(ggtern) ## ggtern breaks ggplot2, load  ggplot2 version 3.2.1 to solve problem
R <- ggtern(droplevels(site_pred[site_pred$div_metric == div, ]), 
            aes(x = dec, y = spruce, z = pine)) +
  geom_mask() +
  geom_point(aes(colour = X50.), size = 10) +
  scale_colour_gradient(low = "yellow", high = "blue") + 
  xlab("") + ylab("") + ggtern::zlab("") +
  labs(colour = ylab, size = 10) +
  theme_classic(40) + 
  theme(legend.position = "bottom",
        legend.key.size = unit(3, 'lines'),
        legend.direction = "horizontal") +
  Tarrowlab("% spruce") + 
  Larrowlab("% deciduous") + 
  Rarrowlab("% pine") +
  theme_showarrows()

png(paste0("figures/hsac_perc_tern_", div, ".png"), 8000/4, 8000/4, "px", res = 600/4)
R
dev.off()

## -------------------------------END-------------------------------------------
