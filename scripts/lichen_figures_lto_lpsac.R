## Make figures for ltr and lpsac
##
## First edit: 20190612
## Last edit: 20200710
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

require(devtools)
# install_version("ggplot2", version = "3.3.0", repos = "http://cran.us.r-project.org")
require(ggplot2)
require(ggpubr)
require(cowplot)
require(rjags)
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

lto <- read.csv("results/lto_tsp.csv")
lto <- lto[lto$identity == "Richness", - 1]
load("clean/export_dec.rda")
load("clean/export_spruce.rda")
load("clean/export_pine.rda")
load("clean/export_tsp.rda")
load("clean/plot_estimates.rda")

## 4. Make graphs for lto predictions ------------------------------------------

## Group mean plot for tree species:

## Make data set:
lto <- t(lto)

lto_r <- sapply(strsplit(lto, split = "(", fixed = TRUE), FUN = "[[", 1)
lto_lu <- sapply(strsplit(lto, split = "(", fixed = TRUE), FUN = "[[", 2)
lto_l <- sapply(strsplit(lto_lu, split = "-"), FUN = "[[", 1)
lto_u <- sapply(strsplit(lto_lu, split = "-"), FUN = "[[", 2)
lto_u <- unlist(strsplit(lto_u, ")"))

d_ls <- data.frame("richness" = as.numeric(lto_r), 
                   "lower" = as.numeric(lto_l),
                   "upper" = as.numeric(lto_u),
                   "species" = rownames(lto))
d_ls$species <- factor(d_ls$species, 
                       levels =  c("tree", "aspen", "birch", "oak", "alder", 
                                   "pine", "spruce"))
levels(d_ls$species) <- c("tree", "Aspen", "Birch spp.", "Oak", "Alder", "Pine", "Spruce")

g1 <- ggplot(d_ls[d_ls$species != "tree", ], aes(x = species, y = richness))
g2 <- geom_point(size = 4)
g3 <- geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3)
g4 <- geom_hline(yintercept = unlist(d_ls[d_ls$species == "tree", 1:3]),
                 linetype = c("solid", "dashed", "dashed"),
                 alpha = 0.4)
g5 <- annotate("text", 
               x = c("Birch spp.", "Oak", "Alder", "Pine", "Spruce"), 
               y = 4, 
               label = c("A", "A,B,C", "B", "C,D", "D"), 
               size = 12)

png("figures/figure_1.png", 10000/4, 7000/4, "px", res = 600/4)
g1 + g2 + g3 + g4 + g5 + 
  theme_classic(45) + ylab("species richness per tree") + xlab("")
dev.off()

## 5. Make graphs for lpsac ----------------------------------------------------

## Nr. tree species:

y_nr_tsp <- rbind(summary(export_tsp[, "r_1tsp"])$quantiles,
                  summary(export_tsp[, "r_2tsp"])$quantiles,
                  summary(export_tsp[, "r_3tsp"])$quantiles,
                  summary(export_tsp[, "r_4tsp"])$quantiles,
                  summary(export_tsp[, "s_1tsp"])$quantiles,
                  summary(export_tsp[, "s_2tsp"])$quantiles,
                  summary(export_tsp[, "s_3tsp"])$quantiles,
                  summary(export_tsp[, "s_4tsp"])$quantiles)

d_lntsp <- data.frame("response" = y_nr_tsp[,3],
                      "lower" = y_nr_tsp[,1],
                      "upper" = y_nr_tsp[,5],
                      "nr_tsp" = rep(1:4, 2),
                      "cat" = c(rep("SAC asymptote (alpha diversity)", 4), 
                                rep("SAC half-saturation (beta diversity)", 4)))

## Annotation data:
ann_text <- data.frame(nr_tsp = c(2, 3), 
                       response = 20, 
                       cat = factor(d_lntsp$cat[1]))

p1 <- ggplot(d_lntsp, aes(x = nr_tsp, y = response))
p2 <- geom_line(size = 1, linetype = "dashed")
p3 <- geom_point(size = 3)
p4 <- geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)
p5 <- facet_grid(cat ~ ., scales = "free_y")
p6 <- geom_text(data = ann_text, 
                label = "A", 
                show.legend = FALSE, 
                size = 11, 
                colour = "black")
  
png("figures/figure_2_new.png", 8000/4, 10000/4, "px", res = 600/4)
p1 + p2 + p3 + p4 + p5 + p6 +
  ylab("        sampled trees                    expected species richness") + 
  xlab("number of tree species") +
  theme_classic(40) 
dev.off()

## All percentages combined:

## Chose diversity index here:
div <- "alpha"
# div <- "beta"

y_all <- rbind(summary(export_dec$r_dec)$quantiles,
               summary(export_spruce$r_spruce)$quantiles,
               summary(export_pine$r_pine)$quantiles,
               summary(export_dec$s_dec)$quantiles,
               summary(export_spruce$s_spruce)$quantiles,
               summary(export_pine$s_pine)$quantiles)

d_all <- data.frame("r" = y_all[,3], "lower" = y_all[,1], "upper" = y_all[,5])
d_all$perc_tree <-  c(export_dec$dec_pred, 
                      export_spruce$spruce_pred,
                      export_pine$pine_pred)
d_all$tree_species <-  c(rep("deciduous", length(export_dec$dec_pred)),
                         rep("spruce", length(export_spruce$spruce_pred)),
                         rep("pine", length(export_pine$pine_pred)))
d_all$diversity <- c(rep("alpha", nrow(d_all)/2), rep("beta", nrow(d_all)/2))

## Maxima:
maxima <- data.frame("value" = c(export_dec$r_dec_max,
                                 export_spruce$r_spruce_max,
                                 export_pine$r_pine_max,
                                 export_dec$s_dec_max,
                                 export_spruce$s_spruce_max,
                                 export_pine$s_pine_max))
maxima$tree_species <-  c(rep("deciduous", length(export_dec$r_dec_max)),
                          rep("spruce", length(export_spruce$r_spruce_max)),
                          rep("pine", length(export_pine$r_pine_max)))
maxima$diversity <- c(rep("alpha", nrow(maxima)/2), rep("beta", nrow(maxima)/2))

## Style for both:
s1 <- scale_color_manual(breaks = c("deciduous", "pine", "spruce"), 
                         values = c("#ffc425", "#d11141", "black"))
s2 <- scale_fill_manual(breaks = c("deciduous", "pine", "spruce"),
                        values = c("#ffc425", "#d11141", "black"))

## Ylab for ggplot & ggtern:
if(div == "alpha"){ylab <- "SAC asymptote (alpha diversity)"} else{
  ylab <- "SAC half-saturation (beta diversity)"
}

## Predictions:
q1 <- ggplot(droplevels(d_all[d_all$diversity == div, ]), 
             aes(x = perc_tree*100, 
                 y = r, 
                 fill = tree_species, 
                 color = tree_species,
                 lty = tree_species))
q2 <- geom_line(size = 2)
q3 <- geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2, colour = NA)

Q <- q1 + q2 + q3 +
     ylab(ylab) + xlab("percentage of trees") +
     s1 + s2 +
     theme_classic(40) +                  
     theme(legend.position = c(0.4, 0.10), 
           legend.title = element_blank(),
           legend.key.size = unit(3, 'lines'),
           legend.direction = "horizontal")

## Maxima:
p1 <- ggplot(data = droplevels(maxima[maxima$diversity == div, ]), 
             aes(x = value*100, 
                 color = tree_species, 
                 fill = tree_species,
                 lty = tree_species))
p2 <- geom_density(alpha = .2)
P <- p1 + p2 + xlim(0, 100) + s1 + s2 + theme0()

PQ <- plot_grid(P, Q, align = "v", nrow = 2, rel_heights = c(1/5, 4/5))

## Ternary plot:

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
