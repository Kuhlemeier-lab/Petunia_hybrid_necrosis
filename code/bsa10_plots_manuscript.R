####    BSRseq analysis    ####

## April 14th, 2022
## Last edited on April 26th, 2022
## 
## Copyright (C) 2022 Marta Binaghi <marta.binaghi at ips.unibe.ch>
## 
## BSRseq plot for manuscript

wdir <- "/xxx/necrosis/"
setwd(wdir)

# libraries -----------------------------------------------------------------
library(ggplot2)


# read data -----------------------------------------------------------
# deltaAF
delta <- read.table("data/raw/bsa_analysis/delta_allele_frequency.csv",
                    sep = ",",
                    header = TRUE,
                    stringsAsFactors = FALSE)
# fix scaffold names
delta$chr_new <- gsub("Scaffold_(\\d+)__\\d+_contigs__length_\\d+", 
                        "Scaffold_\\1",
                        delta$CHROM,
                        perl = TRUE)

# values of thresholds
thresh <- read.table("data/raw/bsa_analysis/thresholds_allele_frequency.txt",
                     sep = "\t",
                     header = TRUE,
                     stringsAsFactors = FALSE)
# genome info
# import chromosomes lengths to calculate cumulative positions from the index
# of the genome
fai <- read.table("data/genomes/Peax403.fasta.fai",
                  header = FALSE,
                  stringsAsFactors = FALSE,
                  colClasses = c("character", "integer", "NULL", "NULL", "NULL"))
# make new positions, where the first scaffold starts at 1, and the next scaffolds
# get added at the end of the previous
fai_scf <- fai[grep("Scaffold", fai$V1), ]
# order scaffolds by decreasing length
fai_scf <- fai_scf[order(fai_scf$V2, decreasing = T), ]
fai_scf$cumulative_start <- c(1,
                              cumsum(fai_scf$V2) + 1)[-length(fai_scf$V2)+1]
# remove Unnecessary info from scf names
fai_scf$chr_new <- gsub("Scaffold_(\\d+)__\\d+_contigs__length_\\d+", 
                        "Scaffold_\\1",
                        fai_scf$V1,
                        perl = TRUE)

# calculate cumulative positions for plotting -----------------------------
# cumulative positions chromosome-wide to concatenate display of un-anchored 
# scaffolds
delta$cum_pos <- rep(0, times = dim(delta)[1])
delta$cum_pos[grep("Scaffold", delta$chr_new)] <- 
  as.integer(delta$POS[grep("Scaffold", delta$chr_new)]) + 
  fai_scf$cumulative_start[match(delta$CHROM[grep("Scaffold", delta$chr_new)], 
                                 fai_scf$V1)]
delta$cum_pos[grep("Chr", delta$chr_new)] <- 
  as.integer(delta$POS[grep("Chr", delta$chr_new)])

# add cumulative chr names where the scaffolds have only one name
delta$chr_name_cum <- gsub("Scaffold.*", "Scaffolds", delta$chr_new)

rm(fai_scf)

# plot deltaAF ------------------------------------------------------------
# supplementary figure
lowQuantile1 <- 0.01
highQuantile1 <- 0.99
lowQuantile2 <- 0.05
highQuantile2 <- 0.95
cutoff1L <- thresh$value[thresh$quantileProb == lowQuantile1]
cutoff1H <- thresh$value[thresh$quantileProb == highQuantile1]
cutoff2L <- thresh$value[thresh$quantileProb == lowQuantile2]
cutoff2H <- thresh$value[thresh$quantileProb == highQuantile2]

# horizontal plot
ggplot(data = delta,
       aes(cum_pos/1000000,
           (deltaAF))) +
  geom_point(col = "gray25", alpha = 0.2, size = 0.2 ) +
  # geom_smooth(method = "loess",
  #             size = 0.4,
  #             span = 0.1,
  #             se = TRUE,
  #             level = 0.95) +
  geom_hline(yintercept = c(cutoff1H, cutoff1L),
             size = 0.3, lty = 5, col = "black") +
  geom_hline(yintercept = c(cutoff2H, cutoff2L),
              size = 0.3, lty = 6) +
  facet_wrap(~chr_name_cum, ncol = 8, scales = "free_x") +
  #scale_x_continuous(limits = c(-100, 190000000), expand=c(0,0)) +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 0.2),
        text = element_text(size = 16),
        axis.text = element_text(size = 10)) +
  #ggtitle(expression(Delta * " SNP frequency")) +
  xlab("Position (Mb)") +
  ylab(expression(Delta * " SNP frequency"))
ggsave(filename = "figures/bsa/MS/plotMS_delta_allele_frequency_horiz_2_10_cutoff.png", 
       width = 10,
       height = 3,
       dpi = 600)
ggsave(filename = "figures/bsa/MS/plotMS_delta_allele_frequency_horiz_2_10_cutoff.pdf", 
       width = 10,
       height = 3,
       dpi = 600)
# plot proportion of deltaAF outside cutoffs ------------------------------

delta_df_wnd_cutoffs <- read.csv("data/raw/bsa_analysis/proportion_deltaAF_outsideCutoffs.csv",
                                 stringsAsFactors = FALSE,
                                 header = TRUE)

# remove scaffolds because they have too few variants and as seen in the delta AF 
# plot, they do not pass the cutoffs
delta_df_wnd_chr <- delta_df_wnd_cutoffs[grep("Chr", delta_df_wnd_cutoffs$chromosome), ]
# plot proportion of above threshold delta AF, only 2%
ggplot(data = delta_df_wnd_chr[delta_df_wnd_chr$cutoff == "2 %", ],
       aes(position/100000000,
           frequency)) +
  geom_line() +
  #geom_point(alpha = 0.2) +
  facet_wrap(~chromosome, nrow = 1) +
  # scale_color_brewer(palette = "Dark2",
  #                    name = "Cutoff") +
  # scale_alpha(guide = F) +
  theme_classic() +
  scale_x_continuous(breaks = c(0, 1, 2)) +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 0.2)) +
  #ggtitle(expression("Proportion of " * Delta * " SNP frequency outside cutoff")) +
  xlab("Position (100 Mb)") +
  ylab(expression("Proportion of " * Delta * " SNP frequency")) 
ggsave(filename = "figures/bsa/MS/plotMS_proportion_delta_allele_frequency_cutoff2_horizontal_10x3.pdf",
       width = 10, height = 3)
ggsave(filename = "figures/bsa/MS/plotMS_proportion_delta_allele_frequency_cutoff2_horizontal_8x3.pdf",
       width = 8, height = 3)
ggsave(filename = "figures/bsa/MS/plotMS_proportion_delta_allele_frequency_cutoff2_horizontal_6x2.pdf",
       width = 6, height = 2)
# plot proportion of above threshold delta AF, only 10%
ggplot(data = delta_df_wnd_chr[delta_df_wnd_chr$cutoff == "10 %", ],
       aes(position/100000000,
           frequency)) +
  geom_line() +
  #geom_point(alpha = 0.2) +
  facet_wrap(~chromosome, nrow = 1) +
  # scale_color_brewer(palette = "Dark2",
  #                    name = "Cutoff") +
  # scale_alpha(guide = F) +
  theme_classic() +
  scale_x_continuous(breaks = c(0, 1, 2)) +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 0.2)) +
  #ggtitle(expression("Proportion of " * Delta * " SNP frequency outside cutoff")) +
  xlab("Position (100 Mb)") +
  ylab(expression("Proportion of " * Delta * " SNP frequency")) 
ggsave(filename = "figures/bsa/MS/plotMS_proportion_delta_allele_frequency_cutoff10_horizontal_10x3.pdf",
       width = 10, height = 3)
ggsave(filename = "figures/bsa/MS/plotMS_proportion_delta_allele_frequency_cutoff10_horizontal_8x3.pdf",
       width = 8, height = 3)
ggsave(filename = "figures/bsa/MS/plotMS_proportion_delta_allele_frequency_cutoff10_horizontal_6x2.pdf",
       width = 6, height = 2)


# session info ------------------------------------------------------------
sessionInfo()
# R version 3.3.3 (2017-03-06)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 9 (stretch)
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
# [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_3.2.0
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.2       withr_2.1.2      assertthat_0.2.1 crayon_1.3.4     dplyr_0.8.3      grid_3.3.3       R6_2.4.0        
# [8] gtable_0.3.0     magrittr_1.5     scales_1.0.0     pillar_1.4.2     rlang_0.4.0      lazyeval_0.2.2   rstudioapi_0.10 
# [15] labeling_0.3     tools_3.3.3      glue_1.3.1       purrr_0.3.2      munsell_0.5.0    yaml_2.2.0       pkgconfig_2.0.2 
# [22] colorspace_1.4-1 tidyselect_0.2.5 tibble_2.1.3    
