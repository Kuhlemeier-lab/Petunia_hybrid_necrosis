####    BSRseq analysis    ####

## July 22nd, 2019
## Last edited on April 14th, 2022
## 
## Copyright (C) 2019 Marta Binaghi <marta.binaghi@ips.unibe.ch>
## 
## BSRseq analysis in necrosis experiment.

wdir <- "necrosis/"
setwd(wdir)

# libraries -----------------------------------------------------------------
library(vcfUtils)
library(tidyr)
library(dplyr)
library(ggplot2)

# my functions ------------------------------------------------------------
# splits the dataset in stepping windows and calculates the frequency in each
# window of values above a cutoff or outside two cutoffs.
stepping_window <- function(data, cutoff, window_size = 20)
{
  # define frequency calculation depending on one or two values for the cutoff
  if (length(cutoff) == 1) {
    myfreq <- function(x) {
      return(
        sum(x > cutoff) / window_size
      )
    }
  } else if (length(cutoff) == 2 ) {
    myfreq <- function(x) {
      return(
        sum(x > max(cutoff) | x < min(cutoff)) / window_size
      )
    }
  }
  AFfreq <- c()
  # sort by position
  data <- data[order(data$POS), ]
  # matrix() recicles the values to fill the last column, if there are not enough
  # elements. Therefore I have to prevent that.
  rmndr <- nrow(data) %% window_size
  trimmed_df <- data[1 : (nrow(data) - rmndr), ]
  # calculate frequency on trimmed df
  AFfreq <- c(AFfreq, apply(matrix(trimmed_df$deltaAF, nrow = window_size),
                            2,
                            myfreq) )
  # calculate position of the window (as middle between first and last position)
  start_pos <- trimmed_df$POS[c(T, rep(F, window_size - 1))]
  #print(start_pos)
  end_pos <- trimmed_df$POS[c(rep(F, window_size - 1), T)]
  if (rmndr != 0) {
    rmndr_df <- data[(nrow(data) - rmndr + 1) : nrow(data), ]
    # calculate freq of remainder
    AFfreq <- c(AFfreq, myfreq(rmndr_df$deltaAF))
    start_pos <- c(start_pos,
                   rmndr_df$POS[1])
    end_pos <- c(end_pos,
                 rmndr_df$POS[nrow(rmndr_df)])
  }
  #print(end_pos)
  wnd_pos <- ((end_pos - start_pos) / 2 ) + start_pos
  return (
    data.frame(
      position = wnd_pos,
      frequency = AFfreq
    )
  )
}


# import genome info ----------------------------------------------------
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
fai_scf <- fai_scf %>%
  mutate(chr_new = gsub("Scaffold_(\\d+)__\\d+_contigs__length_\\d+", "Scaffold_\\1", V1, perl = T))

# Read vcf  -------------------------------------------
# note that sample 23 and 24 IDs correspond to the SRA ID last digits,
# so 23 is healthy individuals, P. axillaris N X P. exserta
# and 24 is necrotic individuals, P. axillaris N X P. exserta

vcf <- read_vcf("data/raw/bsa_variants/BSA_SNP_biallelic_gatkselected_minDP100_thin100.recode.vcf")
colnames(vcf)[1] <- "CHROM"
colnames(vcf)[10] <- "healthy"
colnames(vcf)[11] <- "necrotic"


# Split genotypes of samples
vcf_gtsplit <- vcf %>%
  extract(healthy,
          into = paste0("h_", unlist(strsplit(vcf$FORMAT[1], split = ":"))),
          '(\\d/\\d):(\\d+,\\d+):(\\d+):(\\d+):(\\d+,\\d+,\\d+)',
          remove = FALSE) %>%
  extract(necrotic,
          into = paste0("n_", unlist(strsplit(vcf$FORMAT[1], split = ":"))),
          '(\\d/\\d):(\\d+,\\d+):(\\d+):(\\d+):(\\d+,\\d+,\\d+)',
          remove = FALSE) %>%
  mutate(CHROM_new = gsub("Scaffold_(\\d+)__.*", "Scaffold_\\1", CHROM, perl = T)
  )

vcf_gtsplit$cum_pos <- rep(0, times = dim(vcf_gtsplit)[1])
vcf_gtsplit$cum_pos[grep("Scaffold", vcf_gtsplit$CHROM_new)] <- 
  as.integer(vcf_gtsplit$POS[grep("Scaffold", vcf_gtsplit$CHROM_new)]) + 
  fai_scf$cumulative_start[match(vcf_gtsplit$CHROM[grep("Scaffold", vcf_gtsplit$CHROM_new)], 
                                 fai_scf$V1)]


# calculate proportion of genotypes per chromosome ------------------------
genotypes <- vcf_gtsplit[ , c("CHROM", "POS", "h_GT", "n_GT")]
genotypes$h_GT <- as.factor(genotypes$h_GT)
genotypes$n_GT <- as.factor(genotypes$n_GT)

proportions <- data.frame(hom_ref = numeric(),
                          het = numeric(),
                          hom_alt = numeric(),
                          stringsAsFactors = FALSE)
prop.tmp.names <- data.frame(chromosome = character(),
                          sample = character(),
                          stringsAsFactors = FALSE)

for ( chr in levels(as.factor(genotypes$CHROM))) {
  gt_chr <- genotypes[genotypes$CHROM == chr, ]
  proportions <- rbind(proportions, summary(gt_chr$h_GT), stringsAsFactors = FALSE)
  prop.tmp.names <- rbind(prop.tmp.names, c(chr, "healthy"), stringsAsFactors = FALSE)
  proportions <- rbind(proportions, summary(gt_chr$n_GT), stringsAsFactors = FALSE)
  prop.tmp.names <- rbind(prop.tmp.names, c(chr, "necrotic"), stringsAsFactors = FALSE)
}

proportions <- cbind(prop.tmp.names, proportions)
colnames(proportions) <- c("chromosome", "sample", "hom_ref", "het", "hom_alt")
rm(prop.tmp.names)
proportions$snp_num <- apply(X = proportions[ , 3:5], MARGIN = 1, FUN = sum)
proportions$prop_hom_ref <- proportions$hom_ref / proportions$snp_num
proportions$prop_het <- proportions$het / proportions$snp_num
proportions$prop_hom_alt <- proportions$hom_alt / proportions$snp_num

write.csv(proportions, file = "data/raw/bsa_analysis/genotype_proportions_per_sample.csv",
          quote = FALSE, row.names = FALSE)

# calculate allele frequency ----------------------------------------------
af_wide <- vcf_gtsplit %>% 
  extract(h_AD,
          into = c("h_AD_REF", "h_AD_ALT"),
          '(\\d+),(\\d+)',
          remove = FALSE) %>%
  extract(n_AD,
          into = c("n_AD_REF", "n_AD_ALT"),
          '(\\d+),(\\d+)',
          remove = FALSE) %>%
  select(CHROM, CHROM_new, POS, cum_pos,
         h_AD_REF, h_AD_ALT, h_DP,
         n_AD_REF, n_AD_ALT, n_DP) %>%
  mutate(POS = as.numeric(POS),
         h_AD_REF = as.numeric(h_AD_REF),
         h_AD_ALT = as.numeric(h_AD_ALT),
         h_DP = as.numeric(h_DP),
         n_AD_REF = as.numeric(n_AD_REF),
         n_AD_ALT = as.numeric(n_AD_ALT),
         n_DP = as.numeric(n_DP)) %>%
  mutate(h_AF = h_AD_ALT / h_DP,
         n_AF = n_AD_ALT / n_DP)

af <- gather(af_wide,
             "sample",
             "AF",
             n_AF:h_AF)

af$sample[grep("h_AF", af$sample)] <- "Healthy"
af$sample[grep("n_AF", af$sample)] <- "Necrotic"

af <- af[ , c(1,2,3,4,11,12)]
options(scipen = 999)
af$chr_name <- af$CHROM_new
af$chr_name[grepl("Scaffold", af$chr_name)] <- "Scaffold"
af$pos_plot <- af$POS
af$pos_plot[grepl("Scaffold", af$chr_name)] <- af$cum_pos[grepl("Scaffold", af$chr_name)]

write.csv(af[ , c(1, 3, 5, 6)], 
          file = "data/raw/bsa_analysis/allele_frequency.csv", 
          row.names = FALSE, 
          quote = FALSE)

# plot allele frequency in each sample ------------------------------------
ggplot(data = af,
       aes(pos_plot,
           AF)) +
  geom_point(aes(color = sample), alpha = 0.1) +
  facet_wrap(~chr_name, ncol = 1) +
  scale_color_manual(values = c("#5ab4ac", "#d8b365"),
                     name = "Sample") +
  guides(colour = guide_legend(override.aes = list(alpha=1))) +
  scale_x_continuous(limits = c(-100, 190000000), expand=c(0,0)) +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 0.2)) +
  ggtitle("Allele frequency") +
  xlab("Position") +
  ylab("Allele frequency (alt reads / read depth)")
ggsave(filename = "figures/bsa/plot_allele_frequency_samples.pdf", width = 10, height = 15)


# delta allele frequency --------------------------------------------------
# calculate delta AF as necrotic AF - healthy AF
afH <- af[af$sample == "Healthy", ]
afN <- af[af$sample == "Necrotic", ]
delta_df <- inner_join(afH, afN, by = c("CHROM", "CHROM_new", "POS", "cum_pos", "chr_name", "pos_plot"))
delta_df$deltaAF <- delta_df$AF.y - delta_df$AF.x

write.csv(delta_df[ , c(1, 3, 11)], 
          file = "data/raw/bsa_analysis/delta_allele_frequency.csv", 
          row.names = FALSE,
          quote = FALSE)

# calculate the cutoff thresholds of the deltaAF over all the genome
# on real values (not absolute), with high and low tails of 1, 2.5 and 5%
cutoff10H <- quantile(delta_df$deltaAF, probs = 0.95) 
cutoff10L <- quantile(delta_df$deltaAF, probs = 0.05) 
# on real values, with high and low tails of 2.5%
cutoff5H <- quantile(delta_df$deltaAF, probs = 0.975) 
cutoff5L <- quantile(delta_df$deltaAF, probs = 0.025)
# on real values, with high and low tails of 1%
cutoff2H <- quantile(delta_df$deltaAF, probs = 0.99)
cutoff2L <- quantile(delta_df$deltaAF, probs = 0.01) 

# save them in a text file
fileConn <- file("data/raw/bsa_analysis/thresholds_allele_frequency.txt")
writeLines(c("quantileProb\tvalue",
             paste0("0.01\t", cutoff2L),
             paste0("0.99\t", cutoff2H),
             paste0("0.025\t", cutoff5L),
             paste0("0.075\t", cutoff5H),
             paste0("0.05\t", cutoff10L),
             paste0("0.95\t", cutoff10H)), fileConn)
close(fileConn)

# plot delta allele frequency
ggplot(data = delta_df,
       aes(pos_plot,
           (deltaAF))) +
  geom_point(col = "grey", alpha = 0.2, size = 0.2 ) +
  geom_smooth(method = "loess",
              size = 0.4,
              span = 0.1,
              se = TRUE,
              level = 0.95) +
  geom_hline(yintercept = c(cutoff2H, cutoff2L),
             size = 0.2, lty = 2) +
  geom_hline(yintercept = c(cutoff5H, cutoff5L),
             size = 0.2, lty = 3) +
  geom_hline(yintercept = c(cutoff10H, cutoff10L),
             size = 0.2, lty = 4) +
  facet_wrap(~chr_name, ncol = 1) +
  scale_x_continuous(limits = c(-100, 190000000), expand=c(0,0)) +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 0.2)) +
  ggtitle(expression(Delta * " SNP frequency")) +
  xlab("Position") +
  ylab(expression(Delta * " SNP frequency"))
ggsave(filename = "figures/bsa/plot_delta_allele_frequency.pdf", width = 10, height = 15)


# calculate proportions ---------------------------------------------------
# calculate in windows of 100 SNPs, how many pass the cutoffs

# 10% cutoff
delta_df_wnd <- data.frame(chromosome = character(),
                           position = numeric(),
                           frequency = numeric(),
                           stringsAsFactors = FALSE)
for (chr in levels(as.factor(delta_df$CHROM_new)) ) {
  print(paste0("Chromosome: ", chr))
  data_chr <- delta_df[delta_df$CHROM_new == chr, ]
  if (nrow(data_chr) < 100) {
    wnd <- stepping_window(data_chr, c(cutoff10L, cutoff10H), nrow(data_chr))
    wnd_df <- cbind(chromosome = rep(chr, nrow(wnd)),
                    wnd,
                    stringsAsFactors = F)
    delta_df_wnd <- rbind(delta_df_wnd,
                          wnd_df,
                          stringsAsFactors = F)
  } else {
    wnd <- stepping_window(data_chr, c(cutoff10L, cutoff10H), 100)
    wnd_df <- cbind(chromosome = rep(chr, nrow(wnd)),
                    wnd,
                    stringsAsFactors = F)
    delta_df_wnd <- rbind(delta_df_wnd,
                          wnd_df,
                          stringsAsFactors = F)
  }
}
delta_df_wnd_cutoffs <- cbind(cutoff = c("10 %"),
                              delta_df_wnd)

# 5% cutoff
delta_df_wnd <- data.frame(chromosome = character(),
                           position = numeric(),
                           frequency = numeric(),
                           stringsAsFactors = F)
for (chr in levels(as.factor(delta_df$CHROM_new)) ) {
  print(paste0("Chromosome: ", chr))
  data_chr <- delta_df[delta_df$CHROM_new == chr, ]
  if (nrow(data_chr) < 100) {
    wnd <- stepping_window(data_chr, c(cutoff5L, cutoff5H), nrow(data_chr))
    wnd_df <- cbind(chromosome = rep(chr, nrow(wnd)),
                    wnd,
                    stringsAsFactors = F)
    delta_df_wnd <- rbind(delta_df_wnd,
                          wnd_df,
                          stringsAsFactors = F)
  } else {
    wnd <- stepping_window(data_chr, c(cutoff5L, cutoff5H), 100)
    wnd_df <- cbind(chromosome = rep(chr, nrow(wnd)),
                    wnd,
                    stringsAsFactors = F)
    delta_df_wnd <- rbind(delta_df_wnd,
                          wnd_df,
                          stringsAsFactors = F)
  }
}
delta_df_wnd <- cbind(cutoff = c("5 %"),
                      delta_df_wnd)
delta_df_wnd_cutoffs <- rbind(delta_df_wnd_cutoffs,
                              delta_df_wnd)
#2% cutoff
delta_df_wnd <- data.frame(chromosome = character(),
                           position = numeric(),
                           frequency = numeric(),
                           stringsAsFactors = F)
for (chr in levels(as.factor(delta_df$CHROM_new)) ) {
  print(paste0("Chromosome: ", chr))
  data_chr <- delta_df[delta_df$CHROM_new == chr, ]
  if (nrow(data_chr) < 100) {
    wnd <- stepping_window(data_chr, c(cutoff2L, cutoff2H), nrow(data_chr))
    wnd_df <- cbind(chromosome = rep(chr, nrow(wnd)),
                    wnd,
                    stringsAsFactors = F)
    delta_df_wnd <- rbind(delta_df_wnd,
                          wnd_df,
                          stringsAsFactors = F)
  } else {
    wnd <- stepping_window(data_chr, c(cutoff2L, cutoff2H), 100)
    wnd_df <- cbind(chromosome = rep(chr, nrow(wnd)),
                    wnd,
                    stringsAsFactors = F)
    delta_df_wnd <- rbind(delta_df_wnd,
                          wnd_df,
                          stringsAsFactors = F)
  }
}
delta_df_wnd <- cbind(cutoff = c("2 %"),
                      delta_df_wnd)
delta_df_wnd_cutoffs <- rbind(delta_df_wnd_cutoffs,
                              delta_df_wnd)

write.csv(delta_df_wnd_cutoffs, file = "data/raw/bsa_analysis/proportion_deltaAF_outsideCutoffs.csv",
          quote = FALSE,
          row.names = FALSE)

# remove scaffolds because they have too few variants
delta_df_wnd_chr <- delta_df_wnd_cutoffs[grep("Chr", delta_df_wnd_cutoffs$chromosome), ]
# plot proportion of above threshold delta AF
ggplot(data = delta_df_wnd_chr,
       aes(position/1000000,
           frequency,
           colour = cutoff, alpha = .5)) +
  geom_line() +
  geom_point(alpha = 0.2) +
  facet_wrap(~chromosome, nrow = 7) +
  scale_color_brewer(palette = "Dark2",
                     name = "Cutoff") +
  scale_alpha(guide = F) +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 0.2)) +
  ggtitle(expression("Proportion of " * Delta * " SNP frequency outside cutoff")) +
  xlab("Position (Mb)") +
  ylab(expression("Proportion of extreme " * Delta * " SNP frequency")) 
ggsave(filename = "figures/bsa/plot_proportion_delta_allele_frequency_allcutoffs.pdf", width = 10, height = 15)

# plot proportion of above threshold delta AF each cutoff
ggplot(data = delta_df_wnd_chr[delta_df_wnd_chr$cutoff == "10 %", ],
       aes(position/1000000,
           frequency,
           colour = cutoff, alpha = .5)) +
  geom_line() +
  geom_point(alpha = 0.2) +
  facet_wrap(~chromosome, nrow = 7) +
  scale_color_brewer(palette = "Dark2",
                     name = "Cutoff") +
  scale_alpha(guide = F) +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 0.2)) +
  ggtitle(expression("Proportion of " * Delta * " SNP frequency outside cutoff")) +
  xlab("Position (Mb)") +
  ylab(expression("Proportion of extreme " * Delta * " SNP frequency"))
ggsave(filename = "figures/bsa/plot_proportion_delta_allele_frequency_cutoffs10.pdf", width = 10, height = 15)

ggplot(data = delta_df_wnd_chr[delta_df_wnd_chr$cutoff == "5 %", ],
       aes(position/1000000,
           frequency,
           colour = cutoff, alpha = .5)) +
  geom_line() +
  geom_point(alpha = 0.2) +
  facet_wrap(~chromosome, nrow = 7) +
  scale_color_brewer(palette = "Dark2",
                     name = "Cutoff") +
  scale_alpha(guide = F) +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 0.2)) +
  ggtitle(expression("Proportion of " * Delta * " SNP frequency outside cutoff")) +
  xlab("Position (Mb)") +
  ylab(expression("Proportion of extreme " * Delta * " SNP frequency")) # + # +  ylim(c(0, 0.2))
ggsave(filename = "figures/bsa/plot_proportion_delta_allele_frequency_cutoffs5.pdf", width = 10, height = 15)

ggplot(data = delta_df_wnd_chr[delta_df_wnd_chr$cutoff == "2 %", ],
       aes(position/1000000,
           frequency,
           colour = cutoff, alpha = .5)) +
  geom_line() +
  geom_point(alpha = 0.2) +
  facet_wrap(~chromosome, nrow = 7) +
  scale_color_brewer(palette = "Dark2",
                     name = "Cutoff") +
  scale_alpha(guide = F) +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size = 0.2)) +
  ggtitle(expression("Proportion of " * Delta * " SNP frequency outside cutoff")) +
  xlab("Position (Mb)") +
  ylab(expression("Proportion of extreme " * Delta * " SNP frequency")) # + # +  ylim(c(0, 0.2))
ggsave(filename = "figures/bsa/plot_proportion_delta_allele_frequency_cutoffs2.pdf", width = 10, height = 15)



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
#   [1] ggplot2_3.2.0 tidyr_0.8.3   dplyr_0.8.3   vcfUtils_1.0 
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.2         rstudioapi_0.10    magrittr_1.5       tidyselect_0.2.5   munsell_0.5.0      colorspace_1.4-1  
# [7] R6_2.4.0           rlang_0.4.0        tools_3.3.3        grid_3.3.3         gtable_0.3.0       withr_2.1.2       
# [13] digest_0.6.20      yaml_2.2.0         lazyeval_0.2.2     assertthat_0.2.1   tibble_2.1.3       crayon_1.3.4      
# [19] RColorBrewer_1.1-2 purrr_0.3.2        glue_1.3.1         labeling_0.3       stringi_1.4.3      pillar_1.4.2      
# [25] scales_1.0.0       pkgconfig_2.0.2   
