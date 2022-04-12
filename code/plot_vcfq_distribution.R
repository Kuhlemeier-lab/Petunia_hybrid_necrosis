#!/usr/bin/env Rscript

####    Plot density distribution of a set of quality parameters    ####

## See README.md for detailed info on usage.
## 03 August 2018
## Licence: GNU GPL-3.0 or later
## Copyright (C) 2018  Marta Binaghi <marta.binaghi at ips.unibe.ch>

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.


# Libraries ---------------------------------------------------------------
library("optparse")

# get arguments -----------------------------------------------------------
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--output"), type = "character", default = "./",
              help = "path for the output folder, [default = %default]", metavar = "character"),
  make_option(c("--vcftype"), type = "character", default = "snp",
              help = "Type of variants in the table (snp or indel) [default = %default]",
              metavar = "character"),
  make_option(c("-p", "--parameter"), type="character", default="all",
              help="parameter to plot [default = %default]", metavar="character"),
  make_option(c("--pdf"), type="character", default = "onepage",
              help = "all plots in one page pdf or single plot files (singleplots) [default = %default]",
              metavar = "character")
);

opt_parser = OptionParser(usage = "Usage: %prog -i myvariants.table",
                          option_list=option_list);
opt = parse_args(opt_parser);

# manage NULL input -------------------------------------------------------

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

# check output directory --------------------------------------------------
# check if path ends with a slash "/", add it otherwise
if ( substr(opt$output, start = nchar(opt$output), stop = nchar(opt$output)) != "/") {
  opt$output <- paste0(opt$output, "/")
}
# check if output directory exists, creates it otherwise
if (!dir.exists(opt$output)) {
  dir.create(opt$output,
             recursive = T)
  writeLines(paste0("\nThe directory ", opt$output, " was created to store plot output."))
  } else {
  writeLines(paste0("\nThe output will be saved in ", opt$output))
  }

# import file -------------------------------------------------------------
writeLines(paste0("Reading dataset file '",
             opt$input, "'.",
             "\nThis might take a while..."))

df <- read.table(opt$input,
                 sep = "\t",
                 dec = ".",
                 header = T)

writeLines("...reading done.\n")

# set parameters ----------------------------------------------------------
# store common ranges for each default parameter
# parameters differ if input file contains SNPs or INDELs
opt$vcftype <- tolower(opt$vcftype)
if ( opt$vcftype == "snp" ) {
  ranges <- rbind(min = c(0, 0, 0, 0, -13, -10),
                max = c(60, 50, 12, 100, 8, 10))
  colnames(ranges) <- c("QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum")
  # GATK recommended thresholds
  ranges <- rbind(ranges,
                c(2, 60, 3, 40, -12.5, -8))
} else if ( opt$vcftype == "indel" ) {
  ranges <- rbind(min = c(0, 0, 0, -25), ## to be fixed
                  max = c(60, 50, 12, 10)) ## to be fixed
  colnames(ranges) <- c("QD", "FS", "SOR", "ReadPosRankSum")
  # GATK recommended thresholds
  ranges <- rbind(ranges,
                  c(2, 200, 10, -20))
  } else {
    stop("vcftype must be 'snp' or 'indel'.\n", call.=FALSE)
  }

# define for which quality parameters the plots are done
# warn for automatic set of --pdf="singleplots" if -par != "all"
if (opt$parameter == "all") {
  if ( opt$vcftype == "snp" ) {
    parameters <- c("QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum")
  } else {
    parameters = c("QD", "FS", "SOR", "ReadPosRankSum")
  }
} else {
  parameters <- opt$parameter
  writeLines("With custom setting of parameters to plot, --pdf is forced to 'singleplots'.")
  opt$pdf <- "singleplots"
}

writeLines("The plots will be produced for the quality parameters:")
print(parameters)
# writeLines("\n")

# check if these parameters are really in the table
if (!all(parameters %in% colnames(df))) {
  stop("The input table must have a column for each parameter to plot.\nSee README.md for further details.", call.=FALSE)
  }

# plot function -----------------------------------------------------------
mydensityplot <- function(x_values,
                          main_title = "Density distribution",
                          par_values = NULL) {
  plot(density(x_values,
               na.rm = T),
       main = main_title,
       xlab = "Parameter value",
       bty = "n",
       xaxt = "n",
       zero.line = F,
       xlim = par_values[1:2])
  usr <- par("usr")
  axis(1,
       col = "azure3",
       col.axis="azure3")
  axis(1, labels = F,
       at = round(boxplot.stats(x_values)$stats, digits = 2),
       col = "black",
       las = 2)
  abline(v = par_values[3],
         col="azure3", lwd=.7, lty=1)
  if (!is.null(par_values)) {
    text(x = par_values[3],
         y = usr[4],
         adj = c(-0.1, 1.1),
         cex = .8,
         col = "azure3",
         labels = "GATK hard-filter")
  }
}
#mydensityplot(df$FS, par_values = c(0, 100, 60))

# make plots --------------------------------------------------------------
# check --pdf has a valid value
if (opt$pdf != "singleplots" & opt$pdf != "onepage") {
    stop("--pdf option is invalid. Use 'onepage' (default) or 'singleplots'.\n", call.=TRUE)
  } else if (opt$pdf == "singleplots") {  # each plot saved in a different file
    for (i in parameters) {
      filename <- paste0(opt$output, opt$vcftype, "_", i, ".pdf")
      plot_title <- paste0("Density distribution of ", i)
      cairo_pdf(file = filename,
                width = 6, height = 4.5)
      if (i %in% colnames(ranges)) {
          par_values <- ranges[ , i]
        } else {
          par_values <- c(NULL)
        }
      mydensityplot(df[ , i],
                      main_title = plot_title,
                      par_values = par_values)
      dev.off()
    }
  } else if (opt$pdf == "onepage") { # all plots in one page
    filename <- paste0(opt$output, opt$vcftype, "_quality.pdf")
    cairo_pdf(file = filename,
              onefile = TRUE,
              width = 8.27,
              height = 11.69
              )
    par(mfrow=c(3, 2))
    for (i in parameters) {
      plot_title <- paste0("Density distribution of ", i)
      if (i %in% colnames(ranges)) {
        par_values <- ranges[ , i]
      } else {
        par_values <- c(NULL)
      }
      mydensityplot(df[ , i],
                    main_title = plot_title,
                    par_values = par_values)
     }
     dev.off()
  }

writeLines("\nPlots are plotted!
See README.md for further details on the plots and hard-filtering
recommended values.")
