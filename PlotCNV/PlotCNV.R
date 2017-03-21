#!/usr/bin/env Rscript

##########################################################################################
# MSKCC BergerLab
# Plot CNV
# author: Ronak H Shah
# Date: June 20 2016
# Input: Probes.txt
# Output: Image PDF
##########################################################################################


process_probes <- function(probes, pdfName) {
  pdata <- read.delim(probes)
  chr <-
    do.call('c', lapply(as.character(pdata[, 'region']), function(x) {
      f <- unlist(strsplit(x, '\\:'))
      
      return(f[1])
      
    }))
  
  chr <- factor(chr, levels = c(seq(1, 22, 1)))
  chr.counts <- table(chr)
  
  chr.midpt <- sapply(chr.counts, function(x) {
    return(round(x / 2, 0))
  })
  
  numberpos <-
    chr.midpt + c(0, cumsum(chr.counts)[-length(chr.counts)])
  
  linepos <- cumsum(chr.counts)[-length(chr.counts)]
  
  linepos <- as.list(as.data.frame.character(linepos))
  linepos <- levels(droplevels(linepos$linepos))
  #linepos<-c("480","813","1201","1410","1733","1958","2219","2290","2546","2646","2931","3253","3424","3505","3630","3872","4255","4322","4729","4866","4917")
  linepos <- as.numeric(linepos)
  theme_mine <- function(base_size = 12,
                         base_family = "") {
    # Starts with theme_grey and then modify some parts
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
      theme(
        axis.text.x = element_text(
          size = 14,
          hjust = 2,
          angle = 0,
          face = "bold"
        ),
        axis.text.y = element_text(
          size = 14,
          hjust = 1,
          face = "bold"
        ),
        axis.ticks.x =  element_blank(),
        axis.ticks.y =  element_blank(),
        axis.title.x =  element_blank(),
        axis.title.y = element_text(
          size = 16,
          angle = 90,
          face = "bold"
        ),
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major.x = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major.y = element_line(colour = "grey80")
      )
  }
  ggplot(pdata, aes(x = seq(1, nrow(pdata)), y = lr)) + theme_mine() +
    ylab("LOG2 Tumor/Normal Ratio") +
    geom_jitter(aes(color = factor(sig))) +
    scale_y_continuous(breaks = seq(-4, 4, 1), limits = c(-4, 4)) +
    scale_x_continuous(breaks = linepos, labels = seq(1, 21, 1)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(
      values = c("#0072B2", "#D55E00"),
      labels = c("Insignificant", "Significant"),
      guide = guide_legend(
        title = "",
        label.theme = element_text(
          angle = 0,
          size = 18,
          face = "bold"
        )
      )
    )
  
  ggsave(pdfName, width = 20, height = 5,device="pdf")

}
if (!interactive()) {
  pkgs = c(
    'MASS',
    'plyr',
    'reshape2',
    'ggplot2',
    'RColorBrewer',
    'grid',
    'gridExtra',
    'scales',
    'argparse'
  )
  junk <-
    lapply(pkgs, function(p) {
      suppressPackageStartupMessages(require(p, character.only = T))
    })
  rm(junk)
  currwd <- getwd()
  parser = ArgumentParser()
  parser$add_argument('-p', '--probes', type = 'character', help = 'CNV Probes file')
  parser$add_argument('-o', '--outfile', type = 'character', help = 'Output PDF file')
  args = parser$parse_args()
  
  attach(process_probes(args$probes, args$outfile))
}