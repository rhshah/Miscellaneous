#!/usr/bin/env Rscript

##########################################################################################
# MSKCC CMO
# Decode Mutation Signature
# https://github.com/raerose01/deconstructSigs
# author: Ronak H Shah
# Date: March 9 2016
# Input: MAF
# Output: SignatureDecomposition Matrix & PDF
##########################################################################################


#Read Maf File and Retunr only SNPS

process_maf <- function(maf) {
  maf <- read.delim(maf)
  maf <- maf[maf$Variant_Type == "SNP", ]
  prefix <- "chr"
  inputData <- data.frame(Sample=maf$Tumor_Sample_Barcode,
                         chr=paste(prefix,maf$Chromosome,sep=""),
                         pos=maf$Start_Position,
                         ref=maf$Reference_Allele,
                         alt=maf$Tumor_Seq_Allele2)

  sigs <- mut.to.sigs.input(mut.ref = inputData,
                                  sample.id = "Sample",
                                  chr = "chr",
                                  pos = "pos",
                                  ref = "ref",
                                  alt = "alt")
  sampleList = unique(inputData$Sample)
  return(list(sigs.input = sigs,listOfSamples = sampleList))
}

#Run deconstructSigs
run_deconstructSig <- function(sigs.input, referenceSig, sampleID, contexts, method) {
  outData = whichSignatures(tumor.ref = sigs.input, 
                            signatures.ref = referenceSig, 
                            sample.id = sampleID, 
                            contexts.needed = contexts,
                            tri.counts.method = method)
  return(outData)
}

#Plot overall Signature Weights
plot_weights<-function(weightMatrix,pdfFile){
  theme_mine <- function(base_size = 12, base_family = "") {
    # Starts with theme_grey and then modify some parts
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
      theme(
        strip.text.x = element_text(size=18),
        strip.text.y = element_blank(),
        strip.background = element_rect(colour="black", fill="lightgray"),
        axis.text.x = element_text(size=12,angle=90),
        axis.text.y = element_text(size=14,hjust=1,face="bold"),
        axis.ticks.x =  element_line(colour = "black"), 
        axis.ticks.y =  element_line(colour = "black"), 
        axis.title.x= element_text(size=16,face="bold"),
        axis.title.y= element_text(size=16,angle=90,face="bold"),
        legend.position = "right"
      )
  }
  mweightMatrix = melt(as.matrix(weightMatrix))
  newcolnames = c("Sample", "SignatureType", "Weights")
  colnames(mweightMatrix) <- NULL
  colnames(mweightMatrix) <- newcolnames
  mweightMatrix = mweightMatrix[order(mweightMatrix$Weights,decreasing =TRUE),]
  p = ggplot(mweightMatrix,aes(Sample,Weights)) + 
    geom_bar(stat = "identity",aes(fill = SignatureType)) + 
    scale_fill_manual(values = colorRampPalette(c("chocolate4","yellow","#D55E00","hotpink", "#0072B2","lightgreen","darkorchid4"))(length(unique(mweightMatrix$Sample)))) +
    theme_mine()
  ggsave(pdfFile,width=20, height=10)
}
if (!interactive()) {
  
  pkgs = c('argparse', 'deconstructSigs', 'reshape2', 'ggplot2')
  junk <- lapply(pkgs, function(p){suppressPackageStartupMessages(require(p, character.only = T))})
  rm(junk)
  currwd <- getwd()
  parser=ArgumentParser()
  parser$add_argument('-m', '--maf', type = 'character', help = 'SOMATIC_FACETS.vep.maf file')
  parser$add_argument('-o', '--outfileprefix', type = 'character', help = 'Output file prefix')
  parser$add_argument('-s', '--signatureType', type = 'character', help = 'Can be either cosmic or nature2013 default: nature2013', default='nature2013')
  args=parser$parse_args()
  signatureType = ''
  if(args$signatureType == "cosmic"){
    signatureType = signatures.cosmic
  }
  else{
    signatureType = signatures.nature2013
  }
  
  outfileprefix <- args$outfileprefix
  attach(process_maf(args$maf))
  
  listOfSamples <- sort(listOfSamples)
  txtOut = paste(currwd,"/",outfileprefix,"_deconstructSigs_summary.txt",sep="")
  pdf1Out = paste(currwd,"/",outfileprefix,"_deconstructSigs_summary.pdf",sep="")
  pdf2Out = paste(currwd,"/",outfileprefix,"_deconstructSigs.pdf",sep="")
  pdf(file=pdf2Out,width=20,height=10)  
  decomposeSignatureDF = data.frame()
  count = 0
  for (sample in listOfSamples){
    outputdata = run_deconstructSig(sigs.input, signatureType, sample, TRUE, 'default')
    
    if(count == 0){
      decomposeSignatureDF = outputdata$weights
    }
    else{
      decomposeSignatureDF = rbind(decomposeSignatureDF,outputdata$weights)
      }
    par(mfrow = c(2,1))
    plotSignatures(outputdata)
    count = count + 1
  }
  dev.off()
  decomposeSignatureDF = transform(decomposeSignatureDF, Other=(1 -rowSums(decomposeSignatureDF)))
  decomposeSignatureDFOutPut = (as.matrix(decomposeSignatureDF))
  write.table(decomposeSignatureDFOutPut, txtOut, sep = "\t", quote=FALSE, col.names=NA)
  plot_weights(decomposeSignatureDF, pdf1Out)
  
}