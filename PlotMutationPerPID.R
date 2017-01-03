library(dplyr)
library(ggplot2)
library(BBmisc)
library("plyr")
library("reshape2")
library("RColorBrewer")
library("grid")
library("gridExtra")
library("scales")
library(beyonce)
theme_mine <- function(base_size = 12, base_family = "") {
  # Starts with theme_grey and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text.x = element_text(size=10,hjust=0,vjust=0.5,angle=330,face="bold"),
      axis.text.y = element_text(size=14,hjust = 1,face="bold"),
      axis.ticks.x =  element_blank(),
      axis.ticks.y =  element_blank(),
      axis.title.x =  element_blank(),
      axis.title.y = element_blank(),
      legend.position = "top",
      legend.direction="horizontal"
    )
}
makePlot1 <- function(dataDF,outfile,dflen){
  pal <- beyonce_palette(90, type = "continuous")
  #print (dataDF)
  print(outfile)
  datm=unique(melt(dataDF,id.vars = c("Gene_aa","Chr","Pos","Ref","Alt")))
  #print(datm)
  datm$Gene_aa = factor(datm$Gene_aa,levels=dataDF$Gene_aa)
  datm$variable = factor(datm$variable,levels= unique(datm$variable))
  datm$value = as.numeric(datm$value)
  base_size <- 9
  p <- ggplot(datm,aes(x=factor(Gene_aa), y=variable))
  p + geom_tile(aes(fill=value)) + geom_point(aes(size=value),colour = "black",show_guide=TRUE) + 
    scale_fill_gradient2(limits=c(0,1),
                         low=muted("red"),mid="white",high="steelblue",
                         guide = guide_colorbar(title="Variant Allele Frequency (VAF): ", 
                                                title.theme = element_text(angle=0,size=16,face="bold"),
                                                title.position = "top",
                                                label.theme = element_text(angle=0,size=12,face="bold"))) +
    xlab("") + ylab("") +
    coord_equal() +
    scale_size_continuous(breaks = seq(0,1,0.05),limits = c(0,1)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(limit=(rev(levels(datm$variable))),expand = c(0, 0)) +
    theme_mine() 
  
  ggsave(outfile,width=25, height=15)
  
  
}


makePlot <- function(dataDF,outfile,dflen){
  pal <- beyonce_palette(90, type = "continuous")
  #print (dataDF)
  print(outfile)
  datm=unique(melt(dataDF,id.vars = c("Gene_aa","Chr","Pos","Ref","Alt")))
  #print(datm)
  if(dflen<=10){
    fsize=10
  }
  else if(dflen>11 & dflen<=20){
    fsize=8
  }
  else{
    fsize=1
  }
  datm$Gene_aa = factor(datm$Gene_aa,levels=dataDF$Gene_aa)
  datm$variable = factor(datm$variable,levels= unique(datm$variable))
  datm$value = as.numeric(datm$value)
  base_size <- 9
  p <- ggplot(datm,aes(x=factor(Gene_aa), y=variable))
  p + geom_tile(aes(fill=value)) + geom_tile(aes(fill=value),colour = "black",show_guide=FALSE) + 
    geom_text(aes(label = round(value, 3)),size=fsize) +
    scale_fill_gradient2(limits=c(0,1),
                         low=muted("red"),mid="white",high="steelblue",
                         guide = guide_colorbar(title="Variant Allele Frequency (VAF): ", 
                                                title.theme = element_text(angle=0,size=16,face="bold"),
                                                title.position = "top",
                                                label.theme = element_text(angle=0,size=12,face="bold"))) +
    xlab("") + ylab("") +
    coord_equal() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(limit=(rev(levels(datm$variable))),expand = c(0, 0)) +
    theme_mine() 
  
  ggsave(outfile,width=25, height=15)
  
  
}


outDir = getwd()
dataDF = read.delim("../annotated_exonic_variants_r002.txt",header = TRUE,stringsAsFactors=FALSE)
dataDF <- dataDF[dataDF$Call_Confidence  == "HIGH",]
uPID <- unique(dataDF$PID)
#gDF <- group_by(dataDF,PID) 
allOutDF= list()
for (pId in uPID) {
  pData = dataDF[dataDF$PID  == pId,];
  sampleID = sort(unique(pData$Sample));
  print (paste(pId,":",length(sampleID),sep = " "))
  outfile = paste(outDir,paste(pId,"heatmap.pdf",sep="_"),sep = "/")
  if(length(sampleID) == 3) {
    outDF = NaN;
    #print("In 3");
    a<-sampleID[1];
    b<-sampleID[2];
    c<-sampleID[3];
    numrows = nrow(pData);
    outDF = data.frame(Gene_aa = character(numrows),Chr=character(numrows),Pos=integer(numrows),Ref=character(numrows),Alt=character(numrows),a=double(numrows),b=double(numrows), c=double(numrows),stringsAsFactors = FALSE)
    for (i in 1:nrow(pData)) {
      outDF$Gene_aa[i] = paste(pData[i, "Gene"],pData[i,"AAchange"],sep = "_");
      outDF$Chr[i] = as.character(pData[i,"Chrom"]);
      outDF$Pos[i] = pData[i,"Start"];
      outDF$Ref[i] = pData[i,"Ref"];
      outDF$Alt[i] = pData[i,"Alt"];
      outDF$a[i] = getLast(as.list(strsplit(pData[i,a],"="))[[1]]);
      outDF$b[i] = getLast(as.list(strsplit(pData[i,b],"="))[[1]]);
      outDF$c[i] = getLast(as.list(strsplit(pData[i,c],"="))[[1]]);
    }
    colnames(outDF)<-c("Gene_aa","Chr","Pos","Ref","Alt",a,b,c);
    #outDF = outDF[duplicated(outDF), ]
    dflen=length(outDF$Gene_aa)
    print(paste("Length",length(outDF$Gene_aa),sep = ":"))
    makePlot1(outDF,outfile,dflen);
    allOutDF[[pId]] = outDF;
  }
  if(length(sampleID) == 2) {
    outDF = NaN;
    #print("In 2");
    a<-sampleID[1];
    b<-sampleID[2];
    numrows = nrow(pData);
    outDF = data.frame(Gene_aa = character(numrows),Chr=character(numrows),Pos=integer(numrows),Ref=character(numrows),Alt=character(numrows),a=double(numrows),b=double(numrows),stringsAsFactors = FALSE)
    for (i in 1:nrow(pData)) {
      outDF$Gene_aa[i] = paste(pData[i, "Gene"],pData[i,"AAchange"],sep = "_");
      outDF$Chr[i] = as.character(pData[i,"Chrom"]);
      outDF$Pos[i] = pData[i,"Start"];
      outDF$Ref[i] = pData[i,"Ref"];
      outDF$Alt[i] = pData[i,"Alt"];
      outDF$a[i] = getLast(as.list(strsplit(pData[i,a],"="))[[1]]);
      outDF$b[i] = getLast(as.list(strsplit(pData[i,b],"="))[[1]]);
    }
    colnames(outDF)<-c("Gene_aa","Chr","Pos","Ref","Alt",a,b);
    #outDF = outDF[duplicated(outDF), ]
    dflen=length(outDF$Gene_aa)
    print(paste("Length",length(outDF$Gene_aa),sep = ":"))
    makePlot1(outDF,outfile,dflen);
    allOutDF[[pId]] = outDF;
  }
  if(length(sampleID) == 1) {
    outDF = NaN;
    numrows = nrow(pData);
    #print("In 1");
    a<-sampleID[1];
    outDF = data.frame(Gene_aa = character(numrows),Chr=character(numrows),Pos=integer(numrows),Ref=character(numrows),Alt=character(numrows),a=double(numrows),stringsAsFactors = FALSE)
    for (i in 1:nrow(pData)) {
      outDF$Gene_aa[i] = paste(pData[i,"Gene"],pData[i,"AAchange"],sep = "_");
      outDF$Chr[i] = as.character(pData[i,"Chrom"]);
      outDF$Pos[i] = pData[i,"Start"];
      outDF$Ref[i] = pData[i,"Ref"];
      outDF$Alt[i] = pData[i,"Alt"];
      outDF$a[i] = pData[i,"T_AltFreq"]
    }
    colnames(outDF)<-c("Gene_aa","Chr","Pos","Ref","Alt",a);
    #outDF = outDF[duplicated(outDF), ]
    dflen=length(outDF$Gene_aa)
    print(paste("Length",length(outDF$Gene_aa),sep = ":"))
    makePlot1(outDF,outfile,dflen);
    allOutDF[[pId]] = outDF;
  }
  
}
