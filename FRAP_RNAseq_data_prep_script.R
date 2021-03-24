## Script to organize and plot DESEQ2 output files from RNA-seq experiments
## Generate log2 ratios
## Input: deseq.txt file
## Functions to clean and organize RNA-seq data to plot log2 ratios
## Requires openxlsx if using an annotation file that is in xlsx format

library("openxlsx")

## Function to match Genes to Wormbase IDs for each respective gene
getGenesWbId <- function(c_elegans_annots, genes){
  gene.to.wbid<-read.table(file=c_elegans_annots,header=F,stringsAsFactors=F)
  colnames(gene.to.wbid)<-c('gene','wbid')
  relevant_genes_ix <- match(genes, gene.to.wbid$gene)
  wbid<-gene.to.wbid$wbid[relevant_genes_ix]
}
## Function used to exclude duplications that are in a provided annotation file and keep relevant genes outside of these regions
getCelegansAnnotations <- function(file){
  c_elegans_annots <- read.xlsx(file)
  relevant_cols <- c("Gene.WB.ID", "Sequence.Name.(Gene)", "Chr.Name")
  c_elegans_annots <- c_elegans_annots[,relevant_cols]
  c_elegans_annots <- c_elegans_annots[!duplicated(c_elegans_annots[,"Sequence.Name.(Gene)"]),]
  c_elegans_annots <- c_elegans_annots[complete.cases(c_elegans_annots[,relevant_cols]),]
  return(c_elegans_annots)
}

c_elegans_annots <- getCelegansAnnotations("~/Desktop/Data/RNAseq/WS220_gene_annotations_WormMart_unique.xlsx")
annot<-read.table("~/Desktop/Data/RNAseq/c.elegans.WS220.gene.wbid.chr.start.stop.txt",header=T)
tannot<-annot
tannot$WBID<-as.character(tannot$WBID)
tannot$chr<-as.character(tannot$chr)
g1<-c("WBGene00044083","B0564.1","IV")
g2<-c("WBGene00002977","Y105E8A.7","I")
comp_annot<-rbind(c_elegans_annots,g1,g2)

## Function to append chromosome name to the dataset
get.chr<-function(dataset){
  for(i in 1:length(dataset$wbid)){
    chr.name[i]<-comp_annot$Chr.Name[which(comp_annot$Gene.WB.ID==dataset$wbid[i])]
  }
  output<-cbind(dataset,chr.name)
  return(output)
}
## Function to remove mitochondrial DNA and append start and stop positions for each gene in the dataset
get.pos<-function(dataset){
  dataset<-dataset[which(dataset$chr.name!="MtDNA"),] #removes MtDNA
  chr.start<-""
  chr.stop<-""
  for(i in 1:length(dataset$wbid)){
    chr.start[i]<-tannot$start[which(tannot$WBID==dataset$wbid[i])]
    chr.stop[i]<-tannot$stop[which(tannot$WBID==dataset$wbid[i])]
  }
  output<-cbind(dataset,chr.start,chr.stop)
  return(output)
}