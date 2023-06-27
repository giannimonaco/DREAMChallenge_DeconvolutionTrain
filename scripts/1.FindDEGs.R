# setwd("PATH/TO/DREAMchallengeDirectory")
library(readr)
library(limma )
library(edgeR)

#### Metadata
Metadata <- readr::read_csv("MedatadaDream.csv")
# if there are cancer samples, just label them as "cancer"
Metadata$CoarseLabels[ Metadata$CoarseLabels %in% c("CRC", "BRCA")] <- "cancer"
Metadata$FineLabels[ Metadata$FineLabels %in% c("CRC", "BRCA")] <- "cancer"


#### Use this for counts
Dataset <- readr::read_csv("KallistoCounts.csv")

#### Use this for TPM
# Dataset <- readr::read_csv("KallistoTPM.csv")

### Use the gene names as rownames. Make sure there are not duplicates
Dataset <- as.data.frame(Dataset)
rownames(Dataset) <- Dataset$Genes
Dataset$Genes <- NULL

#### Use this for the Coarse challenge
MetadataCoarse <- Metadata[!is.na(Metadata$CoarseLabels),]
DatasetFilt <- Dataset[, MetadataCoarse$Sample]
colnames(DatasetFilt) <- MetadataCoarse$CoarseLabels
namenorm <- "DEGs/DEGsCoarse" 

#### Use this for the Fine challenge
# MetadataFine <- Metadata[!is.na(Metadata$FineLabels),]
# DatasetFilt <- Dataset[, MetadataFine$Sample]
# colnames(DatasetFilt) <- MetadataFine$FineLabels
# namenorm <- "../DEGs/DEGsFine"

dir.create("DEGs")
dir.create(namenorm)

# make sure that the data are in counts per million
DatasetNorm <- apply(DatasetFilt, 2, function(x) x/sum(x) * 10^6)


Dataset = DatasetNorm[apply(DatasetNorm, 1, function(x) length(x[x>4])>=3), ] 
cpm1 <- DGEList(counts = Dataset, group=factor(colnames(Dataset)))
allSamples <- colnames(Dataset)
GRall <- unique(colnames(Dataset))


DE1 <- data.frame(genes=rownames(cpm1))
rownames(DE1) <- rownames(cpm1)
DE1$genes <- NULL
DE2 <- DE1
DE3 <- DE1
medGE <- data.frame(genes=rownames(Dataset))
rownames(medGE) <- rownames(Dataset)
medGE$genes <- NULL
for(i in GRall){
  idGR<- which(allSamples == i)
  cpm1$samples$group <- "rest"
  cpm1$samples$group[idGR] <- i
  cpm1$samples$group <- factor(cpm1$samples$group, levels= c( "rest", i))
  
  design.mat <- model.matrix(~ 0 + cpm1$samples$group)
  colnames(design.mat) <- levels(cpm1$samples$group)
  v <- voom(cpm1, design.mat, plot = FALSE)
  
  fit = lmFit(v, design.mat)
  fit = contrasts.fit(fit, c(-1, 1))
  fit = eBayes(fit)
  
  #####FROM HERE, select DEG only with adj pvalue<= 0.05
  # just to verify the number of differentially expressed genes:
  DEdt <- decideTests(fit, adjust.method="BH", p.value=0.05)
  
  # save differentially expressed genes with adj.pvalue less than 0.05:
  DEpval <- topTable(fit, sort.by = "p", p.value=0.05, n=nrow(fit))

  
  DEfc <- topTable(fit, sort.by = "p", p.value=1, n=nrow(fit)) # order by pvalue
  # DEfc <- round(DEfc,3)
  DEfcPOS <- DEfc[DEfc[,"logFC"] >0, "adj.P.Val", drop=F]  # now select only positive logFC
  DEfcNEG <- DEfc[DEfc[,"logFC"] <0, "adj.P.Val", drop=F] * -1  # now select only negative logFC
  DEfcTOT <- rbind(DEfcPOS[, ,drop=F], DEfcNEG[, ,drop=F]) # join them
  DE2 <- cbind(DE2, rep(0, nrow(DE2)))  # create new column with zeros
  colnames(DE2)[ncol(DE2)] <- i          # give it a name  
  DE2[rownames(DEfcTOT), i] <- DEfcTOT[,"adj.P.Val"]    # fill it with signed p-value  
  
  DE3 <- cbind(DE3, rep(0, nrow(DE3)))  # create new column with zeros
  colnames(DE3)[ncol(DE3)] <- i          # give it a name  
  DE3[rownames(DEfc), i] <- DEfc[,"logFC"]    # fill it with logFC
  
  med <- apply(Dataset[, idGR], 1, median)
  medGE <- cbind(medGE, med[rownames(medGE)])
  colnames(medGE)[ncol(medGE)] <- i
  print(paste("done with", i))
}


write.table(DE2, paste0(namenorm, "/DEdecFDR.txt"), sep = "\t", col.names=NA)
write.table(DE3, paste0(namenorm, "/DEdecFC.txt"), sep = "\t", col.names=NA)
write.table(medGE, paste0(namenorm, "/cpm_median.txt"), sep = "\t", col.names=NA)


