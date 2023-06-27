
library(dplyr)
library(readr)
library(limma )
library(edgeR)
library(reshape2)

normAbund <- read.table("CellType_RNAabundance.txt",row.names = 1)


#### set for either the coarse or fine challenge
# setwd("DEGs/DEGsCoarse")
setwd("DEGs/DEGsFine")

DEFDR <- read.table( "DEdecFDR.txt", sep="\t", header = T, check.names = FALSE, stringsAsFactors = F,row.names = 1)
DEFDR1 <- data.frame(gene_name = rownames(DEFDR), DEFDR) 
DEFC <- read.table("DEdecFC.txt", sep="\t", header = T, check.names = FALSE, stringsAsFactors = F,row.names = 1)
DEFC <- DEFC[rownames(DEFDR) ,]
DEFC1 <- data.frame(gene_name = rownames(DEFC), DEFC) 
medianGE <- read.table( "cpm_median.txt", sep="\t", header=T, row.names = 1,check.names = FALSE) 
# medianGE <- medianGE[rownames(DEFDR), ]

SIGMAT <- medianGE

DEFDRmelt <- melt(DEFDR1, id.vars = "gene_name")
DEFCmelt <- melt(DEFC1, id.vars = "gene_name")

### FILTER ON FOLD CHANGE 2 > and FDR < 0.05
# arrange them and filter the positive fold change and positive FDR with value less 0.05
FC1 <- DEFCmelt  %>% dplyr:::filter(value > 2 ) %>% arrange(dplyr:::desc(value)) # using the fold change of 2 seems to improve the results respect to smaller fold changes
FDR1 <- DEFDRmelt  %>% dplyr:::filter(value >= 0 & value < 0.05 ) %>% arrange(value) # if you change the threshowld with 0.2 you get worse results
# remove duplicate genes with lower fold change in both FC1 and FDR1
FC1 <- FC1[!duplicated(FC1$gene_name),]
combFC1 <- paste(FC1$gene_name, FC1$variable)
combFDR1 <- paste(FDR1$gene_name, FDR1$variable)
FDR2 <- FDR1[combFDR1 %in% combFC1, ]
FDR2 <- FDR2 %>% group_by(variable) # %>% slice(1:50) 


####################################################################
###############   FILTERS SIGNATURE MATRIX ###########################
### FILTER 1  (HIGH EXPRESSION)
# filter out those genes with very high expression 
filt1 <- apply(SIGMAT, 1, function(x) !any(x > 3000))  # this filtering improves the estimation quite a lot too
SIGMAT1 <- SIGMAT[filt1,]


### FILTER 2 (LOW EXPRESSION)
# filter out those genes with very low expression in all cell type
filt2 <- rowSums(SIGMAT1) > 5  
SIGMAT2 <- SIGMAT1[filt2,]


##### FILTER 3 (Effect size)
SIGMATtemp <- t(apply(SIGMAT2, 1, sort, decreasing =T))
filt3 <- (log2(SIGMATtemp[,1]) - log2(SIGMATtemp[,2])) > 0.1
SIGMAT3 <- SIGMAT2[filt3, ]


FDR3 <- FDR2[FDR2$gene_name %in% rownames(SIGMAT3), ] 

dat <- FDR3 %>% group_by(variable)

##################### SIGNATURE MATRIX
##### create the signature and mixed matrix (create a smaller matrix for the coarse challenge, and a bigger matrix for the fine challenge)
if(grepl("Coarse", getwd())){SliceN <- 70} else { SliceN <- 90 }
dat2 <- dat %>% group_by(variable) %>% dplyr:::slice(1:SliceN)     # nKappa)
genesDec <- unique(as.character(dat2$gene_name))
# SigMat1 <- SIGMAT[genesDec , ]


##################### Normalization for RNA abundance 
## for cancer I just use 1, hence I don't normalize; for other cell type I use the mean value from cell of the same lineage
if(grepl("Coarse", getwd())){
  medianGE <- medianGE[, c("cancer","endothelial.cells","fibroblasts","B.cells","CD4.T.cells","CD8.T.cells",
                       "NK.cells","monocytic.lineage","neutrophils" )]

  normAbund2 <- c(1,1,1,
                  mean(normAbund[c("B Naive","B Memory"), ]),
                  mean(normAbund[c("T CD4 Naive","T CD4 Memory"), ]),
                  mean(normAbund[c("T CD8 Naive","T CD8 Memory"), ]),
                  normAbund["NK", ],
                  normAbund["Monocytes C", ],
                  normAbund["Neutrophils LD", ] )
  medianGEnorm<- sweep(medianGE, MARGIN=2, normAbund2, `*`)
  SigMatCoarse <- data.frame( medianGEnorm[genesDec, ])
  save(SigMatCoarse, file="../../SigMatCoarse.RData")
}else{
    medianGE <- medianGE[, c("cancer","endothelial.cells","fibroblasts",
                         "naive.B.cells","memory.B.cells",
                         "naive.CD4.T.cells", "memory.CD4.T.cells","regulatory.T.cells",
                         "naive.CD8.T.cells","memory.CD8.T.cells", "NK.cells",
                         "monocytes","myeloid.dendritic.cells","macrophages","neutrophils")]
  normAbund2 <- c(1,1,1,
                  normAbund["B Naive", ],
                  normAbund["B Memory", ],
                  normAbund["T CD4 Naive", ],
                  normAbund["T CD4 Memory", ],
                  normAbund["T CD4 Naive", ],
                  normAbund["T CD8 Naive", ],
                  normAbund["T CD8 Memory", ],
                  normAbund["NK", ],
                  normAbund["Monocytes C", ],
                  normAbund["mDCs", ],
                  normAbund["Monocytes C", ],
                  normAbund["Neutrophils LD", ])
  medianGEnorm<- sweep(medianGE, MARGIN=2, normAbund2, `*`)
  SigMatFine <- data.frame(medianGEnorm[genesDec, ])
  save(SigMatFine, file="../../SigMatFine.RData")
}



