
###analysis of human parotid gland scRNA-seq

library(Seurat)
library(cowplot)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)


Normal.data<-Read10X_h5("D://Bioinfo Data/HT2020-22277-1_QC_report20210127/filtered_feature_bc_matrix.h5")
Normal <- CreateSeuratObject(counts = Normal.data, project = "Normal")
Normal[["percent.mt"]] <- PercentageFeatureSet(Normal, pattern = "^MT-")
VlnPlot(Normal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.00001)


Normal <- subset(Normal, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 25 & nCount_RNA<20000)
Normal <- NormalizeData(Normal, normalization.method = "LogNormalize", scale.factor = 10000)
Normal<- FindVariableFeatures(Normal, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Normal)
Normal<- ScaleData(Normal, features = all.genes)

Normal <- RunPCA(Normal)
Normal <- FindNeighbors(Normal, dims = 1:40)
Normal <- FindClusters(Normal, resolution = 1)
Normal <- RunUMAP(Normal, dims = 1:40)


Idents(Normal)<-Normal$integrated_snn_res.1


Normal <- RenameIdents(Normal,"0"="T",`1` = "T", `2` = "T",
                             `3` = "B", `4` = "B", `5` = "T",
                             `6` = "B", `7` = "B", `8` = "T", `9` = "T",
                             `10` = "Serous Acinar Cell", `11` = "B","12"="T",
                       "13"="T","14"="Fibroblast","15"="NK",'16'="Myeloid","17"="T",
                       '18'="Myeloid","19"="Plasma","20"="Ductal Epithelial Cell")

DimPlot(Normal,repel = TRUE,label = TRUE)


Normal.markers <- FindAllMarkers(Normal, only.pos = TRUE,)
Normal.markers <- FindAllMarkers(MinorSG, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
top5 <- Normal.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(Normal, features = top5$gene,label = T)

VlnPlot(Normal, features = c("LPO","KRT19","ACTA2","COL1A1","CDH5","CD3G","KLRF1","CD79A","MZB1","LYZ"), pt.size = 0.000000001, combine = T)

FeaturePlot(Normal, features = c("LPO","KRT19","ACTA2","COL1A1","CDH5","CD3G","KLRF1","CD79A","MZB1","LYZ"))

markers.to.plot <- c("CD3G","LPO","CD79A","KLRF1","COL1A1","MZB1","LYZ","KRT7","CDH5","ACTA2")
p<-DotPlot(Normal, features = markers.to.plot, cols = c("#D5EEF0", "#78A1CA"), dot.scale = 8) 
p
p+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


markers.to.plot <- c( "HLA-DMA","HLA-DMB", "HLA-DQA1","HLA-DRA","HLA-DPA1","HLA-DPB1","PTTG1", "TNPO3","STAT4", "VCAM1" ,"GTF2I")
p<-DotPlot(Normal, features = markers.to.plot, cols = c("#A66FD1", "#FF8451"), dot.scale = 8) 
p
p+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))



####analysis of human minor salivary gland and parotid gland

library(Seurat)
library(cowplot)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)


ParotidGland.data<-Read10X_h5("C:/Users/21578/Desktop/HT2020-22277-1_QC_report20210127/filtered_feature_bc_matrix.h5")
ParotidGland <- CreateSeuratObject(counts = ParotidGland.data, project = "ParotidGland")

MinorSalivaryGland<-readRDS("F:/biodata/warner20_salivary_gland_seurat.rds")

Bothcombined <- merge(ParotidGland, y = c(MinorSalivaryGland), add.cell.ids = c("ParotidGland", "MinorSalivaryGland"), project = "ALL")

Bothcombined[["percent.mt"]] <- PercentageFeatureSet(Bothcombined, pattern = "^MT-")
VlnPlot(Bothcombined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Bothcombined <- subset(Bothcombined, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 25)


Bothcombined.list <- SplitObject(Bothcombined, split.by = "orig.ident")
Bothcombined.list <- lapply(X = Bothcombined.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
anchors <- FindIntegrationAnchors(object.list = Bothcombined.list, dims = 1:30,anchor.features = 2000)
Bothcombined <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(Bothcombined) <- "integrated"
Bothcombined <- ScaleData(Bothcombined, vars.to.regress = c('nCount_RNA'),
                          model.use = 'linear', use.umi = FALSE)
Bothcombined <- RunPCA(object = Bothcombined, pc.genes = VariableFeatures(Bothcombined))
Bothcombined <- FindNeighbors(Bothcombined, dims = 1:20)
Bothcombined <- FindClusters(Bothcombined, resolution = 2)
Bothcombined <- RunUMAP(Bothcombined, dims = 1:40)


DimPlot(Bothcombined, reduction = "umap", split.by = "orig.ident",label = TRUE,repel = TRUE)

Bothcombined <- RenameIdents(Bothcombined,"0"="T","1"="T","2"="T","3"="T","4"="T","5"="Serous Acinar Cell",
                             "6"="T","7"="B","9"="B","10"="NK","11"="B",`12` = "Serous Acinar Cell", 
                             "13"="Fibroblast","14"="T",`15` = "Serous Acinar Cell"
                             ,"16"="Plasma","17"="B","18"="B","19"="Mucinous Acinar Cell"
                             ,"20"="Myeloid","21"="B","22"="Ductal epithelial Cell"
                             ,"23"="Serous Acinar Cell","24"="Serous Acinar Cell"
                             ,"25"="Fibroblast","26"="Mucinous Acinar Cell","27"="Vascular","28"="Serous Acinar Cell"
                             ,"29"="Erythrocyte","30"="Plasma","31"="Myoepithelium","32"="Ductal epithelial Cell","33"="Myeloid")


Idents(Bothcombined)<-Bothcombined$integrated_snn_res.2


Epithelium<-subset(Bothcombined,idents  = c(12,5,15,23,24,26,19,28))
Epithelium <- RenameIdents(Epithelium,"5"="Serous Acinar Cell",`12` = "Serous Acinar Cell", `15` = "Serous Acinar Cell"
                           ,"19"="Mucinous Acinar Cell","22"="Ductal epithelial Cell","23"="Serous Acinar Cell","24"="Serous Acinar Cell","26"="Mucinous Acinar Cell","28"="Serous Acinar Cell","32"="Ductal epithelial Cell")

Epithelium<-subset(Bothcombined,idents  = c("Mucinous Acinar Cell","Serous Acinar Cell"))
Epithelium$celltype<-Idents(Epithelium)


Epithelium <- NormalizeData(Epithelium, normalization.method = "LogNormalize", scale.factor = 10000)
Epithelium<- FindVariableFeatures(Epithelium, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Epithelium)
Epithelium<- ScaleData(Epithelium, features = all.genes)
Epithelium <- RunPCA(Epithelium, features = VariableFeatures(object = Epithelium))

Epithelium <- FindNeighbors(Epithelium, dims = 1:40)
Epithelium <- FindClusters(Epithelium, resolution = 0.5)
Epithelium <- RunUMAP(Epithelium, dims = 1:10)

DimPlot(Epithelium,repel = T)

Idents(Epithelium)<-Epithelium$celltype



FeaturePlot(Epithelium, features = c("HTN1","LTF","AMY1A","MUC16","MUC5B","FUT2"))

Idents(Epithelium) <- Epithelium$celltype
Epithelium.de.markers <- FindMarkers(Epithelium, ident.1 = "Mucinous Acinar Cell", ident.2 = "Serous Acinar Cell")

Epithelium.markers <- FindAllMarkers(Epithelium, only.pos = TRUE, min.pct = -1, logfc.threshold = -1)
Epithelium.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)



# Differential expressed gene analysis

library(DESeq2)
library(limma)
library(pasilla)


database_all <- read.table(file = "C://Users/21578/Downloads/RNA-SEQ GSE159574_Sample_raw_counts_matrix.txt", sep = "\t", header = T,row.names = 1)
database <- database_all[,1:29]

condition <- factor(c(rep("Normal",13), rep("Sjogren",16)))
coldata <- data.frame(row.names = colnames(database), condition)
dds <- DESeqDataSetFromMatrix(countData=database, colData=coldata, design=~condition)

dds <- DESeq(dds)
res <- results(dds)

table(res$padj <1)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)



###Volcano plot

library(ggplot2)
library(ggthemes)
library(ggrepel)
library(dplyr)

data<-Epithelium.de.markers
data$threshold<-as.factor(ifelse(data$adj.P.Val < 0.05 & data$logFC >= 1, 'Up-regulated', ifelse(data$adj.P.Val < 0.05 & data$logFC<= -1, 'Down-regulated', 'Not-significant')))

p<-ggplot(data,aes(x=data$logFC,y=-log10(data$adj.P.Val),colour=threshold))+xlab("log2 Fold Change")+ylab("-log10P-Value")+
  
  geom_point(size=1,alpha=0.6)+
  
  scale_color_manual(values =c("#0072B5","grey","#BC3C28")) +
  
  geom_hline(aes(yintercept=-log10(0.05)),colour="grey",size=1.2 ,linetype=2) + 
  
  geom_vline(aes(xintercept=0), colour="grey",size=1.2 ,linetype=2)+ 
  
  theme_few()+theme(legend.title = element_blank()) 

p


###analysis of human and mouse parotid gland

library(Seurat)
library(cowplot)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(Seurat)
library(tidyverse)

mParotidGland.data<-read.table("C:/Users/21578/Desktop/mouse parotid gland/GSM3895267_PG.P18.filtMatrix.txt",sep = ",",header = T,row.names=1)


library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
class(human)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
class(mouse)


mParotidGland.data<-as.matrix(mParotidGland.data@assays[["RNA"]]@data)
genes = mParotidGland.data$Gene

genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
               values = genes ,mart = mouse,
               attributesL = c("hgnc_symbol"),
               martL = human,
               uniqueRows=T)
gene2<-genes
gene2$Gene<-gene3$Gene.stable.ID
gene2$Gene<-gene3$HGNC.symbol
mParotidGland.data<-merge(gene2,data,by="Gene")

mParotidGland.data<-mParotidGland.data[,c(-1,-2)]

mParotidGland <- CreateSeuratObject(counts = mParotidGland.data, project = "mParotidGland")


hParotidGland.data<-Read10X_h5("C:/Users/21578/Desktop/HT2020-22277-1_QC_report20210127/filtered_feature_bc_matrix.h5")
hParotidGland <- CreateSeuratObject(counts = hParotidGland.data, project = "hParotidGland")

species <- merge(mParotidGland, y = c(hParotidGland), add.cell.ids = c("mParotidGland", "hParotidGland"), project = "ALL")


Bothcombined.list <- SplitObject(Bothcombined, split.by = "orig.ident")
Bothcombined.list <- lapply(X = Bothcombined.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
anchors <- FindIntegrationAnchors(object.list = Bothcombined.list, dims = 1:30,anchor.features = 2000)
Bothcombined <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(Bothcombined) <- "integrated"
Bothcombined <- ScaleData(Bothcombined, vars.to.regress = c('nCount_RNA'),
                          model.use = 'linear', use.umi = FALSE)
Bothcombined <- RunPCA(object = Bothcombined, pc.genes = VariableFeatures(Bothcombined))
Bothcombined <- FindNeighbors(Bothcombined, dims = 1:20)
Bothcombined <- FindClusters(Bothcombined, resolution = 2)
Bothcombined <- RunUMAP(Bothcombined, dims = 1:40)





DimPlot(Bothcombined, reduction = "umap", split.by = "orig.ident",label = T,repel = T)

Bothcombined.markers <- FindAllMarkers(Bothcombined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Bothcombined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)




Bothcombined <- RenameIdents(Bothcombined,"0"="T","2"="T","3"="T","4"="T",
                             "6"="B","7"="T","8"="T","9"="B","10"="B",
                             "13"="B","14"="T",`15` = "T"
                             ,"17"="NK","18"="T","19"="Epithelium"
                             ,"20"="Epithelium","21"="B"
                             ,"23"="Fibroblast","25"="Fibroblast","26"="Plasma"
                             ,"29"="Myeloid","30"="Epithelium","31"="Epithelium","32"="Epithelium","33"="Myeloid",
                             "34"="Myeloid","35"="Vascular")
Bothcombined$celltype<-Idents(Bothcombined)
Bothcombined<-subset(Bothcombined,idents  = c("Fibroblast","Vascular","Myeloid","Epithelium","NK","T","B","Plasma"))



Bothcombined2<-subset(Bothcombined,idents  = c("Fibroblast","Vascular","Epithelium"))
all.genes <- rownames(Bothcombined2)
Bothcombined2<- ScaleData(Bothcombined2, features = all.genes)



Bothcombined2 <- RunPCA(Bothcombined2, features = VariableFeatures(object = Bothcombined2))

VizDimLoadings(Bothcombined2, dims = 1:2, reduction = "pca")
DimPlot(Bothcombined2, reduction = "pca")


Bothcombined2 <- FindNeighbors(Bothcombined2, dims = 1:40)
Bothcombined2 <- FindClusters(Bothcombined2, resolution = 10)
Bothcombined2 <- RunUMAP(Bothcombined2, dims = 1:40)

Idents(Bothcombined2)<-Bothcombined2$integrated_snn_res.10
DimPlot(Bothcombined2,label = TRUE)
Idents(Bothcombined2)<-Bothcombined2$celltype



Bothcombined2.markers <- FindAllMarkers(Bothcombined2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Bothcombined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

Bothcombined2 <- RenameIdents(Bothcombined2,"0"="Serous Acinar Cell","1"="Serous Acinar Cell","2"="Fibroblast","3"="T","4"="Serous Acinar Cell","5"="Serous Acinar Cell",
                              "6"="Serous Acinar Cell","7"="Serous Acinar Cell","8"="T","9"="Fibroblast","10"="Serous Acinar Cell","11"="Serous Acinar Cell","12"="Fibroblast",
                              "13"="T","14"="Ductal epithelial Cell","15" = "T","16"="Serous Acinar Cell"
                              ,"17"="T","18"="Serous Acinar Cell","19"="Serous Acinar Cell"
                              ,"20"="T","21"="T","22"="Ductal epithelial Cell"
                              ,"23"="Serous Acinar Cell","24"="Serous Acinar Cell","25"="Ductal epithelial Cell","26"="Endothelium","27"="Serous Acinar Cell","28"="Fibroblast",
                              "29"="T","30"="Serous Acinar Cell","31"="Myoepithelium","32"="Serous Acinar Cell",
                              "33"="","34"="Serous Acinar Cell","35"="Fibroblast",
                              "36"="Serous Acinar Cell","37"="Serous Acinar Cell","38"="Serous Acinar Cell","39"="Serous Acinar Cell",
                              "40"="Serous Acinar Cell","41"="Serous Acinar Cell","42"="Serous Acinar Cell","43"="Serous Acinar Cell",
                              "44"="Serous Acinar Cell","45"="Serous Acinar Cell","46"="Serous Acinar Cell",
                              "47"="Serous Acinar Cell","48"="Serous Acinar Cell","49"="Serous Acinar Cell","50"="Serous Acinar Cell",
                              "51"="Serous Acinar Cell","52"="Serous Acinar Cell","53"="Serous Acinar Cell","54"="Serous Acinar Cell","55"="Serous Acinar Cell",
                              "56"="Serous Acinar Cell","57"="Fibroblast","58"="Serous Acinar Cell")
Bothcombined2$celltype<-Idents(Bothcombined2)


Bothcombined3<-subset(Bothcombined2,idents  = c("Fibroblast","Endothelium","Serous Acinar Cell","Ductal epithelial Cell","Myoepithelium"))
all.genes <- rownames(Bothcombined3)
Bothcombined3<- ScaleData(Bothcombined3, features = all.genes)
Bothcombined3 <- RunPCA(Bothcombined3, features = VariableFeatures(object = Bothcombined3))
Bothcombined3 <- FindNeighbors(Bothcombined3, dims = 1:4)
Bothcombined3 <- FindClusters(Bothcombined3, resolution = 6)
Bothcombined3 <- RunUMAP(Bothcombined3, dims = 1:4)
Idents(Bothcombined3)<-Bothcombined3$celltype


DimPlot(Bothcombined3,label = TRUE,repel = TRUE,split.by = "orig.ident")



Serous Acinar Cell<-subset(Bothcombined3,idents = c("Serous Acinar Cell"))
Idents(Serous Acinar Cell)<-Serous Acinar Cell$orig.ident
Serous Acinar Cell.de.markers <- FindMarkers(Serous Acinar Cell, ident.1 = "hParotidGland", ident.2 = "mParotidGland")



####integrated analysis of human digestive glands
library(Seurat)
library(cowplot)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)


Pancreas1 <- Read10X(data.dir = "D:/Bioinfo Data/gland single cell data/pancreas/AdjNorm_TISSUE_1/filtered_feature_bc_matrix/")
Pancreas1 <- CreateSeuratObject(counts = Pancreas1, project = "Pancreas1")
orig=c(rep("Pancreas",4969))
Pancreas1@meta.data$orig<-orig

Rectal1<-read.csv("D:/Bioinfo Data/gland single cell data/Rectal single cells GSM4850586_Rectum_Counts.csv/Rectal single cells GSM4850586_Rectum_Counts.csv",header = T)
Rectal1<-CreateSeuratObject(counts = Rectal1,min.cells = 3,project = "Rectal1")
orig=c(rep("Rectal",6280))
Rectal1@meta.data$orig<-orig

Stomach1<-read.csv("D:/Bioinfo Data/gland single cell data/Stomach small amount GSM4850590_Stomach_Counts.csv/Stomach small  GSM4850590_Stomach_Counts.csv",header = T)
Stomach1<-CreateSeuratObject(counts = Stomach1,min.cells = 3,project = "Stomach1")
orig=c(rep("Stomach",5318))
Stomach1@meta.data$orig<-orig

Intestine1<-read.csv("D:/Bioinfo Data/gland single cell data/Small intestine GSM4850588_Small.intestine_Counts.csv/Small intestine  GSM4850588_Small.intestine_Counts.csv",header = T)
Intestine1<-CreateSeuratObject(counts = Intestine1,min.cells = 3,project = "Intestine1")
orig=c(rep("Intestine",4312))
Intestine1@meta.data$orig<-orig

Parotid1<-Read10X_h5("D:/Bioinfo Data/HT2020-22277-1_QC_report20210127/filtered_feature_bc_matrix.h5")
Parotid1<- CreateSeuratObject(counts = Parotid1, project = "Parotid1")
orig=c(rep("Parotid",17736))
Parotid1@meta.data$orig<-orig

Esophagus<-read.csv("D:/Bioinfo Data/gland single cell data/GSM4850580_Esophagus_Counts.csv.gz")
Esophagus1<-CreateSeuratObject(counts = Esophagus,min.cells = 3,project = "Esophagus1")
orig=c(rep("Esophagus",9117))
Esophagus1@meta.data$orig<-orig


###消化腺单细胞合并
AL.combined<-merge (Parotid1,y=c(MinorSalivary,Esophagus1,Stomach1,Pancreas1,Intestine1,Rectal1), 
add.cell.ids = c('Parotid',"SmallSalivary","Esophagus","Stomach","Pancreas","Intestine","Rectal"),
project = "ALL")


AL.combined[["percent.mt"]] <- PercentageFeatureSet(AL.combined, pattern = "^MT-")
VlnPlot(AL.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.000001)
AL.combined <- subset(AL.combined, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 25)
AL.combined.list <- SplitObject(AL.combined, split.by = "orig.ident")
AL.combined.list <- lapply(X = AL.combined.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

anchors <- FindIntegrationAnchors(object.list = AL.combined.list, dims = 1:30)
AL.combined <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(AL.combined) <- "integrated"

AL.combined <- ScaleData(AL.combined, vars.to.regress = c('nCount_RNA'),
                         model.use = 'linear', use.umi = FALSE)
AL.combined <- RunPCA(object = AL.combined, pc.genes = VariableFeatures(AL.combined))
AL.combined <- FindNeighbors(AL.combined, dims = 1:50)
AL.combined <- FindClusters(AL.combined, resolution = 0.3)
AL.combined <- RunUMAP(AL.combined, dims = 1:50)



AL.combined <- RenameIdents(AL.combined,"0"="T","1"="Epithelium","2"="B","3"="NK","4"="Fibroblast","5"="Myeloid",
                            "6"="Epithelium","7"="Smooth Muscle Cell","8"="Plasma","9"="Vasuclar","10"="Epithelium","11"="Epithelium",
                            "12"="Epithelium","13"="Fibroblast","14"="Epithelium","15"="Schwann Cell",
                            "16"="Epithelium")

DimPlot(AL.combined, reduction = "umap",label = TRUE)
DimPlot(AL.combined, reduction = "umap", split.by = "orig.ident",label = T,repel = T)
DimPlot(sce.big, reduction = "umap", split.by = "orig.ident",label = T,repel = T)
p<-DimPlot(AL.combined,label = F,reduction = "umap",repel = T,group.by = "orig.ident")
library(ggsci)
p1<-p + scale_color_simpsons()
p1


Epithelium<-subset(AL.combined,ident=c('Epithelium'))





all.genes <- rownames(Epithelium)
Epithelium<- ScaleData(Epithelium, features = all.genes)

Epithelium <- RunPCA(Epithelium, features = VariableFeatures(object = Epithelium))
Epithelium <- FindNeighbors(Epithelium, dims = 1:30)
Epithelium <- FindClusters(Epithelium, resolution = 0.3)
Epithelium <- RunUMAP(Epithelium, dims = 1:50)

DimPlot(Epithelium,label = TRUE,repel = TRUE)

FeaturePlot(Epithelium, features = c("PIP","MUC1","KRT5","MKI67","ANPEP"))

Epithelium <- RenameIdents(Epithelium,"0"="Mucinous Acinar Cell","1"="Serous Acinar Cell","2"="Serous Acinar Cell","3"="Enterocyte","4"="Basal Epithelial Cell","5"="Serous Acinar Cell",
                       "7"="Squamous Epithelial Cell","8"="Enterocyte","9"="Enterocyte")

DimPlot(Epithelium,label = TRUE,repel = TRUE)
DimPlot(Epithelium,label = TRUE,repel = TRUE,split.by = "orig.ident")



