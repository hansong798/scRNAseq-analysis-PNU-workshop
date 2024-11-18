# scRNAseq-analysis-PNU-workshop using R
#### Compiled by Hansong Lee
### Table of Content
  * [Preparation](#preparation)
  * [Analysis](#analysis)
    * [Step 1. Read files and create seurat objects](#step-1-read-files-and-create-seurat-objects)
    * [Step 2. Quality control](#step-2-quality-control)
    * [Step 3. Integration](#step-3-integration)
    * [Step 4. Run UMAP on a single integrated dataset](#step-4-run-umap-on-a-single-integrated-dataset)
    * [Step 5. Clustering](#step-5-clustering)
    * [Step 6. Annotation](#step-6-annotation)
    * [Step 7. Save the result](#step-7-save-the-result)



 ## Preparation
### Download dataset
Download GEO dataset from [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) with accession number GSE244515.

We will only download only 4 samples, 'healthy control 1', 'healthy control 2', 'PD1', 'PD2', for practice.
Please, Use 'custom' button.

Next, uncompress the .tar file. 
Please do not uncompress .gz files.

!! You should note that
* Create new folder for each sample and move files into each sample folder.
* File name: 'barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz'

### Download packages
In this workshop, we will use Seurat 4.4.0 and SeuratObject 4.1.4.
```R
packageVersion('Seurat')        # Caution with Seurat v5!
packageVersion('SeuratObject')

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_version("SeuratObject", "4.1.4", dependencies = T)   # ~15min
remotes::install_version("Seurat", "4.4.0", upgrade = FALSE)          # ~15min

install.packages('ggplot2')
install.packages("dplyr")
```


 ## Analysis
 ### Step 1. Read files and create seurat objects

```R
library(Seurat); library(ggplot2); library(dplyr)
setwd('D:/0_Workshop')   # <<-------- You need to write your directory
H1 <- Read10X(data.dir = 'H1')
H2 <- Read10X(data.dir = 'H2')
PD1 <- Read10X(data.dir = 'PD1')
PD2 <- Read10X(data.dir = 'PD2')
#H1

seurat_H1 <- CreateSeuratObject(H1,  project = "H1")
seurat_H2 <- CreateSeuratObject(H2,  project = "H2")
seurat_PD1 <- CreateSeuratObject(PD1,  project = "PD1")
seurat_PD2 <- CreateSeuratObject(PD2,  project = "PD2")
# dim(seurat_H1)
# seurat_H1
# seurat_H1@meta.data
seurat_H1$orig.ident; length(seurat_H1$orig.ident)
seurat_H1$nCount_RNA; length(seurat_H1$nCount_RNA)
seurat_H1$nFeature_RNA; length(seurat_H1$nFeature_RNA)

rawdata <- merge(x = seurat_H1, y = c(seurat_H2, seurat_PD1, seurat_PD2),
                      add.cell.ids = c('H1','H2','PD1','PD2'))
# dim(rawdata)
# table(rawdata$orig.ident)
head(rawdata@meta.data)
tail(rawdata@meta.data)

```


### Step 2. Quality control
```R
rawdata[["percent.mt"]] <- PercentageFeatureSet(rawdata, pattern = "^MT-") # "^mt-" in mouse
VlnPlot(rawdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(rawdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

VlnPlot(rawdata, features = c("nFeature_RNA"), pt.size = 0) + geom_hline(yintercept = 5000)
VlnPlot(rawdata, features = c("nCount_RNA"), pt.size = 0) + geom_hline(yintercept = 25000)
VlnPlot(rawdata, features = c("percent.mt"), pt.size = 0) + geom_hline(yintercept = 20)

rawdata@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "blue") +
  theme_classic()
  # geom_hline(yintercept = 600)
```

![image](https://github.com/user-attachments/assets/a1bf64f3-a592-4c12-bf22-27a0c77dafb0)
![image](https://github.com/user-attachments/assets/57583d6c-264a-4b72-869b-7655a2863a39)


We can check a optimal threshold using 'geom_hline' from ggplot2.


```R
filtdata <- subset(rawdata, subset = nFeature_RNA > 600 & nFeature_RNA < 5000 & nCount_RNA < 25000 & percent.mt < 20)
```


### Step 3. Integration

* Option 1 - Log normalization
```R
exdata <- NormalizeData(exdata)
exdata <- FindVariableFeatures(exdata)  # default, nfeatures = 2000
exdata <- ScaleData(exdata)
```

* Option 2 - SCTransform
```R
exdata <- SCTransform(exdata)
```

- Log Normalization: 'simple, fast' vs 'sequencing depth, zero count, asymmetric distribution'

- SCTransform: 'sequencing depth, zero count' vs 'problems with multiple data comparison (information from the other cells to avoid overestimation), random sampling'

We will use Log normalization in this workshop.

#### 3-1 Normalization
```R
splitdata <- SplitObject(filtdata, split.by = "orig.ident")
splitdata <- lapply(X = splitdata, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
```

#### 3-2 select repeatedly variable features across datasets for integration
```R
features <- SelectIntegrationFeatures(object.list = splitdata)
```

#### 3-3 Scale and run pca
```R
splitdata <- lapply(X = splitdata, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
```

#### 3-3 Integrate data
```R
anchors <- FindIntegrationAnchors(object.list = splitdata, anchor.features = features, reduction = "rpca")
intdata <- IntegrateData(anchorset = anchors)
DefaultAssay(intdata)
```


### Step 4. Run UMAP on a single integrated dataset
```R
DefaultAssay(intdata) <- "integrated"
intdata <- ScaleData(intdata, verbose = FALSE)
intdata <- RunPCA(intdata, npcs = 50, verbose = FALSE)
# choose npcs 
ElbowPlot(intdata, ndims = 50)
cumsd <- cumsum(intdata@reductions$pca@stdev)/sum(intdata@reductions$pca@stdev)
which(cumsd > 0.8)

intdata <- RunUMAP(intdata, reduction = "pca", dims = 1:33)
DimPlot(intdata, group.by = 'orig.ident')
```
![image](https://github.com/user-attachments/assets/e7cac226-7ea5-4961-a96e-0ee7bab67f31)


* Before integration

Compare the upper UMAP with below.

![image](https://github.com/user-attachments/assets/e5c776d8-4101-4fa5-b9c7-2df030dcd5c4)



### Step 5. Clustering
```R
intdata <- FindNeighbors(intdata, reduction = "pca", dims = 1:33)
intdata <- FindClusters(intdata, resolution = 0.2)
intdata$integrated_snn_res.0.2
DimPlot(intdata, group.by = 'integrated_snn_res.0.2', label = T) 
```
![image](https://github.com/user-attachments/assets/487c9207-f456-4d12-8520-1d999b12d1fc)




### Step 6. Annotation
```R
intdata@assays
DefaultAssay(intdata) <- 'RNA'
T_marker <- c("CD3E","CD3D","CD3G",'CD4','CD8A','CD8B')
NK_marker <- c('NKG7','PRF1','NCAM1') 
B_marker <- c('MS4A1','FCER2','CD19','CD22')
plasma_marker <- c('MZB1','IGHG1')
monocyte_marker <- c('LYZ','S100A8','S100A9')
pDC_marker <- c('LILRA4','IL3RA','CLEC4C')
mDC_marker <- c('CD1C','ITGAX','CLEC10A')
platelet_marker <- c('PPBP','PF4')
RBC_marker <- c('HBA1','HBA2','HBB')
proliferating_marker <- c('MKI67','PCNA','CDK1', 'CDK2')
progenitor_marker <- c('CDK6', 'CYTL1')

DotPlot(intdata, features = c(T_marker, NK_marker, B_marker,  plasma_marker, monocyte_marker, pDC_marker, mDC_marker,
                              platelet_marker, RBC_marker, proliferating_marker, progenitor_marker))

celltype <- c('0' = 'CD4T',
              '1' = 'CD4T',
              '2' = 'monocyte',
              '3' = 'CD8T',
              '4' = 'NK',
              '5' = 'CD8T',
              '6' = 'NK',
              '7' = 'CD8T',
              '8' = 'B',
              '9' = 'B',
              '10' = 'monocyte',
              '11' = 'NK',
              '12' = 'CD8T',
              '13' = 'mDC',
              '14' = 'CD8T',
              '15' = 'pDC',
              '16' = 'platelet',
              '17' = 'proliferating',
              '18' = 'progenitor'
)

Idents(intdata) <- 'integrated_snn_res.0.2'
intdata <- RenameIdents(intdata, celltype)
celltype <- Idents(intdata)
intdata$celltype <- celltype
DimPlot(intdata, label = T)

ggplot(data = intdata@meta.data, aes(x = orig.ident, fill = celltype)) +
  geom_bar(position = 'fill', colour="black")
```

![image](https://github.com/user-attachments/assets/407c68cd-1835-4802-9448-34991e87e5e3)

![image](https://github.com/user-attachments/assets/a9676148-f868-4286-a8ea-211e9d364fbb)


### Step 7. Save the result
```R
# Way 1
saveRDS(seurat, file="DS1/seurat_obj_all.rds")
saveRDS(seurat_dorsal, file="DS1/seurat_obj_dorsal.rds")
seurat <- readRDS("DS1/seurat_obj_all.rds")
seurat_dorsal <- readRDS("DS1/seurat_obj_dorsal.rds")

# Way 2
save(intdata, file = 'D:/0_Workshop/RData/intdata.RData')
```

