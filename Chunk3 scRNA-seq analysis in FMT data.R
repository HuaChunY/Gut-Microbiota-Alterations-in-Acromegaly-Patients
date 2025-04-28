####################################################################################################################################################################
#      ☆        % Project: Gut Microbiota Alterations in Acromegaly Patients Are Associated with Neutrophil Depletion-Induced Inflammation #
#   ☆ \|/ ☆    % Author: HuaChun Yin                                         
#  ☆  \|/  ☆   % Date: Apr. 4th, 2025                                  
# ☆   \|/   ☆  %                                                          
#  ☆  \|/  ☆   % Environment:   R version 4.4.2           
#  ☆ __|__ ☆   % EPlatform: Mac-IOS(64-bit)                                  
#                % CHUNK3:this script for scRNA-seq of colon from mouse with FMT                   
################################################################################################################################################################### 

packages <- c("dplyr","doParallel",
              "ggplot2","ggpubr",
              "psych","plotly",
              "stringr","corrplot","RColorBrewer",
              "ggsci","microeco",
              "tidyverse","magrittr",
              "data.table","foreach",
              "lme4","nlme",
              "factoextra","vegan",
              "ggalluvial","PMA",
              "ggrepel","scales","hdi",
              "stabs","Rtsne","ropls",
              "Seurat","patchwork",
              "Matrix","harmony","cowplot","ggh4x",
              "readr", "purrr","ReporterScore","org.Mm.eg.db",
              "clusterProfiler") ## package methods is not loaded by default by RScript. 
#check.packages(packages)
lapply(packages, library, character.only = TRUE)
## 
 
GH <- Read10X(data.dir = "outs/")
Ctrl <- Read10X(data.dir = "outs/")
 
## 
## Initialize the Seurat object with the raw (non-normalized data)
GH <- CreateSeuratObject(counts = GH , project = "GH", 
                         min.cells = 10, #gene  
                         min.features = 200)#cell  

Ctrl <- CreateSeuratObject(counts = Ctrl , project = "Ctrl", 
                            min.cells = 10, #gene 
                            min.features = 200)#cell  

########
mergedt <- Ctrl

 
## The operator can add columns to object metadata. This is a great place to stash QC stats
mergedt[["percent.mt"]] <- PercentageFeatureSet(mergedt, pattern = "^mt-")
mergedt[["percent.rp"]] <- PercentageFeatureSet(mergedt, pattern = "^Rp[sl]")# 
mergedt[["percent.hb"]] <- PercentageFeatureSet(mergedt, pattern = "^Hb[^(p)]")# 
 
# Visualize QC metrics as a violin plot
VlnPlot(mergedt, 
        features = c("nFeature_RNA",
                     "nCount_RNA", 
                     "percent.mt",
                     "percent.rp",
                     "percent.hb"),
        ncol = 3,pt.size = 0,# 
        group.by = "orig.ident")
###############



################
 
#GH
mergedt1 <- subset(mergedt, 
                   subset = nFeature_RNA > 500 & nFeature_RNA < 5500 & percent.mt < 10) # 
#Ctrl
mergedt1 <- subset(mergedt, 
                   subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 10) # 
 
VlnPlot(mergedt1, 
        features = c("nFeature_RNA",
                     "nCount_RNA", 
                     "percent.mt",
                     "percent.rp"),
        ncol = 3,pt.size = 0,# 
        group.by = "orig.ident")

dim(mergedt1)

mergedt1 <- NormalizeData(mergedt1, 
                          normalization.method = "LogNormalize")

mergedt1 <- FindVariableFeatures(mergedt1, selection.method = "vst", nfeatures = 2000)
mergedt1 <- ScaleData(mergedt1)
mergedt1 <- RunPCA(mergedt1)
mergedt1 <- RunUMAP(mergedt1, dims = 1:30)


##########################
library(DoubletFinder)
sweep.res.list <- paramSweep(mergedt1, PCs = 1:30)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
sweep.stats

bcmvn <- find.pK(sweep.stats)
mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))

nExp_poi <- round(0.02*ncol(mergedt1))        ## 3000~2.3%
nExp_poi

mergedt1 <- doubletFinder(mergedt1, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
mergedt1

a <- paste("DF.classifications_", "0.25_", mpK, "_", nExp_poi, sep="")

DimPlot(mergedt1, pt.size = 1, label = TRUE, label.size = 5, reduction = "umap", group.by = a)
a
 
## GH: DF.classifications_0.25_0.01_111
#Ctrl:DF.classifications_0.25_0.005_118
mergedt1 <- subset(mergedt1, DF.classifications_0.25_0.005_118 == "Singlet" )# 
 
after.data <- cbind(as.matrix(GetAssayData(object = mergedt1, assay = "RNA",slot = "counts")))

#GHafter <- CreateSeuratObject(counts = after.data, project = "GH")
Ctrlafter <- CreateSeuratObject(counts = after.data, project = "Ctrl")
 
####### 
mergedt <- merge(GHafter, y = c(Ctrlafter), 
                 add.cell.ids = c("GHPA", "Ctrl"))

mergedt[["RNA"]] <- JoinLayers(mergedt[["RNA"]])

dim(mergedt)
table(mergedt$orig.ident)
head(mergedt@meta.data)
 
## The operator can add columns to object metadata. This is a great place to stash QC stats
mergedt[["percent.mt"]] <- PercentageFeatureSet(mergedt, pattern = "^mt-")
mergedt[["percent.rp"]] <- PercentageFeatureSet(mergedt, pattern = "^Rp[sl]")# 
mergedt[["percent.hb"]] <- PercentageFeatureSet(mergedt, pattern = "^Hb[^(p)]")# 
 
# Visualize QC metrics as a violin plot
VlnPlot(mergedt3, 
        features = c("nFeature_RNA",
                     "nCount_RNA", 
                     "percent.mt",
                     "percent.rp",
                     "percent.hb"),
        ncol = 3,pt.size = 0,# 
        group.by = "orig.ident")
###############



################


mergedt <- NormalizeData(mergedt, 
                          normalization.method = "LogNormalize")

mergedt <- FindVariableFeatures(mergedt, 
                                 selection.method = "vst", 
                                nfeatures = 2000)

mergedt <- mergedt %>%
  ScaleData(features = rownames(.)) 

mergedt <- RunPCA(mergedt,npcs = 30, verbose = TRUE) 
 
#########################
mergedt2 <- mergedt %>%
  RunHarmony("orig.ident") %>% 
  FindNeighbors(dims = 1:30,reduction = "harmony") 

mergedt2 <- mergedt2 %>%
  FindClusters(resolution = 0.3) #0.5
  
  mergedt2 <- mergedt2 %>%
  RunUMAP(dims = 1:12,reduction = "harmony")

#DimPlot(mergedt2,label = TRUE, reduction = "umap")
#DimPlot(mergedt2,label = TRUE, reduction = "umap",group.by = "orig.ident")
##################
  col.sc <- c("#9BCFE6",
              "#84CC8C",
              "#0082BD",
              "#29B534",
              "#AE9B66",
              "#FFA631",
              "#FF8200",
              "#C5A5D2",
              "#713A9F",
              "#17becf", 
              "#F2DCAF","#7f7f7f", "#bcbd22", "#F2AFD0")


####################cell annotation
#immune cell: 10.1038/s41421-023-00578-4
#gut cell: 10.1016/j.cell.2022.09.024
#TA from ZZL file10.1016/j.cell.2022.09.024

DotPlot(mergedt2, features = c(
  "Mep1b","Cyp3a13","Cyp4f14",
  "Hsd17b2","Slc26a3","Apol10a",
  #Enterocytes
  "Cd3d",'Cd3e',"Trac",#T.cell
  "Htr4","Kcnq1","Cdo1",#Stem cell
  'Col1a1','Dcn',"Igfbp7","Sparc",#Fibroblasts
  "Agr2","Spink4","Muc2","Ramp1",#Goblet cells
  
  # 'Cd79a','Mzb1',#B.cell.marker
  "Ms4a1","Cd19",'Cd79a', #B
  
  'Mzb1',"Igkc", "Jchain",#Plasma Cells.marker
  
  "Ms4a6c","Lyz2",'Cd68',# Macrophages cell
  
  'Cldn5','Flt1','Ramp2',# Endothelial.cell.marker
  "Il1b","S100a8",#Neutrophils
  "Stmn1","Tubb5",#TA
  "Chgb","Chga","Cpe",#Enteroendocrine
  "Pnliprp2","Pnliprp1","Clps",#Paneth cells 
  "Tppp3","Rhoc","Arl5a"#Tuft cells
)) +
  RotatedAxis(90)+
  theme(text=element_text(family="serif",size=10))

ident <- c("Enterocytes", #0
           "T cell",#1
           "Stem cell",#2
           "Fibroblasts",#3
           "Fibroblasts",#4
           "B cell",#5
           "Goblet cells",#6
           "Macrophages cell",#7
           "B cell",#8
           "Macrophages cell",#9
           "Fibroblasts",#10
           "Endothelial cell",#11
           "Fibroblasts",#12
           "Neutrophils",#13
           "TA",#14
           "Endothelial cell",#15
           "Enteroendocrine",#16
           "Paneth cells",#17
           "Tuft cells"#18
           )
mergedt3 <- mergedt2
names(ident) <- levels(mergedt3)

mergedt3 <- RenameIdents(mergedt3,ident)
###########################
mergedt3$"cell_type" <- Idents(mergedt3)


subset <- subset(mergedt3, idents = "TA")

subset <- subset %>%
  SCTransform(verbose = TRUE) %>% 
  RunPCA(npcs = 15, verbose = TRUE) %>% 
  FindNeighbors(dims = 1:15,reduction = "harmony") %>% 
  FindClusters(resolution = 0.3) %>% 
  RunUMAP(dims = 1:10,reduction = "harmony")

DimPlot(subset, label = TRUE, reduction = "umap")

DotPlot(subset, features = c( "Stmn1","Tubb5",#TA
                              'Cd79a','Mzb1','Cd74'#B.cell.marker
)) +
  RotatedAxis(90)+
  theme(text=element_text(family="serif",size=12))

subset_ident <- c("TA", "B cell","B cell")
names(subset_ident) <- levels(subset)
subset <- RenameIdents(subset, subset_ident)
subset_names <- Idents(subset)


mergedt3$cell_type <- as.character(mergedt3$cell_type)
index <- match(names(subset_names), names(mergedt3$cell_type))
value <- unname(subset_names)
value <- as.character(value)
mergedt3$cell_type[index] <- value
Idents(mergedt3) <- mergedt3$cell_type
#########
 
subset <- subset(mergedt3, idents = "B cell")

subset <- NormalizeData(subset, 
                         normalization.method = "LogNormalize")

subset <- FindVariableFeatures(subset, 
                                selection.method = "vst", 
                                nfeatures = 2000)

subset <- subset %>%
  ScaleData(features = rownames(.)) 

subset <- subset(subset,npcs = 20, verbose = TRUE) 
 
subset <- subset %>%
  RunHarmony("orig.ident") %>% 
  FindNeighbors(dims = 1:20,reduction = "harmony") 

subset <- subset %>%
  FindClusters(resolution = 0.3) #0.5

subset <- subset %>%
  RunUMAP(dims = 1:12,reduction = "harmony")


DimPlot(subset, label = TRUE, reduction = "umap")

DotPlot(subset, features = c(   "Ms4a1","Cd19",'Cd79a', #B
                                
                                'Mzb1',"Igkc", "Jchain"#Plasma Cells.marker
)) +
  RotatedAxis(90)+
  theme(text=element_text(family="serif",size=12))

subset_ident <- c("B cell","Plasma cell","Plasma cell","Plasma cell","Plasma cell" )
names(subset_ident) <- levels(subset)
subset <- RenameIdents(subset, subset_ident)
subset_names <- Idents(subset)


mergedt3$cell_type <- as.character(mergedt3$cell_type)
index <- match(names(subset_names), names(mergedt3$cell_type))
value <- unname(subset_names)
value <- as.character(value)
mergedt3$cell_type[index] <- value
Idents(mergedt3) <- mergedt3$cell_type

#Idents(mergedt3) <- gsub("_GH", "", Idents(mergedt3), fixed = TRUE)
#Idents(mergedt3) <- gsub("_Ctrl", "", Idents(mergedt3), fixed = TRUE)
####################################


#########
mergedt3$cell_type<- factor(mergedt3$cell_type,levels=c("Enterocytes", #0
                                                        "Stem cell",#3
                                                        "Fibroblasts",#4
                                                        "Goblet cells",#5
                                                        "Endothelial cell",#11
                                                        "Paneth cells",#14
                                                        "Enteroendocrine",#16
                                                        "TA",#17
                                                        "Tuft cells",#18
                                                        "T cell",#1
                                                        "Plasma cell",#2
                                                        "B cell",#8
                                                        "Macrophages cell",
                                                        "Neutrophils"
                                                        
))

annosc::plot_dim(mergedt3, fill=NA, col = col.sc[1:14],
                 show.ct = F,group.by = "cell_type",
                 #show.cls=c("Capillary"),
                 label.size = 3,
                 reduction="umap")

annosc::plot_dim(mergedt3, fill=NA, col = c("#FFA631", "#0082BD"),
                 show.ct = F,group.by = "orig.ident",
                 #show.cls=c("Capillary"),
                 label.size = 3,
                 reduction="umap")


####################################

#######

###########03 cell number


###################

table(mergedt3$orig.ident)# 

table(mergedt3$cell_type, mergedt3$orig.ident)# 
Cellratio <- prop.table(table(mergedt3$cell_type, mergedt3$orig.ident), margin = 2)# 
Cellratio <- as.data.frame(Cellratio)

Cellratio$num <- c(table(mergedt3$cell_type, mergedt3$orig.ident)[,1],table(mergedt3$cell_type, mergedt3$orig.ident)[,2])

colnames(Cellratio) <- c("cell_class","group","percentage","n")


library(ggplot2)

ggplot(Cellratio, aes(x = group, y = percentage*100,
                      fill = cell_class,
                      stratum = cell_class, alluvium = cell_class))+
  geom_stratum(width = 0.6, color='white')+
  geom_alluvium(alpha = 0.2,
                width = 0.6,
                color='white',
                linewidth = 1,
                curve_type = "linear")+
 
  scale_y_continuous(expand = c(0,0))+
  labs(y="Percentage (%)")+
  theme_bw()+
  theme(text = element_text(family="serif",size = 10,color = "black"),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        legend.title =  element_blank())+
  guides(fill=guide_legend(keywidth = 1.2, keyheight = 1.2)) +
  scale_fill_manual(values = col.sc[1:14])

#########
  ggplot(Cellratio, aes(x = group, y = cell_class, size = percentage, color = cell_class)) + 
  geom_point(alpha = 1) +
  scale_color_manual(values = col.sc[1:14]) +  #  
  theme_classic() +
  theme(text = element_text(family="serif",size = 10,color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right",
        legend.title = element_text(size = 10),  
        legend.text = element_text(size = 8) 
  ) +
  labs(y = "Cell proportion", x = NULL, size = "Cell proportion", color = "Clusters")
#########

ggplot(Cellratio, aes(x = group, y = n, fill = cell_class)) +
  geom_bar(stat = "identity") +  
  scale_fill_manual(values = col.sc[1:14]) +
  theme_classic()+
  theme(text = element_text(family="serif",size = 10,color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none"
  )+
  labs(y = "Number of cells", x = NULL)

#########
Cellratio$cellgroup <- paste0(Cellratio$cell_class,Cellratio$group)
Cellratio$cellgroup <- factor(Cellratio$cellgroup,
                              levels=c("EnterocytesCtrl","EnterocytesGH",
                                       "Stem cellCtrl","Stem cellGH",
                                       "FibroblastsCtrl","FibroblastsGH",
                                       "Goblet cellsCtrl","Goblet cellsGH" ,
                                       "Endothelial cellCtrl","Endothelial cellGH",
                                       "Paneth cellsCtrl","Paneth cellsGH",
                                       "EnteroendocrineCtrl","EnteroendocrineGH" ,  
                                       "TACtrl","TAGH",           
                                       "Tuft cellsCtrl","Tuft cellsGH",       
                                       "T cellCtrl", "T cellGH",
                                       "Plasma cellCtrl", "Plasma cellGH", 
                                       "B cellCtrl", "B cellGH",         
                                       "Macrophages cellCtrl","Macrophages cellGH",  
                                       "NeutrophilsCtrl","NeutrophilsGH" )
                              
)
ggplot(Cellratio, aes(x = cellgroup, y = n, fill = cellgroup)) +
  geom_bar(stat = "identity",width =0.7) +  # stat="identity"  
  scale_fill_manual(values = col.sc[rep(1:15,each=2)]) +
  
  coord_flip() +
  theme_classic()+
  theme(text = element_text(family="serif",size = 10,color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")+
  labs(y = "Number of cells", x = NULL)

##############################


####### 02 Marker searching


##############################
CELLDEG.list <- list()
for(i in 1:length(unique(mergedt3$cell_type))){
  CELLDEG <- FindMarkers(mergedt3, 
                         ident.1 = unique(mergedt3$cell_type)[i],  
                         verbose = FALSE, 
                         test.use = "wilcox",
                         only.pos = TRUE, 
                         min.pct = 0.3, 
                         logfc.threshold = 1)
  CELLDEG$Cluster <- rep(unique(mergedt3$cell_type)[i],nrow(CELLDEG))
  CELLDEG$Gene <- rownames(CELLDEG)
  CELLDEG.list[[i]] <- CELLDEG
}
CELLDEG_df <- do.call(rbind,CELLDEG.list)# 
#openxlsx::write.xlsx(CELLDEG_df,"Cell marker.xlsx")


##############################
gene.marker <- c(
  "Mep1b","Cyp3a13", "Slc26a3","Apol10a",
  #Enterocytes
  "Htr4","Kcnq1","Cdo1",#Stem cell
  
  'Col1a1','Dcn',#Fibroblasts
  
  "Agr2","Muc2",#Goblet cells
  
  'Cldn5','Flt1','Ramp2',# Endothelial.cell.marker
  "Pnliprp2","Pnliprp1","Clps",#Paneth cells 
  "Chgb","Chga","Cpe",#Enteroendocrine
  "Stmn1","Tubb5",#TA
  "Tppp3","Rhoc","Arl5a",#Tuft cells
  "Cd3d",'Cd3e',"Trac",#T.cell
  'Mzb1',"Igkc", "Jchain",#Plasma Cells.marker
  "Ms4a1","Cd19",'Cd79a', #B
  
  "Ms4a6c","Cd68", # Macrophages cell
  "Il1b","S100a8"#Neutrophils
  
)

gene.cluster <- c(rep("Enterocytes",4),
                  rep("Stem cell",3),
                  rep("Fibroblasts",2),
                  rep("Goblet cells",2),
                  rep("Endothelial cell",3),
                  rep("Paneth cells",3), 
                  rep( "Enteroendocrine",3), 
                  rep("TA",2),
                  rep("Tuft cells",3),
                  rep("T cell",3),
                  rep("Plasma cell",3),
                  rep("B cell",3),
                  rep("Macrophages cell",2),
                  rep("Neutrophils",2)
                  
)


p <-DotPlot(mergedt3,group.by = "cell_type",
            features=split(gene.marker,
                           gene.cluster),
            cols=c("#ffffff","#448444"))+
  RotatedAxis()+#来自Seurat
  theme(text=element_text(family="serif",size=8),
        panel.border=element_rect(color="black"),
        panel.spacing=unit(1,"mm"),
        axis.title=element_blank(),
        axis.text.x = element_text(angle=90, size=8),
        axis.text.y=element_blank()
  )
p

##############################
 

##############################

#############DEG

##############################
 
##################

mergedt3@meta.data$celltype.group <- paste(mergedt3@meta.data$cell_type, mergedt3@meta.data$orig.ident, sep = "_")

Idents(mergedt3) <- "celltype.group"

cellfordeg<-unique(mergedt3@meta.data$cell_type)[-11]
 

CELLDEG.list <- list()

for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(mergedt3, 
                         ident.1 = paste0(cellfordeg[i],"_GH"), 
                         ident.2 = paste0(cellfordeg[i],"_Ctrl"), 
                         verbose = FALSE, 
                         test.use = 'wilcox',
                         min.pct = 0.2)
  CELLDEG$Cluster <- rep(cellfordeg[i],nrow(CELLDEG))
  CELLDEG$Gene <- rownames(CELLDEG)
  CELLDEG.list[[i]] <- CELLDEG
}

combined_df <- do.call(rbind,CELLDEG.list)# 
#combined_df <- openxlsx::read.xlsx("combined_df0318.xlsx")
################
combined_deg <- combined_df[abs(combined_df$avg_log2FC)>0.5&combined_df$pct.1>0.3&
                              combined_df$p_val<0.05,]
 
#################

#10.1038/s41467-023-36707-6

############ 
########Calculate the percentage of DEG
up_deg <- combined_df[combined_df$avg_log2FC>0.5& combined_df$pct.1>0.3&
                        combined_df$p_val<0.05,]
down_deg <- combined_df[abs(combined_df$avg_log2FC)>0.5& combined_df$pct.1>0.3&
                          combined_df$avg_log2FC<0&
                          combined_df$p_val<0.05,]

percentage.up <- table(up_deg$Cluster)/table(combined_df$Cluster)

percentage.down <- table(down_deg$Cluster)/table(combined_df$Cluster)

percentage.Genes <- 1-percentage.up-percentage.down
###############
percentage.up <- as.data.frame(percentage.up)
percentage.down <- as.data.frame(percentage.down)
percentage.Genes <- as.data.frame(percentage.Genes )

percentage.up$label <- rep("GH",nrow(percentage.up))
percentage.down$label <- rep("Ctrl",nrow(percentage.down))
percentage.Genes$label <- rep("Genes",nrow(percentage.Genes))

percentage.dt <- rbind(percentage.up,percentage.down,percentage.Genes)

View(percentage.dt)
 
library(ggalluvial)

ordered_categories <- c("Enterocytes", #0
                        "Stem cell",#3
                        "Fibroblasts",#4
                        "Goblet cells",#5
                        "Endothelial cell",
                        "Enteroendocrine",#16
                        "TA",#17
                        "Tuft cells",#18
                        "T cell",#1
                        "Plasma cell",#2
                        "B cell",#8
                        "Macrophages cell",
                        "Neutrophils"
                        
)

 ggplot(percentage.dt, aes(x = Var1, y = Freq*100,
               fill = label,
               stratum = label, alluvium = label))+ 
  geom_stratum(width = 0.6, color='white')+
  geom_alluvium(alpha = 0.2,
                width = 0.6,
                color='white',
                linewidth = 1,
                curve_type = "linear")+
  geom_text(aes(label=Freq*100),  
            position = position_stack(vjust =0.5), #  
            color="white", family="serif",
            size=3)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(limits = ordered_categories) +
  labs(y="Percentage (%)")+
  theme_bw()+
  theme(text = element_text(family="serif",size = 9),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title = element_text(size=12))+
  guides(fill=guide_legend(keywidth = 1.2, keyheight = 1.2)) +
  scale_fill_manual(values = c("#000000","#F6B75D" ,"#99C7E1"))
 
Colors2 <- col.sc[c(1:10,12:14)]
 
combined_df <- combined_df %>% 
  mutate(
    label=case_when(
      combined_df$avg_log2FC>0.5& combined_df$pct.1>0.3&
        combined_df$p_val<0.05 ~ "DEGs in GH",
      combined_df$avg_log2FC< -0.5& combined_df$pct.1>0.3&
        combined_df$p_val<0.05 ~ "DEGs in Ctrl",
      TRUE ~ "Genes")
  )

 
up10 <-  combined_deg %>%
  group_by(Cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) 
down10 <-  combined_deg %>%
  group_by(Cluster) %>%
  slice_min(order_by = avg_log2FC, n = 5) 
up10$label <- rep("DEGs in GH",nrow(up10))
down10$label <- rep("DEGs in Ctrl",nrow(down10))

top10 <- rbind(up10,down10)
 
combined_df <- combined_df %>% 
  mutate(
    size=case_when(
      Gene%in%top10$Gene ~ 2,
      TRUE ~ 1)
  )
#提取非Top10的基因表格；
combined_df_2<- filter(combined_df,size==1)

###########
dfbar<-data.frame(x=ordered_categories,
                  y=c(15,15,15,8,14,12,12,10,12,13,8,8,9))
dfbar1<-data.frame(x=ordered_categories,
                   y=c(-5,-5,-9,-9,-5,-5,-5,-9,-9,-5,-9,-9,-6))
#绘制背景柱：
p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1
#######################
dfcol<-data.frame(x=ordered_categories,
                  y=0,
                  label=ordered_categories)
mycol<- Colors2

p <- p1+
  geom_jitter(data = top10,
              aes(x = Cluster, y = avg_log2FC, color = label),
              size = 0.35,
              width =0.4)+ 
  geom_jitter(data = combined_df_2,
              aes(x = Cluster, y = avg_log2FC, color = label),
              size = 0.25,
              width =0.4)+
  geom_tile(data = dfcol,
  aes(x=x,y=y),
  height=1,
  color = "black",
  fill = mycol,
  #alpha = 0.6,
  show.legend = F)+
  scale_color_manual(name=NULL,
                     values = c("#0082BD","#FFA631" ,"#000000"))+
  labs(x="Cluster",y="Average log2FC")+
  scale_x_discrete(limits = ordered_categories) +
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =6,
            family="serif",
            color ="white")+
  theme(text = element_text(family="serif",size = 10),#"Times New Roman","serif"
        panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_blank()
  )+
  geom_text_repel(
    data=top10,
    aes(x = Cluster, y = avg_log2FC,label=Gene),
    force = 1.2,
    size=6,family="serif",
    max.overlaps = Inf,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last")
  )
p
####6*12
ggsave("DEG火山图0318.pdf", 
       plot = p, 
       device = "pdf", 
       width = 12, height = 6)



#############pathway function and KEGG list

pathway_gene <- custom_modulelist_from_org(
  org = "mmu",
  feature = c("ko", "gene", "compound")[2]
)
pathway_gene$id <- gsub("mmu", "map",pathway_gene$id)
KEGGpathway <- ReporterScore::load_KO_htable()

pathway_gene <- merge(pathway_gene,
                      KEGGpathway[,c("level1_name" ,"level2_name","level3_id","level3_name")],
                      by.x="id",by.y="level3_id")
pathway_gene <- pathway_gene[!duplicated(pathway_gene$KOs),]#KEGG all pathway and gene

##############################
 

########################
 
########Enrichment


#########
background.num<- dim(mergedt3)[[1]]
Gene <- combined_deg[which(combined_deg$Cluster%in%"Neutrophils"),] 


enrichkeggall <- enrichment.kegg(df = Gene$Gene, #gene list
                                 background.num=background.num, 
                                 
                                 class1=c("Cellular Processes",  "Organismal Systems"),
                                 class2=c("Immune system","Cell growth and death"),
                                 #class2="Amino acid metabolism",
                                 adj.p="BH",#c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                 #   "fdr", "none")
                                 ORA="hypergeometric", 
                                 #hypergeometric test, Fisher test
                                 min_exist_KO=0)



enrichkeggall$ratio <- enrichkeggall$exist_k/enrichkeggall$K_num
enrichkeggall <- enrichkeggall[enrichkeggall$adj.p<0.05,]

enrichkeggall$Description<- gsub(" - Mus musculus \\(house mouse\\)", "",enrichkeggall$Description)


###########GSEA

gmt <- pathway_gene[which(pathway_gene$level3_name%in%enrichkeggall$Description),]
gmt$Description<- gsub(" - Mus musculus \\(house mouse\\)", "",gmt$Description)
 
gmt1 <- gmt %>%
  separate_rows(KOs, sep = ",") %>%
  rename(KO = KOs) 

gmt1 <- gmt1[,c("level3_name","KO")]

colnames(gmt1) <- c("term","gene")

Gene <- Gene %>% 
  arrange(desc(avg_log2FC))

geneList = Gene$avg_log2FC  
names(geneList) <- Gene$Gene  

gsea_custom_result <-  clusterProfiler::GSEA(geneList, 
                                             TERM2GENE = gmt1, 
                                             minGSSize = 1, 
                                             pvalueCutoff = 0.99, 
                                             verbose = F)

gsea.out.df <- gsea_custom_result@result

################

pathway <- merge(enrichkeggall,gsea.out.df[,c("ID", "enrichmentScore","NES")],by.x ="Description",by.y="ID" )

pathway <- pathway%>%select(-"level3_name")

colnames(pathway)[7] <- c("count")

pathway$count <- as.numeric(pathway$count)
pathway$adj.p <- as.numeric(pathway$adj.p)
 

 ggplot(pathway,aes(NES,reorder(Description,NES)))+
  geom_point(aes(size=count,color=-log10(adj.p)))+
  scale_color_gradient(low = "#c7e9c0", high = "#006d2c") +
  labs(color=expression(-log10(adj.p)),
       size="Number",
       y=NULL,
       x="NES Score",
       fill="-log10(adj P)")+
  #geom_vline(xintercept = c(1,-1), linetype = "dashed", color = "red")+
  labs(title="Pathway enrichment in Neutrophils cells")+ 
  theme_bw()+
  theme(axis.text.x = element_text(size=10,color = "black"),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=10),
        text=element_text(family="serif",size=10))

############GO enrichment
 
 #human gene mapping
# stemDEG<- combined_deg[which(combined_deg$Cluster%in%"Stem cell"),]
 
 stemDEG<- openxlsx::read.xlsx("stem cell.xlsx")
 
 homologene <- homologene::homologene(stemDEG$Gene, inTax = 10090, outTax = 9606)
 
 #openxlsx::write.xlsx(homologene,"stem cell_DEG_homologene.xlsx")
 
 stemhomologene <- merge(stemDEG,homologene,by.x="Gene",by.y="10090")
 
 ######################
 a<- openxlsx::read.xlsx("GOBP stem cell adj0.05.xlsx")
 
###################
 
 for(i in 1:length(a$ID)){
   ID<- a$ID[i]
   ANCESTOR <- GO.db::GOBPANCESTOR[[ID]]
  
   
   if (length(ANCESTOR) != 0) {
     result <- AnnotationDbi::select(
       GO.db,
       keys = ANCESTOR,
       columns = c("GOID", "TERM"),
       keytype = "GOID"
     )
    
     a[which(a$ID%in%ID),"Category"]<- paste(result$TERM, collapse = "; ") 
    
   } else {
     a[which(a$ID%in%ID),"Category"] <-NA
   }
   
 }
 
 a$ratio <- a$Hit.Count.in.Query.List/a$Hit.Count.in.Genome

 GOimmune_system <- a[grep("\\bimmune system process\\b", a$Category), ]
 
 a1 <- a[-which(a$ID%in%GOimmune_system$ID),]
 
 GOmetabolic <- a1[grep("\\bmetabolic process\\b", a1$Category), ]
 
 a2 <- a1[-which(a1$ID%in%GOmetabolic$ID),]
 
 GOstimulus <- a2[grep("\\bresponse to stimulus\\b", a2$Category), ]
 
 a3 <- a2[-which(a2$ID%in%GOstimulus$ID),]
 
 GOdeath<- a3[grep("\\bcell death\\b", a3$Category), ]
 
 a4 <- a3[-which(a3$ID%in%GOdeath$ID),]
 
 GOlocalization <- a4[grep("\\blocalization\\b", a4$Category), ]
 a5 <- a4[-which(a4$ID%in%GOlocalization$ID),]
 
 GOdevelopment <- a5[c(grep("\\bdevelopment\\b", a5$Category),
                       grep("\\bdevelopment\\b", a5$Name)), ]
 
 a6 <- a5[-which(a5$ID%in%GOdevelopment$ID),]
 
 GOimmune_system$Class <- rep("immune system",nrow(GOimmune_system))
 
 GOmetabolic$Class <- rep("metabolic process",nrow(GOmetabolic))
 
 GOstimulus$Class <- rep("response to stimulus",nrow(GOstimulus))
 
 GOdeath$Class<- rep("cell death",nrow(GOdeath))
 
 GOlocalization$Class<- rep("localization",nrow(GOlocalization))
 
 GOdevelopment$Class<- rep("development",nrow(GOdevelopment))
 
 a6$Class<- rep("other",nrow(a6))

 a_class <- rbind(GOimmune_system,
       GOmetabolic,
       GOstimulus,
       GOdeath,
       GOlocalization,
       GOdevelopment,
       a6)
 
 
 a_class <- unique(a_class)
 
 table(a_class$Class)
 
 
 gmt1 <- openxlsx::read.xlsx("GOimmune.xlsx")
 
 #############################
 
 
 stemhomologene <- stemhomologene %>% 
   arrange(desc(avg_log2FC))
 
 
 geneList = stemhomologene$avg_log2FC # 
 
 names(geneList) <- stemhomologene$`9606_ID` # 
 
 
 gsea_custom_result <-  clusterProfiler::GSEA(geneList, 
                                              TERM2GENE = gmt1, 
                                              minGSSize = 0, 
                                              pvalueCutoff = 1, 
                                              verbose = F)
 
 gsea.out.df <- gsea_custom_result@result
 
 gsea.out.df <- unique( gsea.out.df)
 
 GOimmune_system <- merge(GOimmune_system,gsea.out.df[,c("ID","NES")],by.x="Name",by.y="ID")
 
 #############
 colnames(GOimmune_system)[c(1,7,9)] <- c("Description","adj.p","count")
 
 GOimmune_system <- GOimmune_system[1:10,]
 
 ggplot(GOimmune_system,aes(NES,reorder(Description,NES)))+
   geom_point(aes(size=count,color=-log10(adj.p)))+
   scale_color_gradient(low = "#c7e9c0", high = "#006d2c") +
   labs(color=expression(-log10(adj.p)),
        size="Number",
        y=NULL,
        x="NES Score",
        fill="-log10(adj P)")+
   #geom_vline(xintercept = c(1,-1), linetype = "dashed", color = "red")+
   labs(title="Immune system in stem cells")+ 
   theme_bw()+
   theme(axis.text.x = element_text(size=10,color = "black"),
         panel.grid = element_blank(),
         axis.text.y = element_text(size=10),
         text=element_text(family="serif",size=10))
 
###########################################
 
 
 P1 <- FeaturePlot(mergedt3,features = "Piezo2",order = TRUE, 
                   cols=c("lightgray", "#D53E4F"),
                   reduction = "umap",pt.size = 0.1)+
   scale_x_continuous("")+scale_y_continuous("")
 
 VlnPlot(mergedt3, features =c('Piezo2'), split.by = "orig.ident", group.by = "cell_type", 
         pt.size = 1, ncol = 1,cols =c("#FFA631","#9BCFE6"))+
   
   theme(text = element_text(family="serif",size = 10,color = "black"),
         axis.text.x.bottom = element_text(family="serif",size = 10,color = "black"),
         legend.position = "right",
         panel.border=element_rect(color="black"),
         panel.spacing=unit(1,"mm"),
         axis.text.y=element_text(size=7))
 
 
 features <-c("Piezo2") 
 expr <- as.matrix(mergedt3@assays$RNA$data[features,])
 
 posi <- match(rownames(expr),rownames(mergedt3@meta.data))
 
 #提取出来的ts是一个行名为细胞，列名为目标基因、condition信息、celltype注释信息的表格
 ts <- cbind(expr,mergedt3@meta.data[posi,c("orig.ident","cell_type")])
 
 P2 <- ggplot(ts,aes(x = cell_type,y=expr,fill=orig.ident))+
   geom_violin(scale="width",alpha=0.8,width=0.6,size=0.8)+ 
   geom_jitter(alpha=0.5,aes(color=orig.ident),
               position=position_jitterdodge(jitter.width = 0.1, 
                                             jitter.height = 0, 
                                             dodge.width = 0.4
               ))+
   
   scale_fill_manual(values=c("#A40545", "#4B65AF"))+
   scale_color_manual(values = c("#A40545", "#4B65AF"))+
   stat_compare_means(aes(group=orig.ident),
                      method="t.test",
                      label="p.signif",
                      label.y=by(ts[,"expr"],ts$cell_type, max) + 0.02, size=4.5)+
   xlab("")+
   ylab("Expression Level")+
   theme_bw()+
   theme(text = element_text(family="serif",size = 10,color = "black"),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         panel.border=element_rect(size=1.2),
         axis.text.x=element_text(angle=45,size=10,vjust=1,hjust =1,color="black"),
         axis.text.y=element_text(size =10)) 
 P1/P2
 
 
 