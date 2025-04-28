####################################################################################################################################################################
#      ☆        % Project: Gut Microbiota Alterations in Acromegaly Patients Are Associated with Neutrophil Depletion-Induced Inflammation #
#   ☆ \|/ ☆    % Author: HuaChun Yin                                         
#  ☆  \|/  ☆   % Date: Apr. 4th, 2025                                  
# ☆   \|/   ☆  %                                                          
#  ☆  \|/  ☆   % Environment:   R version 4.4.2           
#  ☆ __|__ ☆   % EPlatform: Mac-IOS(64-bit)                                  
#                % CHUNK6:this script for scRNA-seq of colon from mouse with Bilophila.w                      
################################################################################################################################################################### 

packages <- c("dplyr","doParallel","harmony",
              "ggplot2","ggpubr","DoubletFinder",
              "psych","plotly",
              "stringr","corrplot",
               "RColorBrewer",
              "ggsci","microeco",
              "tidyverse","magrittr",
              "data.table","foreach",
              "lme4","nlme",
              "factoextra","vegan",
              "ggalluvial","PMA",
              "ggrepel","scales","hdi",
              "stabs","Rtsne","ropls",
              "Seurat","patchwork","Matrix"
)  
#check.packages(packages)
lapply(packages, library, character.only = TRUE)


#########################################################

#########01 Data filtering and doublet detection 

#########################################################

##data
EX <- Read10X(data.dir = "EX/")
CR <- Read10X(data.dir = "CR/")
 
## Initialize the Seurat object with the raw (non-normalized data)
EX <- CreateSeuratObject(counts = EX , project = "EX", 
                         min.cells = 10, #gene  
                         min.features = 50)#cell  

CR <- CreateSeuratObject(counts = CR , project = "CR", 
                            min.cells = 10, #gene  
                            min.features = 50)#cell  

########
mergedt <- CR

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
        ncol = 3,pt.size = 0, 
        group.by = "orig.ident")
###############

mergedt1 <- subset(mergedt, 
                   subset = nFeature_RNA > 50 & nFeature_RNA < 7000 & percent.mt < 10) # 


VlnPlot(mergedt1, 
        features = c("nFeature_RNA",
                     "nCount_RNA", 
                     "percent.mt",
                     "percent.rp"),
        ncol = 3,pt.size = 0, 
        group.by = "orig.ident")

dim(mergedt1)

mergedt1 <- NormalizeData(mergedt1, 
                          normalization.method = "LogNormalize")

mergedt1 <- FindVariableFeatures(mergedt1, selection.method = "vst", nfeatures = 2000)
mergedt1 <- ScaleData(mergedt1)
mergedt1 <- RunPCA(mergedt1)
mergedt1 <- RunUMAP(mergedt1, dims = 1:30)


##########################

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
 
## 
#EX: DF.classifications_0.25_0.24_168
#CR: DF.classifications_0.25_0.23_138

mergedt1 <- subset(mergedt1, DF.classifications_0.25_0.23_138 == "Singlet" )# 

 
after.data <- cbind(as.matrix(GetAssayData(object = mergedt1, assay = "RNA",slot = "counts")))

#EX <- CreateSeuratObject(counts = after.data, project = "EX")
CR <- CreateSeuratObject(counts = after.data, project = "CR")

mergedt <- merge(EX, y = CR, 
                 add.cell.ids = c("EX", "CR"))

mergedt[["RNA"]] <- JoinLayers(mergedt[["RNA"]])
 
table(mergedt$orig.ident)
head(mergedt@meta.data)

########################################################################### 

####################02 Cell annotation

########################################################################### 
DotPlot(mergedt2, features = c(
  
  'Col1a1','Dcn',"Igfbp7","Sparc",#Fibroblasts
  
  "Ms4a6c","Cd68", "Csf1r","C1qb", "C1qa", "C1qc",# Macrophages cell
  
  "Ms4a1","Cd19",'Cd79a', #B
  
  "Cd3d",'Cd3e',"Trac",#T.cell
  
  'Cldn5','Flt1','Ramp2',# Endothelial 
  
  "Agr2","Spink4","Muc2","Ramp1",#Goblet cell
  
  'Mzb1',"Igkc", "Jchain",#Plasma cell
  
  "Htr4","Kcnq1","Cdo1",#Stem cell
  
  "Mep1b","Cyp3a13","Cyp4f14",
  "Hsd17b2","Slc26a3","Apol10a",
  #Enterocytes
  
  "Il1b","Itgam","Itgb2","Itgax","Alox5ap", 
  "Olfm1", "Ncf4", #Neutrophils Itgam(CD11b)
  
  "Stmn1","Tubb5",#TA
  
  "Chgb","Chga","Cpe",#Enteroendocrine
  
  "Pnliprp2","Pnliprp1","Clps",#Paneth cells 
  
  "Tppp3","Rhoc","Arl5a"#Tuft cells
)) +
  RotatedAxis(90)+
  theme(text=element_text(family="serif",size=12))


ident <- c( "Fibroblast", #0
            "Macrophage cell",#1
            "B cell",#2
            "Fibroblast", #3
            "T cell",#4
            "T cell",#5
            "Endothelial",#6
            "Goblet cell",#7
            "Fibroblast",#8
            "Endothelial",#9
            "Plasma cell", #10
            "Stem cell" ,#11
            "Enterocyte",#12
            "Goblet cell",#13
            "Fibroblast",#14
            "Stem cell" ,#15
            "Fibroblast",#16
            "Neutrophil cell",#17
            "Tuft cell",#18
            "un1",#19
            "Fibroblast", #20
            "TA",#21
            "Fibroblast", #22
            "Plasma cell", #23
            "un2",#24
            "un3",#25
            "Enteroendocrine" 
)


mergedt3 <- mergedt2

names(ident) <- levels(mergedt3)

mergedt3 <- RenameIdents(mergedt3,ident)

mergedt3$"cell_type" <- Idents(mergedt3)

Idents(mergedt3) <- mergedt3$cell_type

mergedt3 <- subset(mergedt3, idents = c( "Fibroblast",
                                         "Macrophage cell", 
                                         "B cell", 
                                         "T cell", 
                                         "Endothelial",
                                         "Goblet cell", 
                                         "Plasma cell", 
                                         "Stem cell" , 
                                         "Enterocyte", 
                                         "Neutrophil cell", 
                                         "Stem cell" , 
                                         "Tuft cell",
                                         "TA",
                                         "Enteroendocrine"))
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
                         logfc.threshold = 0.5)
  CELLDEG$Cluster <- rep(unique(mergedt3$cell_type)[i],nrow(CELLDEG))
  CELLDEG$Gene <- rownames(CELLDEG)
  CELLDEG.list[[i]] <- CELLDEG
}
CELLDEG_df <- do.call(rbind,CELLDEG.list) 
CELLDEG_df <- CELLDEG_df[CELLDEG_df$p_val<0.05,]

#openxlsx::write.xlsx(CELLDEG_df,"Cell marker.xlsx")

###########################
col.sc <- c("#9BCFE6", "#84CC8C", "#0082BD", "#29B534",  "#AE9B66", "#FF8200",
            "#C5A5D2",  "#713A9F","#17becf", "#F2DCAF","#7f7f7f", "#bcbd22", "#F2AFD0")

mergedt3$cell_type<- factor(mergedt3$cell_type,levels=c("Enterocyte", 
                                                        "Stem cell", 
                                                        "Fibroblast", 
                                                        "Goblet cell", 
                                                        "Endothelial", 
                                                        "Enteroendocrine", 
                                                        "TA", 
                                                        "Tuft cell", 
                                                        "T cell", 
                                                        "Plasma cell", 
                                                        "B cell", 
                                                        "Macrophage cell",
                                                        "Neutrophil cell"
                                                        
))

annosc::plot_dim(mergedt3, fill=NA, col = col.sc[1:13],
                 show.ct = F,group.by = "cell_type",
                 label.size = 3,
                 reduction="umap")


######################################

coords <- Embeddings(mergedt3, "umap")[, 1:2]

coords <- as.data.frame(coords)

colnames(coords) <- c("UMAP_1", "UMAP_2")

coords$sample <- mergedt3$orig.ident   
 
ggplot(coords, aes(x = UMAP_1, y = UMAP_2, color = sample)) +
  geom_point(size = 0.1, alpha = 0.6) +
  scale_color_manual(values = c("#FFA631","#9BCFE6")) +
  theme_classic() +
  facet_wrap(~ sample)+   
  guides(color = guide_legend(override.aes = list(size = 3)))  


############################################################################
 
###########03 Cell number

############################################################################
 
table(mergedt3$orig.ident) 

table(mergedt3$cell_type, mergedt3$orig.ident) 

Cellratio <- prop.table(table(mergedt3$cell_type, mergedt3$orig.ident), margin = 2)

Cellratio <- as.data.frame(Cellratio)

Cellratio$num <- c(table(mergedt3$cell_type, mergedt3$orig.ident)[,1],table(mergedt3$cell_type, mergedt3$orig.ident)[,2])

colnames(Cellratio) <- c("cell_class","group","percentage","n")
 
ggplot(Cellratio, aes(x = group, y = percentage*100,
                      fill = cell_class,
                      stratum = cell_class, alluvium = cell_class))+
  geom_stratum(width = 0.6, color='white')+
  geom_alluvium(alpha = 0.2,
                width = 0.6,
                color='white',
                linewidth = 1,
                curve_type = "linear")+
  geom_text(aes(label=round(percentage*100,2)),  
            position = position_stack(vjust =0.5),  
            color="black", family="serif",
            size=3)+
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
  scale_fill_manual(values = col.sc[1:13])

#########
ggplot(Cellratio, aes(x = group, y = cell_class, size = percentage, color = cell_class)) + 
  geom_point(alpha = 1) +
  scale_color_manual(values = col.sc[1:14]) +  
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
                              levels=c("EnterocyteCR","EnterocyteEX",
                                       "Stem cellCR","Stem cellEX",
                                       "FibroblastCR","FibroblastEX",
                                       "Goblet cellCR","Goblet cellEX" ,
                                       "EndothelialCR","EndothelialEX",
                                       "Paneth cellCR","Paneth cellEX",
                                       "EnteroendocrineCR","EnteroendocrineEX" ,  
                                       "TACR","TAEX",           
                                       "Tuft cellCR","Tuft cellEX",       
                                       "T cellCR", "T cellEX",
                                       "Plasma cellCR", "Plasma cellEX", 
                                       "B cellCR", "B cellEX",         
                                       "Macrophage cellCR","Macrophage cellEX",  
                                       "Neutrophil cellCR","Neutrophil cellEX")
                              
)
ggplot(Cellratio, aes(x = cellgroup, y = n, fill = cellgroup)) +
  geom_bar(stat = "identity",width =0.7) +  # stat="identity" 使用原始数据的值
  scale_fill_manual(values = col.sc[rep(1:15,each=2)]) +
  
  coord_flip() +
  theme_classic()+
  theme(text = element_text(family="serif",size = 10,color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")+
  labs(y = "Number of cells", x = NULL)

##############################

gene.marker <- c(
  "Mep1b","Cyp3a13", "Slc26a3","Apol10a",
  #Enterocytes
  "Htr4","Kcnq1","Cdo1",#Stem cell
  'Col1a1','Dcn',#Fibroblasts
  "Agr2","Muc2",#Goblet cells
  'Cldn5','Flt1','Ramp2',# Endothelial.cell.marker 
  "Chgb","Chga","Cpe",#Enteroendocrine
  "Stmn1","Tubb5",#TA
  "Tppp3","Rhoc","Arl5a",#Tuft cells
  "Cd3d",'Cd3e',"Trac",#T.cell
  'Mzb1',"Igkc", "Jchain",#Plasma Cells.marker
  "Ms4a1","Cd19",'Cd79a', #B
  "Cd68", # Macrophages cell
  "C1qb", "C1qa", "C1qc",
  "Il1b", "Itgb2","Itgax","Alox5ap", #Neutrophils Itgam(CD11b)
  "Olfm1" 
)

gene.cluster <- c(rep("Enterocyte",4),
                  rep("Stem cell",3),
                  rep("Fibroblast",2),
                  rep("Goblet cell",2),
                  rep("Endothelial cell",3),
                  rep( "Enteroendocrine",3), 
                  rep("TA",2),
                  rep("Tuft cells",3),
                  rep("T cell",3),
                  rep("Plasma cell",3),
                  rep("B cell",3),
                  rep("Macrophage cell",4),
                  rep("Neutrophil cell",5)
                  
)


DotPlot(mergedt3,group.by = "cell_type",
        features=split(gene.marker,
                       gene.cluster),
        cols=c("#ffffff","#448444"))+
  RotatedAxis()+ 
  theme(text=element_text(family="serif",size=8),
        panel.border=element_rect(color="black"),
        panel.spacing=unit(1,"mm"),
        axis.title=element_blank(),
        axis.text.x = element_text(angle=90, size=8),
        axis.text.y=element_blank()
  )

##############################


#############DEG


##################

mergedt3@meta.data$celltype.group <- paste(mergedt3@meta.data$cell_type, mergedt3@meta.data$orig.ident, sep = "_")

Idents(mergedt3) <- "celltype.group"

cellfordeg<-unique(mergedt3@meta.data$cell_type) 


CELLDEG.list <- list()

for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(mergedt3, 
                         ident.1 = paste0(cellfordeg[i],"_EX"), 
                         ident.2 = paste0(cellfordeg[i],"_CR"), 
                         verbose = FALSE, 
                         test.use = 'wilcox',
                         min.pct = 0.1)
  CELLDEG$Cluster <- rep(cellfordeg[i],nrow(CELLDEG))
  CELLDEG$Gene <- rownames(CELLDEG)
  CELLDEG.list[[i]] <- CELLDEG
}

combined_df <- do.call(rbind,CELLDEG.list) 
# openxlsx::write.xlsx(combined_df,"combined_df.xlsx")
combined_df <- openxlsx::read.xlsx("combined_df.xlsx")
################
combined_deg <- combined_df[abs(combined_df$avg_log2FC)>0.5&  
                              combined_df$p_val<0.05,]

# openxlsx::write.xlsx(combined_deg,"combined_deg.xlsx")

#################
combined_df2 <- combined_df[which(combined_df$Cluster%in%"Neutrophil cell"),]
 
markers <- combined_df2%>%
  mutate(Difference = pct.1-pct.2)
  
markers$sign <- ifelse(
  markers$avg_log2FC > 0.5 & markers$p_val < 0.05, "up",
  ifelse(
    markers$avg_log2FC < -0.5 & markers$p_val < 0.05, "down", "ns"
  )
)
 
ggplot(markers, aes(x=Difference, y=avg_log2FC, color = sign)) + 
  geom_point(size=1.2) + 
  scale_color_manual('sign',labels=c(paste0("down(",table(markers$sign)[[1]],')'),
                                     'ns',
                                     paste0("up(",table(markers$sign)[[3]],')' )),
                     values=c("#0082BD", "grey","#F5AF5B" ))+
  geom_label_repel(data=subset(markers, 
                               avg_log2FC >= 1 &p_val <= 0.05), 
                   aes(label=Gene),  
                   color="black",  
                   segment.colour = "black", 
                   label.padding = 0.2, 
                   segment.size = 0.3,  
                   size=4,
                   max.overlaps= 20) +   
  geom_label_repel(data=subset(markers, 
                               avg_log2FC <= -1 & p_val <= 0.05), 
                   aes(label=Gene), 
                   label.padding = 0.2, 
                   color="black",
                   segment.colour = "black",
                   segment.size = 0.3, size=4,
                   max.overlaps = 20) + 
  geom_vline(xintercept = 0,linetype = 2) +
  geom_hline(yintercept = 0,linetype = 2) +
  labs(x="Percentage difference",y="Avg_log2FC") + 
  theme_bw()+
  theme(text = element_text(family="serif",size = 8,color = "black"),#"Times New Roman","serif"
        panel.grid = element_blank()
  )

#################

#10.1038/s41467-023-36707-6

#############

########Pathway enrichment

#########
pathway_gene <- ReporterScore::custom_modulelist_from_org(
  org = "mmu",
  feature = c("ko", "gene", "compound")[2]
)
pathway_gene$id <- gsub("mmu", "map",pathway_gene$id)
KEGGpathway <- ReporterScore::load_KO_htable()

pathway_gene <- merge(pathway_gene,
                      KEGGpathway[,c("level1_name" ,"level2_name","level3_id","level3_name")],
                      by.x="id",by.y="level3_id")
pathway_gene <- pathway_gene[!duplicated(pathway_gene$KOs),]#KEGG all pathway and gene

background.num<- dim(mergedt3)[[1]]
#20368

#Gene <- combined_deg[which(combined_deg$Cluster%in%"Enteroendocrine"),] 

Gene <- combined_deg[which(combined_deg$Cluster%in%"Neutrophil cell"),] 


enrichkeggall <- enrichment.kegg(df = Gene$Gene, #gene list
                                 background.num=background.num, ,
                                 class1=c("Cellular Processes",  "Organismal Systems"),
                                 class2=c("Immune system","Cell growth and death"),
                                 adj.p="BH",
                                 ORA="hypergeometric",  
                                 min_exist_KO=0)

enrichkeggall$ratio <- enrichkeggall$exist_k/enrichkeggall$K_num

enrichkeggall <- enrichkeggall[enrichkeggall$adj.p<0.05,]

enrichkeggall$Description<- gsub(" - Mus musculus \\(house mouse\\)", "",enrichkeggall$Description)
 

#####################
 
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
#openxlsx::write.xlsx(pathway,"GSEA&ORA_Neutrophils.xlsx")

ggplot(pathway,aes(NES,reorder(Description,NES)))+
  geom_point(aes(size=count,color=-log10(adj.p)))+
  scale_color_gradient(low = "#c7e9c0", high = "#006d2c") +
  labs(color=expression(-log10(adj.p)),
       size="Number",
       y=NULL,
       x="NES Score",
       fill="-log10(adj P)")+
  labs(title="Neutrophil cell")+ 
  theme_bw()+
  theme(axis.text.x = element_text(size=10,color = "black"),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=10),
        text=element_text(family="serif",size=10))



