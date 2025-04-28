####################################################################################################################################################################
#      ☆        % Project: Gut Microbiota Alterations in Acromegaly Patients Are Associated with Neutrophil Depletion-Induced Inflammation #
#   ☆ \|/ ☆    % Author: HuaChun Yin                                         
#  ☆  \|/  ☆   % Date: Apr. 4th, 2025                                  
# ☆   \|/   ☆  %                                                          
#  ☆  \|/  ☆   % Environment:   R version 4.4.2           
#  ☆ __|__ ☆   % Platform: Mac-IOS(64-bit)                                  
#                % CHUNK5:this script includes metaproteomics  mapping, missing value imputation and normalization             
################################################################################################################################################################### 
 
packages <- c("dplyr","doParallel","multiUS","ReporterScore","vsn", 
              "ggplot2","ggpubr","tidyr","ggfun",
              "psych","plotly",
              "stringr","corrplot",
               "RColorBrewer",
              "ggsci","microeco",
              "tidyverse","magrittr",
              "data.table","foreach",
              "lme4","nlme",
              "factoextra","vegan",
              "ggalluvial","PMA",
              "ggrepel","scales","hdi","ropls",
              "stabs","Rtsne","ropls","MicrobiotaProcess")  
#check.packages(packages)

lapply(packages, library, character.only = TRUE)

#########################################################

#########01 Data filtering

#########################################################

Na.num <- 5 ## 70% were detected at least  

protein<- openxlsx::read.xlsx("宏蛋白/protein.xlsx")

KEGGIDtransall<- openxlsx::read.xlsx("宏蛋白/KEGGIDtrans.xlsx")

KEGGIDtrans <- KEGGIDtransall[,1:2]

protein1 <- protein[,c(2,7:26)]

rownames(protein1) <- protein1$PG.ProteinAccessions

protein1 <- protein1[,2:21]

#filtering
keeprow <- c() 

removerow <- c()

for (i in 1:nrow(protein1)){
  
  if( 20-sum(nchar(gsub("NaN", "", protein1[i,]))==0) >Na.num){
    
    keeprow <- c(keeprow,i)
  }else{
    removerow <- c(removerow,i)
  }
  
}

clean_df <- protein1[keeprow, ]

dim(clean_df)

reove_df <- protein1[removerow, ]

dim(reove_df)

rownames <- rownames(clean_df)

####################################################################

######02 KNN miss value imputation

################################################################

clean_df[] <- lapply(clean_df, as.numeric)

e.na <-grep("A",colnames(clean_df))

c.na <- grep("B",colnames(clean_df))

clean_df1 <- multiUS::KNNimp(clean_df, k = 5) #k from DOI:10.1038/s41467-020-17916-9

clean_df2 <- as.data.frame(clean_df1)

clean_df2$ID <- rownames(clean_df2)

#####################

for (i in seq_len(nrow(clean_df2))) {
  
  row_name <- clean_df2$ID[i]
  
  if (grepl(";", row_name)) {  
    
    split_names <- strsplit(row_name, ";")[[1]] 
    
    new_df <- matrix(0, ncol = ncol(clean_df2), nrow = length(split_names))  
    
    new_df <- as.data.frame(new_df)
    
    for (ith in 1:length(split_names)) {
      
      new_df[ith,] <-clean_df2[i,]
      
    }
    colnames(new_df) <- colnames(clean_df2)
    
    new_df[,"ID"] <- split_names
 
    df.cb <- rbind(clean_df2, new_df)
    
  }else{
    df.cb<-clean_df2 
  }
  clean_df2 <- df.cb
}    

rows_with_semicolon <- grep(";", clean_df2$ID)

nonseicolondf <- clean_df2[-rows_with_semicolon,]# 

 

####################################################################

######03 Data preparation

################################################################

#######################001 KEGG ID mapping

KEGGdf <- merge(nonseicolondf,KEGGIDtrans,by.x= "ID",by.y = "Accession")

e.p <-grep("A",colnames(KEGGdf))

c.p <- grep("B",colnames(KEGGdf))

KEGGdf[,c(e.p,c.p)] <- lapply(KEGGdf[,c(e.p,c.p)], as.numeric)


#######################002 Normalization
 
par(mfrow = c(1, 2))

boxplot(log2(KEGGdf[,c(e.p,c.p)]),main = "Before normalization",  outline =F,names=F)

fit<- vsn2(as.matrix(KEGGdf[,c(e.p,c.p)]))
 
ynorm = predict(fit, as.matrix(KEGGdf[,c(e.p,c.p)]))

KEGGdfnorm = justvsn(ynorm)

boxplot(KEGGdfnorm,main = "After normalization", outline =F,names=F)

par(mfrow = c(1, 1))  #  

pmad <- mean(mad(KEGGdfnorm))#less than 0.14 

print(pmad)
 
data_long <- gather(as.data.frame(KEGGdfnorm), key = "Variable", value = "Value") 

data_long$Variable <- rep("A",nrow(data_long))

ggplot(data_long, aes(x = Variable, y = Value)) +
  geom_jitter(width=0.3,size=0.5, 
              color = "#5E7BBB", 
              alpha = 0.6)+
  labs(title = paste("PMAD:",round(pmad,3))) +
  theme(text = element_text(family="serif",size = 8),#"Times New Roman","serif"
        panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6),
        legend.position = "none"
  )

#DOI: 10.1093/bib/bbz061
#The RLA plots before and after the VSN normalization and the PMAD for data sets. 
#(A) The RLA plots for unnormalized intensities; (
#(B) the RLA plots for normalized intensities; 
#(C) the distributions of PMAD values. 

#######################003 Statistical significance of the difference between groups

KEGGdf1 <- KEGGdf

KEGGdf1[,c(e.p,c.p)] <- KEGGdfnorm

KEGGdf1$"Pvalue"<- apply(KEGGdf1[,c(e.p,c.p)], 1, function(row){
  
  wilcox.test(row[1:10], row[11:20])$p.value
  
})

#Statistical significance of the difference between groups was evaluated using Mann–Whitney U test (nonparametric test of the null hypothesis which is suitable for data that does not pass normality test), unless otherwise indicated.

#DOI:10.1038/s41467-020-17916-9

View(KEGGdf1)

KEGGdfduplicated_rows <- KEGGdf1[duplicated(KEGGdf1$KEGG) | duplicated(KEGGdf1$KEGG, fromLast = TRUE), ]

KEGGdfduplicated_p <- KEGGdfduplicated_rows[KEGGdfduplicated_rows$Pvalue<0.05,]

KEGGdfduplicated_pdu <- KEGGdfduplicated_p[duplicated(KEGGdfduplicated_p$KEGG) | duplicated(KEGGdfduplicated_p$KEGG, fromLast = TRUE), ]

length(unique(KEGGdfduplicated_p$KEGG))
#####################

#######################004 Fold change

rawFC <- clean_df

rawFC$Num_GHPA <- rep(NA,nrow(rawFC))

rawFC$Num_Ctrl <- rep(NA,nrow(rawFC))

rawFC$Agv_GHPA <- rep(NA,nrow(rawFC))

rawFC$Agv_Ctrl <- rep(NA,nrow(rawFC))

for (i in 1:nrow(rawFC)){
  
  epsun <- 10-sum(nchar(gsub("NaN", "", clean_df[i,e.na]))==0) 
  
  cpsun <- 10-sum(nchar(gsub("NaN", "", clean_df[i,c.na]))==0) 
  
  allsum <- 20-sum(nchar(gsub("NaN", "", clean_df[i,]))==0) 
  
  if (cpsun<3){
    ep.value <- as.numeric(gsub("NaN", "0", clean_df[i,e.na]))
    
    epagv <- mean(sum(ep.value[ep.value>0]))
    
    cpagv <- 0
  }
  if (epsun<3){
    cp.value <- as.numeric(gsub("NaN", "0", clean_df[i,c.na]))
    
    cpagv <- mean(sum(cp.value[cp.value>0]))
    
    epagv <- 0
  } 
  if(cpsun>2&epsun>2){
    ep.value <- as.numeric(gsub("NaN", "0", clean_df[i,e.na]))
    
    epagv <- mean(sum(ep.value[ep.value>0]))
    
    cp.value <- as.numeric(gsub("NaN", "0", clean_df[i,c.na]))
    
    cpagv <- mean(sum(cp.value[cp.value>0]))
  }
  
  rawFC$Num_GHPA[i] <- epsun
  
  rawFC$Num_Ctrl[i] <- cpsun
  
  rawFC$Agv_GHPA[i] <- epagv
  
  rawFC$Agv_Ctrl[i] <-cpagv
}

View(rawFC)

rawFC$FC<- rep(NA,nrow(rawFC))

for (ig in 1:nrow(rawFC)){
  Agv_GHPA <- rawFC[ig,"Agv_GHPA"]
  
  Agv_Ctrl <- rawFC[ig,"Agv_Ctrl"]
  
  if(Agv_GHPA==0){
    
    rawFC$FC[ig] <- 0.0000001
  }
  if(Agv_Ctrl==0){
    
    rawFC$FC[ig] <- 16
    
  }
  if(Agv_Ctrl!=0&Agv_GHPA!=0){
    
    rawFC$FC[ig] <-rawFC$Agv_GHPA[ig]/rawFC$Agv_Ctrl[ig]
  }
}

rawFC$log2FC<- log2(rawFC$FC)

rawFC1 <- rawFC

rawFC1$ID <- rownames(rawFC1)

for (i in seq_len(nrow(rawFC1))) {
  
  row_name <- rawFC1$ID[i]
  
  if (grepl(";", row_name)) {   
    
    split_names <- strsplit(row_name, ";")[[1]]  
    
    new_df <- matrix(0, ncol = ncol(rawFC1), nrow = length(split_names))  
    
    new_df <- as.data.frame(new_df)
    
    for (ith in 1:length(split_names)) {
      new_df[ith,] <- rawFC1[i,]
    }
    colnames(new_df) <- colnames(rawFC1)
    
    new_df[,"ID"] <- split_names
     
    df.cb <- rbind(rawFC1, new_df)
    
  }else{
    df.cb<-rawFC1
  }
  rawFC1 <- df.cb
}    

View(rawFC1)

rawFC1_semicolon <- grep(";", rawFC1$ID)

nonseicolonrawFC1 <- rawFC1[-rawFC1_semicolon,]

KEGGFC <- merge(nonseicolonrawFC1,KEGGIDtrans,by.x= "ID",by.y = "Accession")

#######################################

reslutdu <- merge(nonseicolonrawFC1,KEGGdf1,by= "ID")
 
reslutdu <- merge(reslutdu,KEGGIDtransall[,c(1,3)],by.x ="ID",by.y = "Accession")
 
pandFC <- reslutdu[abs(reslutdu$log2FC)>0.5 & reslutdu$Pvalue<0.05,]

pandFC0.5 <- reslutdu[reslutdu$Pvalue<0.05,]

#############

reslutdu1 <- reslutdu[,c("ID","log2FC","Pvalue","KEGG","Taxonomy")]

reslutdu1 <- reslutdu1 %>%
      mutate(Highest=dplyr::case_when(
        reslutdu1$log2FC>0.5 ~ "Up",
        reslutdu1$log2FC < -0.5 ~ "Down",
        TRUE ~ "ins")
      )

######################
dim(pandFC)

unique(unlist(strsplit(na.omit(pandFC$KEGG), ";")))  

######261 proteins（mapping KEGG ko）from 194 genus 

uppandFC <- pandFC[pandFC$log2FC>0,]

unique(uppandFC$Taxonomy)

unique(unlist(strsplit(na.omit(uppandFC$KEGG), ";")))  

#154 Ko from 131 genus

downpandFC <- pandFC[pandFC$log2FC<0,]

unique(downpandFC$Taxonomy)

unique(unlist(strsplit(na.omit(downpandFC$KEGG), ";")))

#139 proteins（mapping KEGG ko）from 106 genus
#139 Ko from 106 genus

#openxlsx::write.xlsx(pandFC,"宏蛋白/宏蛋白pandFC1.xlsx")

genomicDE <- openxlsx::read.xlsx("宏基因组/lefse_result_cut.xlsx")

genomicDEG <- paste0("g__",genomicDE$Taxa)

genomicDEG.pos <- sapply(genomicDEG, function(x) which(grepl(x, pandFC$Taxonomy)))


############################ 

data1 <- pandFC[unlist(genomicDEG.pos),]

data3 <- reslutdu1[which(reslutdu1$ID%in%data1$ID),]


ggplot(reslutdu1) + 
  geom_point(aes(x=log2FC,y=-log10(Pvalue),
                 color=log2FC,
                 size=-log10(Pvalue)))+
  geom_point(data=data3%>%
               tidyr::drop_na()%>%
               dplyr::filter(Highest!="ins")%>%
               dplyr::arrange(desc(abs(log2FC)))%>%
               dplyr::slice(1:150),
             aes(x=log2FC,y=-log10(Pvalue),
                 fill=log2FC,
                 size=-log10(Pvalue)),
             shape=21, show.legend=F, color="#000000")+
  geom_text_repel(data =data3%>%
                    tidyr::drop_na() %>% 
                    dplyr::filter(Highest !="ins")%>%
                    dplyr::arrange(desc(-log2FC))%>%
                    dplyr::slice(1:150)%>%
                    dplyr::filter(Highest=="Up"),
                  aes(x=log2FC,y=-log10(Pvalue), label =ID),
                  box.padding=0.5,
                  nudge_x=0.5,
                  nudge_y=0.2,
                  segment.curvature=-0.1,
                  segment.ncp=3,
                  segment.angle=10,
                  max.overlaps = Inf,
                  direction="y",
                  hjust="left"
  )+
  geom_text_repel(data=data3%>% 
                    tidyr::drop_na() %>% 
                    dplyr::filter(Highest !="ins")%>%
                    dplyr::arrange(desc(log2FC))%>%
                    dplyr::slice(1:150)%>%
                    dplyr::filter(Highest == "Down"),
                  aes(x = log2FC, y = -log10(Pvalue), label =ID),
                  box.padding = 0.5,
                  nudge_x = -0.2,
                  nudge_y = 0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20,
                  direction = "y", 
                  max.overlaps = Inf,
                  hjust = "right"
  ) + 
  scale_color_gradientn(colours = c("#3288bd",  # -24
                                    "#66c2a5",  # -11
                                    "#ffffbf",   # 0
                                    "#f46d43",   # 2.5
                                    "#9e0142"),  # 5
                        values = rescale(c(-24, -3, 0, 2.5, 5), to = c(0, 1))) +
  scale_fill_gradientn(colours = c("#3288bd", 
                                   "#66c2a5", 
                                   "#ffffbf", 
                                   "#f46d43", 
                                   "#9e0142"),
                       values = rescale(c(-24, -3, 0, 2.5, 5), to = c(0, 1)))  +
   
  scale_x_continuous(breaks = c(7,0.5,-0.5, -6, -10, -21, -24), 
                     labels = c("7","0.5", "-0.5","-6", "-10", "", "-24"),
                     limits = c(-24, 7)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 4) + 
  scale_size(range = c(1,7)) + 
  ylim(c(-1, 8)) + 
  theme_bw() + 
  theme(text = element_text(family="serif",size = 8),
        panel.grid = element_blank(),
        legend.background = element_roundrect(color = "#808080", linetype = 1),
        axis.text = element_text(size = 13, color = "#000000"),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
  ) + 
  
  coord_cartesian(clip = "off") + 
  annotation_custom(
    grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
      gp = grid::gpar(lwd = 3, col = "#74add1")
    ), 
    xmin = -2, 
    xmax = -8,
    ymin = 7,
    ymax = 7
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "Ctrl",
      gp = grid::gpar(col = "#74add1")
    ),
    xmin = -2, 
    xmax = -8,
    ymin = 7,
    ymax = 7
  ) +
  annotation_custom(
    grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
      gp = grid::gpar(lwd = 3, col = "#d73027")
    ), 
    xmin = 2, 
    xmax = 7,
    ymin = 7,
    ymax = 7
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "GHPA",
      gp = grid::gpar(col = "#d73027")
    ),
    xmin = 2, 
    xmax = 7,
    ymin = 7,
    ymax = 7
  ) 

 
#########################################################

#########04 Pathway enrichment

#########################################################

Background.num <- length(unique(KEGGdf$ID))

#######################001 Pathway preparation

KO_htable <- load_KO_htable()

mappingpathway <- KO_htable[KO_htable$KO_id%in%unique(KEGGdf$KEGG),]

length(unique(mappingpathway$level3_name))

View(mappingpathway)

Metabolism.table <- mappingpathway[mappingpathway$level1_name%in%"Metabolism",]

length(unique(Metabolism.table$KO_id))

unique(Metabolism.table$level2_name)

table(Metabolism.table$level2_name)

length(unique(Metabolism.table$level3_name))

###########################
 
KEGGpathway <- ReporterScore::load_KO_htable()
 
enrichment.kegg <- function(
    df = df, #gene list
    background.num=9501, #DOI:10.1038/s41467-020-17916-9
    class1="Metabolism", 
    #"Metabolism","Genetic Information Processing" ,"Environmental Information Processing","Cellular Processes"                  
    #"Organismal Systems","Human Diseases","Brite Hierarchies","Not Included in Pathway or Brite"    
    class2="Carbohydrate metabolism",
    adj.p="BH",#c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
    #   "fdr", "none")
    ORA="hypergeometric", #hypergeometric test, Fisher test
    min_exist_KO=3
){df = df #gene list
class1=class1
class2=class2
adj.p=adj.p
ORA=ORA
min_exist_KO=min_exist_KO

KEGGpathway <- ReporterScore::load_KO_htable()

KEGGpathwaysum <- as.data.frame(table(KEGGpathway$level3_name))
colnames(KEGGpathwaysum) <- c("Name","Number")
Komapping <- intersect(df,KEGGpathway$KO_id)
pathwaymapping <- KEGGpathway[which(KEGGpathway$KO_id%in%Komapping),]

if(is.null(class1)&is.null(class2)){
  pathwaymapping <- pathwaymapping
}

if(is.null(class1)==FALSE&is.null(class2)){
  pathwaymapping <- pathwaymapping[which(pathwaymapping$level1_name%in%class1),]
}

if(is.null(class1)&is.null(class2)==FALSE){
  pathwaymapping <- pathwaymapping[which(pathwaymapping$level2_name%in%class2),]
}
if(is.null(class1)==FALSE&is.null(class2)==FALSE){
  pathwaymapping <- pathwaymapping[which(pathwaymapping$level1_name%in%class1),]
  
  pathwaymapping <- pathwaymapping[which(pathwaymapping$level2_name%in%class2),]
}

enrichres <- as.data.frame(table(pathwaymapping$level3_name))

colnames(enrichres) <- c("Name","exist_k")
enrichres <- merge(KEGGpathwaysum,enrichres,by="Name")

if(is.null(min_exist_KO)){
  min_exist_KO <- 3
  enrichres <- enrichres[enrichres$exist_k>min_exist_KO,]
}else{
  enrichres <- enrichres[enrichres$exist_k>min_exist_KO,]
}


M <-background.num
n <- length(df) #   

enrichres$Map_ID <- c(rep("NA",nrow(enrichres)))
enrichres$hit_ko <- c(rep("NA",nrow(enrichres)))
enrichres$p.value <- c(rep("NA",nrow(enrichres)))
for(i in 1:nrow(enrichres)){
  hitko <- pathwaymapping$KO_id[which(pathwaymapping$level3_name%in%enrichres$Name[i],)]
  
  mapid<- pathwaymapping$level3_id[which(pathwaymapping$level3_name%in%enrichres$Name[i],)]
  enrichres$Map_ID[i] <- unique(mapid) 
  enrichres$hit_ko[i] <- paste(hitko, collapse = ",")
  N <- enrichres[i,"Number"]#  
  k <- enrichres[i,"exist_k"]
  
  if(ORA=="hypergeometric"){
    p_value <- phyper(k - 1, N, M - N, n, lower.tail = FALSE)
  }
  if(ORA=="Fisher"){
    fisher_result <- fisher.test(matrix(c(k, N-k, n-k, M-N-n+k), nrow=2))
    p_value <- fisher_result$p.value
  }
  enrichres$p.value[i] <- p_value
  enrichres <- enrichres
}
enrichres$adj.p <- p.adjust(enrichres$p.value,method = adj.p)

return(enrichres)
}
 
###########################002 KEGG enrichment

genelist <- unique(unlist(strsplit(na.omit(pandFC0.5$KEGG), ";"))) 


enrichkegg <- enrichment.kegg(df = genelist, #gene list
                              
                              background.num=Background.num, 
                              
                              class1=c("Metabolism","Organismal Systems","Human Diseases"), 
                              
                              class2=c("Carbohydrate metabolism",
                                       "Energy metabolism",
                                       "Lipid metabolism",
                                       "Nucleotide metabolism",
                                       "Amino acid metabolism",
                                       "Metabolism of other amino acids",
                                       "Glycan biosynthesis and metabolism",
                                       "Metabolism of cofactors and vitamins",
                                       "Metabolism of terpenoids and polyketides",
                                       "Biosynthesis of other secondary metabolites",
                                       "Xenobiotics biodegradation and metabolism",
                                       "Immune system",
                                       "Endocrine system",
                                       "Immune disease",
                                       "Endocrine and metabolic disease"
                              ),
                              
                              adj.p="BH",#c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",   "fdr", "none")
                              ORA="hypergeometric", 
                              #hypergeometric test, Fisher test
                              min_exist_KO=0)

enrichkegg$ratio <- enrichkegg$exist_k/enrichkegg$Number

enrichkegg <- enrichkegg[enrichkegg$adj.p<0.05,]

b <- unique(KEGGpathway[,c("level1_name","level2_name","level3_id")])
enrichkegg <- merge(enrichkegg,b,by.y="level3_id",by.x="Map_ID")
View(enrichkegg)
 
######### 
 
KEGGdf2 <- merge(KEGGdf1,KEGGFC[,c("ID","Num_GHPA","Num_Ctrl","Agv_GHPA","Agv_Ctrl","log2FC")],by="ID")

dupID <- unique(KEGGdf2$KEGG[duplicated(KEGGdf2$KEGG)]) 

nondupID <- setdiff(unique(KEGGdf2$KEGG),dupID)

dupdf <- KEGGdf2[which(KEGGdf2$KEGG%in%dupID),]
###
dupdf1 <- dupdf[dupdf$Pvalue<0.05,]

dupname <- unique(dupdf1$KEGG[duplicated(dupdf1$KEGG)])

dupdf1_non <- dupdf1[which(dupdf1$KEGG%in%setdiff(dupdf1$KEGG,dupname)),]

for(ip in seq(dupname)){
  
  temp <- dupdf1[which(dupdf1$KEGG%in%dupname[ip]),]
  
  if(any(temp$Num_GHPA<3)|any(temp$Num_Ctrl<3)){
    
    keeprow <- temp[temp$Num_GHPA>2&temp$Num_Ctrl>2,]
    
    keeprow2 <- keeprow[which.max(abs(keeprow$log2FC)),]
    
    dupdf1 <- rbind(dupdf1_non,keeprow2)
  }
  if(all(temp$Num_GHPA<3)|all(temp$Num_Ctrl<3)){
    
    keeprow2 <- temp[which.max(abs(temp$log2FC)),]
    
    dupdf1 <- rbind(dupdf1_non,keeprow2)
  }
  if(all(temp$Num_GHPA>2)|all(temp$Num_Ctrl>2)){
    
    keeprow2 <- temp[which.max(abs(temp$log2FC)),]
    
    dupdf1 <- rbind(dupdf1_non,keeprow2)
  }
}
#####
dupdf2 <- dupdf[dupdf$Pvalue>0.05,] 

a <- which(dupdf2$KEGG%in%intersect(dupdf1$KEGG,dupdf2$KEGG))

dupdf2 <- dupdf2[-a,]

dupdf2 <- dupdf2 %>%
  group_by(KEGG) %>%
  slice(which.min(Pvalue))# Keep the row with min p Value 

uniqueKEGGdf <- rbind(dupdf1,dupdf2,KEGGdf2[which(KEGGdf2$KEGG%in%nondupID),] )

uniqueKEGGdf <- as.data.frame(uniqueKEGGdf)

uniqueKEGGdf1 <- uniqueKEGGdf %>% 
  mutate(Highest=dplyr::case_when(
    
    uniqueKEGGdf$log2FC>0 ~ "GHPA",
    
    uniqueKEGGdf$log2FC<0 ~ "ctrl")
  )

colnames(uniqueKEGGdf1)[22:29] <- c("KO_id","p.value","Num_GHPA","Num_Ctrl","average_GHPA","average_Ctrl","diff_mean", "Highest")

View(uniqueKEGGdf1[uniqueKEGGdf1$p.value<0.05,])

ko_stat=pvalue2zs(uniqueKEGGdf1[,22:29],mode="directed")

rownames(ko_stat) <- ko_stat$KO_id

reporter_s=get_reporter_score(ko_stat,
                              #modulelist=mydat,
                              min_exist_KO =2)
modulelist <- get_modulelist(feature = "ko", 
                             type = "pathway", 
                             verbose = FALSE, 
                             chr = TRUE)

rownames(uniqueKEGGdf1) <- uniqueKEGGdf1$KO_id

kodf <- uniqueKEGGdf1[,2:21]

Group_tab <- data.frame(Group=c(rep("GHPA",10),rep("Ctrl",10)))

rownames(Group_tab) <- colnames(KO_abundance)

Group_tab$Group <- factor(Group_tab$Group)

######## create modulelist

KEGGpathway <- ReporterScore::load_KO_htable()

modulelist <- KEGGpathway[,c("level1_name","level3_id","level3_name","KO_id" )]

modulelist1 <- unique(modulelist[,c("level1_name","level3_id","level3_name")])

for (ii in seq(modulelist1$level3_id)) {
  
  mapid <- modulelist1$level3_id[ii]
  
  konum <- get_features(map_id = mapid, ko_stat = NULL, modulelist = NULL)
  
  modulelist1$KOs[ii] <- paste(konum, collapse = ",")
  
  modulelist1$K_num[ii] <- length(konum)
}

modulelist1 <- as.data.frame(modulelist1)

colnames(modulelist1) <- c("Class","id","Description","KOs","K_num")

################create a list

reporter_res <- GRSA(kodf, "Group", Group_tab,
                     mode = "directed",
                     method = "wilcox.test", perm = 999,
                     type = "pathway", feature = "ko",
                     modulelist=modulelist1
)

reporter_res1 <- reporter_res

reporter_res1$kodf <- kodf

reporter_res1$ko_stat <- ko_stat

reporter_res1$reporter_s <- reporter_s

reporter_res1$modulelist <- modulelist1

reporter_res1$metadata <- Group_tab

plot_report(reporter_res1,
            rs_threshold=c(1.64,-1.64),#1.64
            y_text_size=10,str_width=40)

plot_report_bar(reporter_res, rs_threshold = c(0.3,-0.3), facet_level = TRUE)

plot_KOs_network(reporter_res1,
                 map_id = c("map00220","map00430"),#enrichkegg$Map_ID,
                 kos_color = c(Depleted = "orange", 
                               Enriched = "orange", 
                               None = "grey", 
                               Significant ="red2", 
                               Pathway = "#80b1d3"),
                 pathway_description = T, 
                 kos_description = F,
                 main = "", mark_module = TRUE
)

plot_KOs_box(reporter_res1,map_id = "map00430",only_sig = TRUE)

koins <- NULL

koins <- plot_KOs_box.YHC(reporter_res,map_id = c("map00220","map00430"),only_sig = TRUE)

a <- c("ID","Num_GHPA","Num_Ctrl","Agv_GHPA","Agv_Ctrl","FC","log2FC","Pvalue","KEGG.y")
 
pandFC1 <- pandFC%>%select(a)

pathP <- pandFC1[which(pandFC1$KEGG.y%in%koins),]

set.seed(1234)

gsea_res <- KO_gsea(reporter_res, weight = "Z_score")

plot_KOs_box(reporter_res,map_id = mapid,only_sig = TRUE,
             
             KO_description = T,str_width = 70)

KOINS <- plot_KOs_box.YHC(reporter_res,map_id = mapid,only_sig = TRUE,
                          
                          KO_description = T)
 
#########################################
data4 <- data3 

data4$Taxonomy <- sapply(str_split(data3$Taxonomy,";"), "[",7)

data4$Taxonomy <- sapply(str_split(data4$Taxonomy,"g__"), "[",2)

data4 <- data4 %>%
  mutate(annotation = paste(Taxonomy,":",KEGG))
 
data4KEGG <- unique(unlist(strsplit(na.omit(data4$KEGG), ";"))) 


##########

KEGGpathway <- ReporterScore::load_KO_htable()
###########################KEGG enrichment
data4kegg <- enrichment.kegg(df = data4KEGG, #gene list
                              background.num=backgroud.num, 
                              #DOI:10.1038/s41467-020-17916-9
                              class1=c("Metabolism"), 
                              class2=NULL,
                              #class2="Amino acid metabolism",
                              adj.p="BH",#c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
                              ORA="hypergeometric", 
                              #hypergeometric test, Fisher test
                              min_exist_KO=0)


meta_data4 <- unique(unlist(strsplit(data4kegg$hit_ko, ","))) 

data4kegg <- data4kegg%>% 
  dplyr::arrange(desc(-adj.p))
 
############################# 
data4kegg1 <- data4kegg[data4kegg$adj.p<0.05,]

split(data4kegg1$hit_ko,",")
data4keggenrich <- unlist(strsplit(data4kegg1$hit_ko, ","))   

data <- data4[which(data4$KEGG%in%unique(data4keggenrich)),]# 

################GSEA
 
gmt <- KEGGpathway[,c("level3_name","KO_id")]
colnames(gmt) <- c("term","gene")
data <- data %>% 
  arrange(desc(log2FC))

geneList = data$log2FC  
names(geneList) <- data$KEGG  
library(clusterProfiler)

gsea_custom_result <-  clusterProfiler::GSEA(geneList, 
                            TERM2GENE = gmt, 
                            minGSSize = 1, 
                            pvalueCutoff = 0.99, 
                            verbose = F)

gsea.out.df <- gsea_custom_result@result
View(gsea.out.df)
gsea.out.df$core_enrichment

data4score <- merge(data4kegg1,gsea.out.df[,c("ID","enrichmentScore","NES")],by.x="Name",by.y="ID")

colnames(data4score)[c(6,7)] <- c("hypergeometric P value","BH adj P")
 
###########################

#########################

library(networkD3)
 
soure <- KEGGpathway[which(KEGGpathway$KO_id%in%data4keggenrich),]
soure <- soure[which(soure$level1_name%in%"Metabolism"),]
soure <- merge(soure,data4,by.x = "KO_id",by.y="KEGG")
 

soure1 <- soure[,c("level3_name","annotation")]

soure1$W <- rep(1,nrow(soure1))

nodes <- data.frame(name = c(unique(soure1$level3_name),unique(soure1$annotation)))
links <- data.frame(source = match(soure1$annotation, nodes$name) - 1,
                    target = match(soure1$level3_name, nodes$name) - 1,
                    value = soure1$W)

 
sankey <- sankeyNetwork(Links = links, Nodes = nodes, 
              Source = "target", Target = "source", 
              Value = "value", NodeID = "name", 
              fontSize = 12, nodeWidth = 30)

sankey 

htmlwidgets::saveWidget(sankey, "sankey_3.html")

library(ggalluvial)

 ggplot(soure1, aes(axis1 = annotation, axis2 = level3_name, y = W)) +
   geom_alluvium(aes(fill = annotation), width = 1/12) +
   geom_stratum(aes(fill = level3_name),width = 1/12, fill = "black", color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()

 ###################
 ALD <- pandFC[which(pandFC$ID%in%"A_4_k97_147147_1_1"),2:21]
 
 ALD <- t(ALD)
ALD[which(grepl("NaN",ALD)),]<- 0
 
 group <- c(rep("GHPA",10),rep("Ctrl",10))
 
 ALD <- as.data.frame(ALD)
 ALD$Group <- group
 
 ALD$Group <- factor(group, levels=c("Ctrl","GHPA"))
 
 colnames(ALD)[1] <- "ALD"
 
 library(ggpubr)
 
 sig<-compare_means(ALD~Group, 
                    data=ALD,
                    method = "wilcox.test",
                    p.adjust.method = "BH",) #t.test
#                            p    p.adj   p.format    p.signif method  
#  
#   1 ALD   Ctrl   GHPA   0.0433 0.043 0.043    *        Wilcoxon
 
 y.max <-max(ALD$ALD)+1.3
 
 lab.max <- y.max-0.2
 
 ggplot(data=ALD, 
        aes(x=Group,y=ALD))+
   annotate("text", x=2, y=lab.max, label= sig$p.signif)+
   geom_boxplot(aes(color=Group),size=0.8, outliers=F, position=position_dodge(1))+
   geom_jitter(alpha=0.5,aes(color=Group),
               position=position_jitterdodge(jitter.width = 0.3, 
                                             jitter.height = 0, 
                                             dodge.width = 0.4
               ))+
   scale_color_manual(values = c("#5E7BBB","#D95201"))+
   # facet_wrap(~ Name, scales="free")+
   theme_bw()+
   ylab("Abundance levels")+
   labs(title="Bilophila:K00259")+
   scale_x_discrete(labels=c("Ctrl (N=10)","GHPA (N=10)"))+
   theme(text = element_text(family="serif",size = 8),#"Times New Roman","serif"
         panel.grid = element_blank(),
         axis.ticks.x.bottom = element_blank(),
         axis.title.x = element_blank(),
         plot.title = element_text(hjust=0.5),
         axis.text.x = element_text(size=8),
         axis.text.y = element_text(size=10),
         legend.position = "none"
   )
 
 
 
#########################################################################################################

#####################05 spearman based on common genuse of metagenomics and metaproteomics
 
######################################################################################################### 
 
 alldata <- clean_df
 
 for (i in seq_len(nrow(alldata))) {
   row_name <- rownames(alldata)[i]
   if (grepl(";", row_name)) {  
     
     split_names <- strsplit(row_name, ";")[[1]] 
     
     new_df <- matrix(0, ncol = ncol(alldata), nrow = length(split_names))   
     
     new_df <- as.data.frame(new_df)
     
     for (ith in 1:length(split_names)) {
       new_df[ith,] <-alldata[i,]
     }
     colnames(new_df) <- colnames(alldata)
     
     rownames(new_df) <- split_names
 
     df.cb <- rbind(alldata, new_df)
     
   }else{
     df.cb<-alldata 
   }
   alldata <- df.cb
 }  
 
 
 DEdata <- alldata[which(rownames(alldata)%in%intersect(data$ID,rownames(alldata))) ,]
 
 DEdata$ID <- rownames(DEdata)
 
 DEdata <- merge(data[c(1,5,7)],DEdata,by="ID")
 
 DEdata1 <- DEdata[,-c(1:2)]
 
 rownames(DEdata1) <- DEdata1$annotation
 
 DEdata1 <- DEdata1[,-1]#### 
 
 genus.pt <- read.table("tax_g_RPKM.xls",header = T)
 
 genus.pt$genus <-unlist( lapply(1:nrow(genus.pt), function(t) {
   
   a <- unlist(strsplit(genus.pt$Taxonomy[t], split = ";"))[7]
   
 }))
 
 genus.p.matrix <- genus.pt[,2:21]
 
 rownames(genus.p.matrix) <- genus.pt$genus
 
 b <- grep("g__unclassified", rownames(genus.p.matrix) )
 
 genus.p.matrix <- genus.p.matrix[-b,]
 
 rownames(genus.p.matrix) <- sapply(str_split(rownames(genus.p.matrix),"g__"), "[",2)
 
 DEgenus.matrix <- genus.p.matrix[which(rownames(genus.p.matrix)%in%data$Taxonomy),]
 
 DEgenus.matrix <- DEgenus.matrix[,colnames(DEdata1)]
 
 P.cor<- corr.test(t(DEgenus.matrix),
                   t(DEdata1),
                   method = "spearman",adjust="BH")
 
 par(family= "serif")# 
 corrplot(t(P.cor$r), 
          p.mat = round(t(P.cor$p),2), 
          sig.level = 0.05,
          pch.cex=1,
          tl.srt=45,
          addgrid.col="grey90",
          insig = 'label_sig', 
          tl.col = "black",
          method = 'square',  
          col = rev(COL2(n=100)), 
          number.cex = 0.4, 
          addCoef.col = NULL,
          tl.cex = 0.7,
          cl.cex = 0.7)
 
 



