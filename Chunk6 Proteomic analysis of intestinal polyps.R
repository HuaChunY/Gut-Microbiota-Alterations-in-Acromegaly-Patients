####################################################################################################################################################################
#      ☆        % Project: Gut Microbiota Alterations in Acromegaly Patients Are Associated with Neutrophil Depletion-Induced Inflammation #
#   ☆ \|/ ☆    % Author: HuaChun Yin                                         
#  ☆  \|/  ☆   % Date: Oct. 14th, 2024                            
# ☆   \|/   ☆  %                                                          
#  ☆  \|/  ☆   % Environment:   R version 4.4.2           
#  ☆ __|__ ☆   % Platform: Mac-IOS(64-bit)                                  
#                % CHUNK5:this script includes proteomics mapping, missing value imputation and normalization             
################################################################################################################################################################### 


packages <- c("dplyr","doParallel","multiUS","vsn", 
              "ggplot2","ggpubr","tidyr","ggfun",
              "psych","plotly",
              "stringr","corrplot",
              "RColorBrewer",
              "ggsci","microeco",
              "tidyverse","magrittr",
              "data.table","foreach", 
              "ggalluvial","PMA",
              "ggrepel","scales","hdi","ropls",
              "stabs","Rtsne","ropls")  
#check.packages(packages)

lapply(packages, library, character.only = TRUE)

#########################################################

#########01 Data filtering

#########################################################

Na.num <- 3 #蛋白出现在30%(3, 至少4个样本)

protein<- openxlsx::read.xlsx("protein.xlsx")

KEGGIDtrans1<- openxlsx::read.xlsx("proteins_anno_detail.xlsx")

KEGGIDtrans <- KEGGIDtrans1[,c(1,9)]

protein1 <- protein[,c(1,7:10,14:15,11:13,16:18)]

rownames(protein1) <- protein1$PG.ProteinGroups

protein1 <- protein1[,2:13]


#filtering
keeprow <- c()

removerow <- c()

for (i in 1:nrow(protein1)){
  
  if( 12-sum(nchar(gsub("NaN", "", protein1[i,]))==0) >Na.num){# 70% were detected at least
    
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
 
clean_df[] <- lapply(clean_df, as.numeric)


####################################################################

######02 KNN miss value imputation

################################################################ 

clean_df1 <- multiUS::KNNimp(clean_df, k = 5) #k from DOI:10.1038/s41467-020-17916-9

clean_df2 <- as.data.frame(clean_df1)

clean_df2$ID <- rownames(clean_df2)

#####################
for (i in seq_len(nrow(clean_df2))) {
  
  row_name <- clean_df2$ID[i]
  
  if (grepl(";", row_name)) {  #  
    split_names <- strsplit(row_name, ";")[[1]]  # 
    
    new_df <- matrix(0, ncol = ncol(clean_df2), nrow = length(split_names))  #  
    
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

nonseicolondf <- clean_df2[-rows_with_semicolon,] 
 
length(intersect(rownames(nonseicolondf),KEGGIDtrans$Accession_id))

#######################KEGG ID mapping

KEGGdf <- merge(nonseicolondf,KEGGIDtrans,by.x= "ID",by.y="Accession_id")

background.num <- length(unique(KEGGdf$ID))

background.num
#[1] 8290
 
length(unique(KEGGdf$KO_id))

KO_htable <- load_KO_htable()

mappingpathway <- KO_htable[KO_htable$KO_id%in%unique(KEGGdf$KO_id),]

length(unique(mappingpathway$level3_name))

View(mappingpathway)

Metabolism.table <- mappingpathway[mappingpathway$level1_name%in%"Metabolism",]

length(unique(Metabolism.table$KO_id))

unique(Metabolism.table$level2_name)

table(Metabolism.table$level2_name)

length(unique(Metabolism.table$level3_name))

colnames(KEGGdf)

e.p <-2:7
c.p <- 8:13

KEGGdf[,c(e.p,c.p)] <- lapply(KEGGdf[,c(e.p,c.p)], as.numeric)
##################TSNE
delet.col <- c("ID","KO_id")
 
KEGGdf1 <- KEGGdf %>%
  select(-ID, -KO_id)


group <- c(rep("GHPA",6),rep("Ctrl",6))

##########

initial_value<-10
theta_value <- 0.1 #[0,1]
plex_value <- 3 # < 1/3 samples


#######

set.seed(5)
tsne_out<-Rtsne(t(KEGGdf1),
                dims = 2,initial_dims = initial_value,
                pca = FALSE,
                perplexity = plex_value,
                theta = theta_value, 
                max_iter = 1000
)

tsne <- data.frame(tSNE1=tsne_out$Y[,1],
                   tSNE2=tsne_out$Y[,2],
                   Group=group,
                   ID=colnames(KEGGdf1))
 
ggplot(tsne,aes(tSNE1,tSNE2))+
  geom_point(aes(color=group),size=1.5)+
  scale_color_lancet()+
  scale_color_manual(values=c("#D95201","#5E7BBB"))+
  theme_bw()+
  xlim(-220,220)+#300
  ylim(-210,210)+#150
  theme(plot.margin = unit(rep(1.5,4),"lines"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.9,0.2), # 
        legend.background = element_rect(size = 1, colour = "white"))

tsne_out$Y[1,]


tsne_df <- data.frame(tsne_out$Y, group = group)

colnames(tsne_df) <- c("Dim1", "Dim2", "group")

# Calculate p-value using t-test for one of the dimensions (e.g., Dim1)
group1 <- tsne_df %>% filter(group == "GH") %>% pull(Dim1)

group2 <- tsne_df %>% filter(group == "nonGH") %>% pull(Dim1)

# Perform t-test
t_test_result <- t.test(group1, group2)

# Extract p-value
p_value <- t_test_result$p.value

p_value
#[1] 0.02120031
#################


#########Normalization
 
par(mfrow = c(1, 2))
boxplot(log2(KEGGdf[,c(e.p,c.p)]),main = "Before normalization",  outline =F,names=F)

fit<- vsn2(as.matrix(KEGGdf[,c(e.p,c.p)]))

#meanSdPlot(fit)
 
ynorm = predict(fit, as.matrix(KEGGdf[,c(e.p,c.p)]))

KEGGdfnorm = justvsn(ynorm)

boxplot(KEGGdfnorm,main = "After normalization", outline =F,names=F)
par(mfrow = c(1, 1))  

pmad <- mean(mad(KEGGdfnorm))#less than 0.14 

print(pmad)
#0.022
 
data_long <- gather(as.data.frame(KEGGdfnorm), key = "Variable", value = "Value")

data_long$Variable <- rep("A",nrow(data_long))

library(ggplot2)
ggplot(data_long, aes(x = Variable, y = Value)) +
  #geom_boxplot()+
  geom_jitter(width=0.3,size=0.5, color = "#5E7BBB", alpha = 0.6)+
  labs(title = paste("PMAD:",round(pmad,3))) +
  theme(text = element_text(family="serif",size = 10),#"Times New Roman","serif"
        panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        legend.position = "none"
  )
#DOI: 10.1093/bib/bbz061
#The RLA plots before and after the VSN normalization and the PMAD for six benchmark data sets. 
#(A) The RLA plots for unnormalized intensities; (
#(B) the RLA plots for normalized intensities; 
#(C) the distributions of PMAD values. 
KEGGdf1 <- KEGGdf

KEGGdf1[,c(e.p,c.p)] <- KEGGdfnorm

KEGGdf1$"Pvalue"<- apply(KEGGdf1[,c(e.p,c.p)], 1, function(row){
  wilcox.test(row[1:6], row[7:12])$p.value
})

#Statistical significance of the difference between groups was evaluated using Mann–Whitney U test (nonparametric test of the null hypothesis which is suitable for data that does not pass normality test), unless otherwise indicated.
#DOI:10.1038/s41467-020-17916-9
 
KEGGdfduplicated_rows <- KEGGdf1[duplicated(KEGGdf1$KO_id) | duplicated(KEGGdf1$KO_id, fromLast = TRUE), ]

KEGGdfduplicated_p <- KEGGdfduplicated_rows[KEGGdfduplicated_rows$Pvalue<0.05,] 

KEGGdfduplicated_pdu <- KEGGdfduplicated_p[duplicated(KEGGdfduplicated_p$KO_id) | duplicated(KEGGdfduplicated_p$KO_id, fromLast = TRUE), ]

length(unique(KEGGdfduplicated_pdu$KO_id))
#####################

####################################################################

######03 Differentially abundant proteins 

################################################################

rawFC <- clean_df
rawFC$Num_GHPA <- rep(NA,nrow(rawFC))

rawFC$Num_Ctrl <- rep(NA,nrow(rawFC))

rawFC$Agv_GHPA <- rep(NA,nrow(rawFC))

rawFC$Agv_Ctrl <- rep(NA,nrow(rawFC))

e.na <-1:6
c.na <- 7:12

for (i in 1:nrow(rawFC)){
  
  epsun <- 6-sum(nchar(gsub("NaN", "", clean_df[i,e.na]))==0) 
  
  cpsun <- 6-sum(nchar(gsub("NaN", "", clean_df[i,c.na]))==0) 
  
  allsum <- 12-sum(nchar(gsub("NaN", "", clean_df[i,]))==0) 
  
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
 
rawFC1_semicolon <- grep(";", rawFC1$ID)

nonseicolonrawFC1 <- rawFC1[-rawFC1_semicolon,]

KEGGFC <- merge(nonseicolonrawFC1,KEGGIDtrans,by.x= "ID",by.y = "Accession_id")

#######################################
pandFC <- merge(KEGGdf1,KEGGFC,by= "ID")

pandFC <- merge(pandFC,KEGGIDtrans1,by.x= "ID",by.y="Accession_id")

pandFC <- pandFC %>% 
  mutate(Highest=dplyr::case_when(
    
    pandFC$log2FC>log2(2) ~ "Up",
    
    pandFC$log2FC<log2(0.5) ~ "Down",
    
    TRUE ~ "ins")
  )
 
pandFC.res <- pandFC[-which(pandFC$Highest%in%"ins"),]

pandFC.res <- pandFC.res[pandFC.res$Pvalue<0.05,]
 
##################
 
S100A8.df <- pandFC.res[which(pandFC.res$symbol%in%"S100A8"),c(2:13)]
 
S100A8.df <- as.data.frame(t(S100A8.df))
 
S100A8.df$Group <- factor(group, levels=c("Ctrl","GHPA"))

colnames(S100A8.df)[1] <- "S100A8"
  
sig<-compare_means(S100A8~Group, 
                   data=S100A8.df,
                   method = "wilcox.test",
                   p.adjust.method = "BH",) #t.test

y.max <-max(S100A8.df$S100A8)+10000

lab.max <- y.max-0.2

p1 <-  ggplot(data=S100A8.df, 
                    aes(x=Group,y=S100A8))+
  annotate("text", x=2, y=lab.max, label= sig$p.signif)+
  geom_boxplot(aes(color=Group),size=0.8, outliers=F, position=position_dodge(1))+
  geom_jitter(alpha=0.5,aes(color=Group),
              position=position_jitterdodge(jitter.width = 0.3, 
                                            jitter.height = 0, 
                                            dodge.width = 0.4
              ))+
  scale_color_manual(values =c("#5E7BBB","#EF9703"))+
  # facet_wrap(~ Name, scales="free")+
  theme_bw()+
  ylab("Abundance levels")+
  labs(title="S100A8")+
  scale_x_discrete(labels=c( "Ctrl (N=6)" , "GHPA (N=6)"))+
  theme(text = element_text(family="serif",size = 8),#"Times New Roman","serif"
        panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=10),
        legend.position = "none"
  )

 
 #########################################################
 
 ##################
 
 PADI4.df <- pandFC.res[which(pandFC.res$symbol%in%"PADI4"),c(2:13)]
 
 PADI4.df <- as.data.frame(t(PADI4.df))
 
 PADI4.df$Group <- factor(group, levels=c("Ctrl","GHPA"))
 
 colnames(PADI4.df)[1] <- "PADI4"
 
 PADI4.df <- as.data.frame(na.omit(PADI4.df))
 
 sig<-compare_means(PADI4~Group, 
                    data=PADI4.df,
                    method = "wilcox.test",
                    p.adjust.method = "BH",) #t.test
 
 y.max <-max(PADI4.df$PADI4)+1000
 
 lab.max <- y.max-0.2
 
p2 <- ggplot(data=PADI4.df, 
        aes(x=Group,y=PADI4))+
   annotate("text", x=2, y=lab.max, label= sig$p.signif)+
   geom_boxplot(aes(color=Group),size=0.8, outliers=F, position=position_dodge(1))+
   geom_jitter(alpha=0.5,aes(color=Group),
               position=position_jitterdodge(jitter.width = 0.3, 
                                             jitter.height = 0, 
                                             dodge.width = 0.4
               ))+
   scale_color_manual(values =c("#5E7BBB","#EF9703"))+
   # facet_wrap(~ Name, scales="free")+
   theme_bw()+
   ylab("Abundance levels")+
   labs(title="PADI4")+
   scale_x_discrete(labels=c( "Ctrl (N=6)" , "GHPA (N=6)"))+
   theme(text = element_text(family="serif",size = 8),#"Times New Roman","serif"
         panel.grid = element_blank(),
         axis.ticks.x.bottom = element_blank(),
         axis.title.x = element_blank(),
         plot.title = element_text(hjust=0.5),
         axis.text.x = element_text(size=8),
         axis.text.y = element_text(size=10),
         legend.position = "none"
   )
 p1/p2
 #########################################################
 
 #########04 Pathway enrichment
 
 #########################################################
 
gmt <- openxlsx::read.xlsx("GSEApathway.xlsx",sheet=1)

#pandFC.res <- openxlsx::read.xlsx("pandFC.xlsx")

data <- pandFC.res %>% 
  arrange(desc(log2FC))

geneList = data$log2FC  

names(geneList) <- data$symbol



gsea_custom_result <-  clusterProfiler::GSEA(geneList[-which(is.na(names(geneList) ))], 
                            TERM2GENE = gmt, 
                            minGSSize = 1, 
                            pvalueCutoff = 0.99, 
                            verbose = F)

gsea.out.df <- gsea_custom_result@result

View(gsea.out.df)

gsea.out.df$core_enrichment
################

pathway <- openxlsx::read.xlsx("toppgene.xlsx",sheet = 1)

pathway <- pathway[1:10, ]

CC <- openxlsx::read.xlsx("CC.xlsx",sheet = 2)

CC <- CC[1:10, ]

pathway <- merge(pathway,gsea.out.df[,c("ID", "enrichmentScore","NES")],by.x ="Name",by.y="ID" )

CC <- merge(CC,gsea.out.df[,c("ID", "enrichmentScore","NES")],by.x ="Name",by.y="ID" )

colnames(pathway)[c(6:7)] <- c("AdjP","count")

colnames(CC)[c(5:6)] <- c("AdjP","count")
  
P1 <- ggplot(CC,aes(NES,reorder(Name,count)))+
  geom_point(aes(size=count,color=-log10(AdjP)))+
  scale_color_gradient(low = "#c7e9c0", high = "#006d2c") +
  labs(color=expression(-log10(AdjP)),
       size="Number",
       y=NULL,
       x="NES Score",
       fill="-log10(Adj P)")+
  #geom_vline(xintercept = c(1,-1), linetype = "dashed", color = "red")+
  labs(title="GO: Cellular Component")+ 
  theme_bw()+
  theme(axis.text.x = element_text(size=10,color = "black"),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=10),
        text=element_text(family="serif",size=10))

P2 <- ggplot(pathway,aes(NES,reorder(Name,count)))+
  geom_point(aes(size=count,color=-log10(AdjP)))+
  scale_color_gradient(low = "#c7e9c0", high = "#006d2c") +
  labs(color=expression(-log10(AdjP)),
       size="Number",
       y=NULL,
       x="NES Score",
       fill="-log10(Adj P)")+
  #geom_vline(xintercept = c(1,-1), linetype = "dashed", color = "red")+
  labs(title="Pathway")+ 
  theme_bw()+
  theme(axis.text.x = element_text(size=10,color = "black"),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=10),
        text=element_text(family="serif",size=10))
 
#####################
result <- openxlsx::read.xlsx("network.xlsx",sheet = 3)
 
mytheme2 <- theme_classic() +
  theme(axis.text.y = element_text(family = "serif",size = 10,color = "#050505"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.border = element_blank(),
        legend.title.position = "left",
        legend.key.height=unit(0.5, "cm"),
        legend.key.width=unit(0.4, "cm"),
        legend.title = element_text(family = "serif",size = 9,hjust = 0.5, angle = 90),
        axis.line = element_line(linewidth = 0.5,lineend = "square"),
        plot.margin=unit(x=c(0,0.2,0.1,0.2),units="inches"))
 
result$Pathway <- factor(result$Pathway,levels =rev(result$Pathway))

f1 <- ggplot(result[1:10,], aes(x = Count,y = Pathway,fill = -log10(P.value))) +
  geom_col(color = "white",width=0.85,linewidth=0.7) +
  labs(x="Count",y="",fill="-log10(P value)")+
  scale_color_manual(values = "#D95201")+
  scale_fill_gradient(low = "#ffffff",high = "#D64515",
                      limits = c(1,2.6),breaks = c(1,1.5,2,2.5))+
  scale_x_continuous(breaks = seq(0, 8, by = 2),
                     limits = c(0, 8),expand = c(0,0)) +
  
  mytheme2
f1

f2 <- ggplot(result[11:20,], aes(x = Count,y = Pathway,fill = -log10(P.value))) +
  geom_col(color = "white",width=0.85,linewidth=0.7) +
  labs(x="Count",y="",fill="-log10(P value)")+
  scale_color_manual(values ="#5E7BBB")+
  scale_fill_gradient(low = "#ffffff",high = "#425C9E",
                      limits = c(1,7.1),breaks = c(1,2,3,4,5,6,7))+
  scale_x_continuous(breaks = seq(0, 8, by = 2),
                     limits = c(0, 8),expand = c(0,0)) +
  
  mytheme2
f2


##############
 

