####################################################################################################################################################################
#      ☆        % Project: Gut Microbiota Alterations in Acromegaly Patients Are Associated with Neutrophil Depletion-Induced Inflammation #
#   ☆ \|/ ☆    % Author: HuaChun Yin                                         
#  ☆  \|/  ☆   % Date: Apr. 4th, 2025                                  
# ☆   \|/   ☆  %                                                          
#  ☆  \|/  ☆   % Environment:   R version 4.4.2           
#  ☆ __|__ ☆   % Platform: Mac-IOS(64-bit)                                  
#                % CHUNK4:this script includes metagenomics and metaproteomics                       
################################################################################################################################################################### 

packages <- c("dplyr","doParallel",  "patchwork",
              "ggplot2","ggpubr",
              "psych","plotly","tidyr",
              "stringr","corrplot",
              "RColorBrewer",
              "ggsci","microeco",
              "tidyverse","magrittr",
              "data.table","foreach",
              "lme4","nlme",
              "factoextra","vegan",
              "ggalluvial","PMA",
              "ggrepel","scales","hdi","ropls",
              "stabs","Rtsne","ropls","MicrobiotaProcess","dbplyr") ## package methods is not loaded by default by RScript. 
#check.packages(packages)

lapply(packages, library, character.only = TRUE)
#######################

###########################################################################
 
############# 01 Spearman's rank correlation coefficient

###########################################################################


#####001 Metaproteomic based on genus

protein<- openxlsx::read.xlsx("metaprotein/protein.xlsx")

protein1 <- protein[,c(2,7:26)]


for (i in seq_len(nrow(protein1))) {
  row_name <- protein1$PG.ProteinAccessions[i]
  if (grepl(";", row_name)) {  n
    split_names <- strsplit(row_name, ";")[[1]]  
    new_df <- matrix(0, ncol = ncol(protein1), nrow = length(split_names)) 
    new_df <- as.data.frame(new_df)
    for (ith in 1:length(split_names)) {
      new_df[ith,] <-protein1[i,]
    }
    colnames(new_df) <- colnames(protein1)
    new_df[,"PG.ProteinAccessions"] <- split_names
    
   n
    df.cb <- rbind(protein1, new_df)
    
  }else{
    df.cb<-protein1 
  }
  protein1 <- df.cb
}    

rows_with_semicolon <- grep(";", protein1$PG.ProteinAccessions)

protein2 <- protein1[-rows_with_semicolon,]#keep the uniqe ID 

######################################################################################

KEGGIDtrans<- openxlsx::read.xlsx("metaprotein/KEGGIDtrans.xlsx") #genus ananotation file

split_taxonomy <- strsplit(KEGGIDtrans$Taxonomy, split = ";")

taxonomy_expanded <- do.call(rbind, lapply(split_taxonomy, function(x) {
  length(x) <- max(sapply(split_taxonomy, length))   
  return(x)
}))

taxonomy_expanded <- as.data.frame(taxonomy_expanded)

taxonomy_expanded$Accession <- KEGGIDtrans$Accession

KEGGIDtrans1 <- merge(KEGGIDtrans,taxonomy_expanded,by="Accession")

colnames(KEGGIDtrans1)[4:11] <- c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")

protein_g <- merge(KEGGIDtrans1[c(1,10)],protein2,by.x = "Accession",by.y = "PG.ProteinAccessions")

protein_g <- protein_g[,-1]

g_search <- unique(protein_g$Genus)

all_non_zero<- data.frame()

for(i in 1:length(g_search)){

g_pos <- g_search[i]

protein_gdu <- protein_g[which(protein_g$Genus%in%g_pos),]

non_zero_count <- sapply(protein_gdu[, -1], function(x) sum(x != "NaN"))

non_zero_count <- t(as.data.frame(non_zero_count))

rownames(non_zero_count) <- g_pos

all_non_zero<-rbind(all_non_zero, non_zero_count)
}

protein_g_num <-  all_non_zero

##########################################################################################################################

##########002 Metagenomic based on genus

##########################################################################################################################

genus.pt <- read.table("metagene/tax_g_RPKM.xls",header = T)

genus.pt$genus <-unlist( lapply(1:nrow(genus.pt), function(t) {
  a <- unlist(strsplit(genus.pt$Taxonomy[t], split = ";"))[7]
}))

genus.p.matrix <- genus.pt[,2:21]
rownames(genus.p.matrix) <- genus.pt$genus

##########################################################################################################################

##########003 correlation coefficient

a <- intersect(rownames(protein_g_num),genus.pt$genus )

###### 
nonGH.p <- grep("B_", colnames(protein_g_num))

GH.p<- grep("A_", colnames(protein_g_num))

protein_g_num <- protein_g_num[which(rownames(protein_g_num)%in%a),c(nonGH.p,GH.p)]

nonGH.p <- grep("B_", colnames(genus.p.matrix))
GH.p<- grep("A_", colnames(genus.p.matrix))

gene.abundance <- genus.p.matrix[which(rownames(genus.p.matrix)%in%a),colnames(protein_g_num)]

P.cor<- corr.test(t(protein_g_num),
                  t(protein_g_num),
                  method = "spearman",adjust="BH")

 

##############

cor.dt.all <- data.frame(genus="NA",
                     cor="NA",
                     pvalue="NA")

for(i in 1:length(a)){
  
g_postion <- rownames(protein_g_num)[i]

P.cor<- cor.test(as.matrix(gene.abundance[which(rownames(gene.abundance)%in%g_postion),]),
                 
                  as.matrix(protein_g_num[which(rownames(protein_g_num)%in%g_postion),]),
                 
                  method = "spearman"
                 )
cor.dt <- data.frame(genus="NA",
                         cor="NA",
                         pvalue="NA")
cor.dt$cor <- P.cor$estimate

cor.dt$pvalue <- P.cor$p.value

cor.dt$genus <- g_postion

cor.dt.all<-rbind(cor.dt.all, cor.dt)

}
View(cor.dt.all)

cor.dt.all <- na.omit(cor.dt.all)

cor.dt.all <- cor.dt.all[-grep("NA",cor.dt.all$genus),]

cor.dt.all <- cor.dt.all[-grep("g__unclassified",cor.dt.all$genus),]


cor.dt.all$adjp <- p.adjust(cor.dt.all$pvalue,"BH")

cor.cut <- cor.dt.all[cor.dt.all$adjp<0.01 &cor.dt.all$cor>0.7,]
 

row_sums <- rowSums(gene.abundance)
 
genus.p.matrix_with_sums <- cbind(gene.abundance, row_sums)

row_sums.p <- rowSums(protein_g_num)

protein_g_num_with_sums <- cbind(protein_g_num, row_sums.p)


genus.p.matrix_with_sums <- genus.p.matrix_with_sums[which(rownames(genus.p.matrix_with_sums)%in%cor.cut$genus),]

protein_g_num_with_sums <- protein_g_num_with_sums[which(rownames(protein_g_num_with_sums)%in%cor.cut$genus),]

genus.p.matrix_with_sums$genus <- rownames(genus.p.matrix_with_sums)

protein_g_num_with_sums$genus <- rownames(protein_g_num_with_sums)

df <- merge(genus.p.matrix_with_sums[,c("row_sums","genus" )],protein_g_num_with_sums[,c("row_sums.p","genus" )],by="genus")

df <- merge(df,KEGGIDtrans1[,c("Phylum","Genus")],by.x="genus",by.y = "Genus")
df1 <- df[!duplicated(df),]
table(df1$Phylum)
 

Colors <- c("#038FC4", "#4FAD9F", "#7FD165", "#A7CBCC", 
  "#D09978", "#F7979C", "#F6B75D", "#D9190B", "#FF8600",
  "#F19462", "#B893CB", "#713A9F"
)
ggplot(df1, aes(x = log10(row_sums.p), y = log10(row_sums))) + 
  geom_point(aes(color=Phylum),size = 4) +   
  geom_text(aes(label = genus), vjust = -0.5, family="serif",size = 3) +   
  scale_color_manual(values = Colors[1:length(names(table(df1$Phylum)))])+
  theme_bw()+
  labs(x = "Number of meta-proteins from each genus (log10)", 
       y = "Genus abundance based on metagenomics (log10)", title = "R >0.7, adj P <0.05") +
  theme(text = element_text(family="serif",size = 8))  #  
################################################################################
 
#############02 calculate alpha index and visualization

################################################################################
 

anno_overview <- read.table("metagene/anno_overview.csv",sep = ",", header = TRUE, fileEncoding = "utf-8")

anno_overview1 <- anno_overview[,c(1,3:10)]

rownames(anno_overview1) <- anno_overview1$GeneID

anno_overview1 <- anno_overview1[,-1]

#######

data <- read.table("metagene/reads_number.xls",header = T)

rownames(data) <- data$GeneID
data1 <- data[,-c(1,22)]

data1 <- data1%>%
  mutate(across(everything(), as.numeric))

nonGH.p <- grep("B_", colnames(data1))
GH.p<- grep("A_", colnames(data1))

data1 <- data1[,c(GH.p,nonGH.p)]

data1 <- t(data1)

#######samplegroup
samplegroup <- data.frame(Sample=rownames(data1),
                          Group=c(rep("GHPA",10),rep("Ctrl",10)))


mpse1 <- mp_import_dada2(seqtab=data1, taxatab=as.matrix(anno_overview1), 
                         sampleda=samplegroup %>% 
                           tibble::column_to_rownames("Sample"))

mpse1 %<>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance,
    chunks = 200
  )

mpse1 %>% print(width=180)


############

mpse1 %<>% 
  mp_cal_alpha(.abundance=RareAbundance,force=TRUE)

 mpse1 %>% 
  mp_plot_alpha(
    .group=Group, 
    .alpha=c(Observe, Chao1, ACE, Shannon, Simpson)##Pielou
  ) +
  scale_fill_manual(values=c("#5E7BBB","#EF9703"), guide="none") +
  scale_color_manual(values=c("#5E7BBB","#EF9703"), guide="none")+
  theme(text = element_text(family="serif"))
 
#openxlsx::write.xlsx(alpha_re,"metagene/alpha_re.xlsx")

    
####################################################################################
 
alpha_re <- alpha_re[,c("Group","Observe","Chao1","ACE","Shannon","Simpson")]

data <- reshape2::melt(
  data = alpha_re,
  id.vars = c("Group", "Observe"),  
  measure.vars = c("Chao1", "ACE", "Shannon", "Simpson"),  
  variable.name = "Name",    
  value.name = "Value"      
)

sig<-compare_means(Value~Group, 
                   data=data,
                   method = "wilcox.test",group.by = "Name", 
                   p.adjust.method = "BH") #t.test

ggplot(data=data, aes(x=Group,y=Value))+
  geom_boxplot(aes(color=Group),size=0.8, outliers=F, position=position_dodge(1))+
  geom_jitter(alpha=0.5,aes(color=Group),
              position=position_jitterdodge(jitter.width = 0.3, 
                                            jitter.height = 0, 
                                            dodge.width = 0.4
              ))+ 
  scale_color_manual(values = c("#5E7BBB","#EF9703"))+
  stat_compare_means(aes(group = Group),
                     method = "wilcox.test",
                     label = "p.signif",
                     label.x.npc = "centre")+
  facet_wrap(~ Name, scales="free")+
  ylab("Value")+ #
  theme_bw()+
  theme(text = element_text(family="serif",size = 8),#"Times New Roman","serif"
        panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 90),
        axis.text.y = element_text(size=8),
        legend.position = "none"
  )

 

################################################################################

#############03 Beta diversity analysis

################################################################################

#################

mpse1 %<>% 
  mp_decostand(.abundance=Abundance)

mpse1 %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")

p <- mpse1 %>% mp_plot_dist(.distmethod = bray)

p3<- mpse1 %>% mp_plot_dist(.distmethod = bray, .group = Group, 
                            group.test=TRUE, 
                            test = "wilcox.test",
                            textsize=2) +
  scale_color_manual(values=c("#EF9703","#3C5488FF","#5E7BBB")) +
  scale_fill_manual(values=c("#EF9703","#3C5488FF","#5E7BBB"),guide="none")+
  theme(text = element_text(family="serif"))  
p3 
##############


################################################################################

#############04 Linear discriminant analysis

################################################################################

##########
#####LDA
 
Gene <- read.table("metagene/tax_g_RPKM.xls",header = T)

taxonomy <- Gene$Taxonomy

Gene_filter <- Gene[,1:21]
rownames(Gene_filter) <- Gene_filter$Taxonomy

Gene_filter <- Gene_filter[,-1]

nonGH.p <- grep("B_", colnames(Gene_filter))
GH.p<- grep("A_", colnames(Gene_filter))

Gene_filter1 <- Gene_filter[,c(GH.p,nonGH.p)]

Gene_filter2 <- Gene_filter1 %>%
  mutate(across(everything(), as.numeric))

rownames(Gene_filter2) <- rownames(Gene_filter1)

###############filtering 


##############taxa file

split_taxonomy <- strsplit(taxonomy, split = ";")

 
taxonomy_expanded <- do.call(rbind, lapply(split_taxonomy, function(x) {
  length(x) <- max(sapply(split_taxonomy, length))   
  return(x)
}))

colnames(taxonomy_expanded) <- c("Domain","Kingdom","Phylum","Class","Order","Family","Genus")

rownames(taxonomy_expanded) <- taxonomy

sample_table <-data.frame(sample_ID=colnames(Gene_filter2),
                          Group=c(rep("GH",10),rep("Ctrl",10)))

rownames(sample_table) <- sample_table$sample_ID

tax_table <-data.frame("ID"=taxonomy_expanded[,7],
                       "genus"=taxonomy_expanded[,7])

rownames(Gene_filter2) <- taxonomy_expanded[,7]

dataset <- microtable$new(sample_table = sample_table,
                          otu_table = as.data.frame(Gene_filter2), 
                          tax_table = tax_table)

#The feature abundance table; rownames are features (e.g. OTUs/ASVs/species/genes); column names are samples.
 
lefse <- trans_diff$new(dataset = dataset,  
                        method = "lefse",  
                        group = "Group",  
                        alpha = 0.05, # significance threshold to select taxa 
                        add_sig_label = "Significance",
                        taxa_level = "genus", # phylum
                        lefse_min_subsam = 8,
                        p_adjust_method = "BH"
) 

 
View(lefse$res_diff)

index <- which(lefse$res_diff$LDA > 2.5 &lefse$res_diff$P.adj<0.05)

cut.les <- as.data.frame(lefse$res_diff)[index,]
 
cut.les$Group <- factor(cut.les$Group,levels = c("GHPA", "Ctrl"))

cut.les<- arrange(cut.les,desc(P.adj))

b <- grep("g__unclassified", cut.les$Taxa)

cut.les <- cut.les[-b,]

cut.les$Taxa <- sapply(str_split(cut.les$Taxa,"\\|"), "[",1)

cut.les$Taxa <- sapply(str_split(cut.les$Taxa,"g__"), "[",2)

cut.les <- cut.les %>%  
  arrange(Group,P.adj)

cut.les$Taxa <- factor(cut.les$Taxa,levels = rev(cut.les$Taxa),ordered = T)

 ggplot(cut.les, aes(x = LDA,y = Taxa,fill = -log10(P.adj))) +
  geom_col(color = "white",width=0.85,linewidth=0.7) +
  labs(y="",fill="-log10 P.adj")+
  scale_color_manual(values =  "#D95201")+
  scale_fill_gradient(low = "#ffffff",high = "#D95201",
                      limits = c(0,4),breaks = c(0,1,2,3,4))+
  scale_x_continuous(breaks = seq(0, 5, by = 2),
                     limits = c(0, 5),expand = c(0,0)) +
  theme_classic() +
  theme(axis.text.y = element_text(family = "serif",size = 9),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.border = element_blank(),
        legend.title.position = "left",
        legend.key.height=unit(0.5, "cm"),
        legend.key.width=unit(0.4, "cm"),
        legend.title = element_text(family = "serif",size = 9,hjust = 0.5, angle = 90),
        axis.line = element_line(linewidth = 0.5,lineend = "square"),
        plot.margin=unit(x=c(0,0.2,0.1,0.2),units="inches"))

##################################
 
###################################
 
 rownames(Gene_filter2) <- sapply(str_split( rownames(Gene_filter2),"g__"), "[",2)
 
 cut_abundance <- Gene_filter2[which(rownames(Gene_filter2)%in%cut.les$Taxa),]
 
 cut_abundance1 =as.data.frame(t(t(cut_abundance)/colSums(cut_abundance,na=T)*100))

 #####################
 Bilophila.plot <- as.data.frame(t(cut_abundance[which(rownames(cut_abundance)%in%"Bilophila"),]))
 
 group <- c(rep("GHPA",10),rep("Ctrl",10))
 
 Bilophila.plot$Group <- factor(group, levels=c("Ctrl","GHPA"))
 
 sig<-compare_means(Bilophila~Group, 
                    data=Bilophila.plot,
                    method = "wilcox.test",
                    p.adjust.method = "BH",) #t.test
 
 y.max <-max(Bilophila.plot$Bilophila)+1.3
 
 lab.max <- y.max-0.2
 
  ggplot(data=Bilophila.plot, 
                     aes(x=Group,y=Bilophila))+
   annotate("text", x=2, 
            y=lab.max, 
            label= sig$p.signif)+
    geom_boxplot(aes(color=Group),
                 size=0.8, 
                 outliers=F, 
                 position=position_dodge(1))+
   geom_jitter(alpha=0.5,aes(color=Group),
               position=position_jitterdodge(jitter.width = 0.3, 
                                             jitter.height = 0, 
                                             dodge.width = 0.4
               ))+
   scale_color_manual(values = c("#5E7BBB","#EF9703"))+
   theme_bw()+
   ylab("Abundance levels")+
   labs(title="Bilophila")+
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
 
 
 
 