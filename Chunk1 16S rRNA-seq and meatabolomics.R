####################################################################################################################################################################
#      ☆        % Project: Gut Microbiota Alterations in Acromegaly Patients Are Associated with Neutrophil Depletion-Induced Inflammation #
#   ☆ \|/ ☆    % Author: HuaChun Yin                                         
#  ☆  \|/  ☆   % Date: Apr. 4th, 2024                                
# ☆   \|/   ☆  %                                                          
#  ☆  \|/  ☆   % Environment: R version 4.4.2           
#  ☆ __|__ ☆   % EPlatform: Mac-IOS(64-bit)                                  
#                % CHUNK1:this script includes differential metabolites identification, rare ASV removeal                        
################################################################################################################################################################### 


packages <- c("dplyr","doParallel","Biostrings",
              "ggdendro",
              "ggplot2","ggpubr","ggsci","ggrepel",
              "psych","plotly","FactoMineR",
              "stringr","corrplot",
              "RColorBrewer",
              "ggsci","microeco",
              "tidyverse","magrittr",
              "data.table","foreach",
              "lme4","nlme","Rtsne",
              "factoextra","vegan",
              "ggalluvial","PMA",
              "ggrepel","scales","hdi","ropls", "glmnet","data.table","hdi","stabs",
              "stabs","Rtsne","ropls","MicrobiotaProcess","philentropy","Hmisc","cowplot") ##  

#check.packages(packages)

lapply(packages, library, character.only = TRUE)

################################################################################

#############01 Calculate alpha index and visualization

################################################################################

ASV <- read.table("Data/ASV_table_tax.xls",header = T)

taxonomy <- ASV[,c("ASV_ID","taxonomy")]

ASV_filter <- ASV[,1:61]

rownames(ASV_filter) <- ASV_filter$ASV_ID

rownames(taxonomy) <- taxonomy$ASV_ID
 
taxonomy_filter <- taxonomy[rownames(ASV_filter),]

taxonomy_filter$phylum <-unlist( lapply(1:nrow(taxonomy_filter), function(t) {
  a <- unlist(strsplit(taxonomy_filter$taxonomy[t],"\\|"))[2]
  unlist(strsplit(a,"_"))[3]
}))

taxonomy_filter$genus <- unlist(lapply(1:nrow(taxonomy_filter), function(t) {
  a <- unlist(strsplit(taxonomy_filter$taxonomy[t],"\\|"))[6]
 
}))

taxonomy_filter <- as.matrix(taxonomy_filter)

names(table(na.omit(as.character(taxonomy_filter[,2])))) #200 microbial taxa

ASV_filter1 <- ASV_filter[,-1]

rownames(ASV) <- ASV$ASV_ID


################

Species_relative <- openxlsx::read.xlsx("Data/Species_relative.xlsx",rowNames = TRUE)

Species_relative <- t(Species_relative)

Species_relative <- as.data.frame(Species_relative)

Species_relative <- Species_relative %>%
  rowwise() %>%
  mutate(sum = sum(c_across(everything()) > 0))

Species_sum <- Species_relative$sum

data <- data.frame(Species=Species_sum,
                   Group=c(rep("GHPA",35),rep("Ctrl",25)))

sig<-compare_means(Species~Group, 
                   data=data,
                   method = "wilcox.test",
                   p.adjust.method = "BH") #t.test

  ggplot(data=data, aes(x=Group,y=Species))+
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
  ylab("Number of observed species")+ #
  theme_bw()+
  theme(text = element_text(family="serif",size = 8),#"Times New Roman","serif"
        panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 90),
        axis.text.y = element_text(size=8),
        legend.position = "none"
  )

 
##########
 
# input FASTA file

fasta_file <- "Data/ASV_seq.fasta"  

fasta_data <- readDNAStringSet(fasta_file)

 
seq_names <- names(fasta_data)

sequences <- as.character(fasta_data)

seq_matrix <- cbind(seq_names, sequences)

ASVseq_matrix <- as.data.frame(seq_matrix, stringsAsFactors = FALSE)


##############taxa file

View(taxonomy)

taxonomy <- merge(taxonomy,ASVseq_matrix,by.x="ASV_ID",by.y="seq_names" )

rownames(taxonomy) <- taxonomy$sequences

split_taxonomy <- strsplit(taxonomy$taxonomy, split = "\\|")

 
taxonomy_expanded <- do.call(rbind, lapply(split_taxonomy, function(x) {
  length(x) <- max(sapply(split_taxonomy, length))   
  return(x)
}))

rownames(taxonomy_expanded) <- rownames(taxonomy)

colnames(taxonomy_expanded) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
#############

#########ASV sequence

ASV <- ASV[,-62]

ASVseq <- merge(taxonomy[,c(1,3)],ASV,by="ASV_ID")

rownames(ASVseq) <- ASVseq$sequences

ASVseq <- ASVseq[,-c(1:2)]

ASVseq <- ASVseq%>%
  mutate(across(everything(), as.numeric))
ASVseq <- t(ASVseq)

#######samplegroup
samplegroup <- data.frame(Sample=rownames(ASVseq),
                          Group=c(rep("GHPA",35),rep("Ctrl",25)))

#write.csv(samplegroup,"Data/samplegroup.csv")

samplegroup <- "Data/samplegroup.csv"

mpse1 <- mp_import_dada2(seqtab=ASVseq, 
                         taxatab=taxonomy_expanded, 
                         sampleda=samplegroup)

View(mpse1)

mpse1 %<>% mp_rrarefy()# 

mpse1 %<>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance,
    chunks = 400
  )

mpse1 %>% print(width=180)

###################
p1 <- mpse1 %>% 
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve, 
    .alpha = Observe,
  )+
  theme(text = element_text(family="serif"))

p2 <- mpse1 %>% 
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve, 
    .alpha = Observe, 
    .group = Group
  ) +
  scale_color_manual(values=c("#5E7BBB","#EF9703")) +
  scale_fill_manual(values=c("#5E7BBB","#EF9703"), guide="none")+
  theme(text = element_text(family="serif"))
# combine the samples belong to the same groups if 
# plot.group=TRUE
p3 <- mpse1 %>% 
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve, 
    .alpha = "Observe", 
    .group = Group, 
    plot.group = TRUE
  ) +
  scale_color_manual(values=c("#5E7BBB","#EF9703")) +
  scale_fill_manual(values=c("#5E7BBB","#EF9703"),guide="none")+
  theme(text = element_text(family="serif"))  

p1 + p2 + p3
#save as readsNums1-3
############

mpse1 %<>% 
  mp_cal_alpha(.abundance=RareAbundance)

f1 <- mpse1 %>% 
  mp_plot_alpha(
    .group=Group, 
    .alpha=c(Observe, Chao1, ACE, Shannon, Simpson)#Pielou
  ) +
  scale_fill_manual(values=c("#5E7BBB","#EF9703"), guide="none") +
  scale_color_manual(values=c("#5E7BBB","#EF9703"), guide="none")+
  theme(text = element_text(family="serif"))

f2 <- mpse1 %>%
  mp_plot_alpha(
    .alpha=c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
  )+
  theme(text = element_text(family="serif"))

alpha_re <- mp_extract_sample(mpse1)

#openxlsx::write.xlsx(alpha_re,"alpha_re.xlsx")

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


############
p1 <- mpse1 %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Group, 
    taxa.class = Phylum, 
    topn = 20,
    relative = TRUE
  )
# visualize the abundance (rarefied) of top 20 phyla for each sample.
p2 <- mpse1 %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Group,
    taxa.class = Phylum,
    topn = 20,
    relative = FALSE
  )

p3 <- mpse1 %>%
  mp_plot_abundance(
    .abundance=RareAbundance, 
    .group=Group,
    taxa.class = Phylum,
    topn = 20,
    plot.group = TRUE
  )

# visualize the abundance of top 20 phyla for each .group (time)
p4 <- mpse1 %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group= Group,
    taxa.class = Phylum,
    topn = 20,
    relative = FALSE,
    plot.group = TRUE
  )
p3 / p4

################################################################################

#############02 Beta diversity analysis

################################################################################

#################
mpse1 %<>% 
  mp_decostand(.abundance=Abundance)

mpse1 %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")
 
p5<- mpse1 %>% mp_plot_dist(.distmethod = bray, .group = Group, 
                            group.test=TRUE, 
                            test = "wilcox.test",
                            textsize=2) +
  scale_color_manual(values=c("#5E7BBB","#3C5488FF",,"#EF9703")) +
  scale_fill_manual(values=c("#5E7BBB","#3C5488FF",,"#EF9703"),guide="none")+
  theme(text = element_text(family="serif"))  
 
###########################################################################

###########03 rare ASV removeal

###########################################################################

###########################################################################
# read and filtering the ASV
###########################################################################

A=t(ASV_filter1)  #row=Sample col=OTU

C=A/rowSums(A)

genus=t(C)

 
genus <- as.data.frame(genus)

genus$ASV_ID <-rownames(genus)

genus <- merge(taxonomy_filter[,c("ASV_ID","genus")], genus)

genus <- aggregate(.~genus$genus,genus[,3:ncol(genus)], sum)

rownames(genus) <- genus$`genus$genus`

genus <- genus[,-1]

#microbial taxa had to be present in an sum relative abundance of greater than 0.01% across those samples.

genus <- genus[which(rowSums(genus) >= 0.0001), ] 

genus1 <- genus
 
genus1[genus1>0] <- 1

#microbial taxa had to be present in at least 5% of samples

genus1 <- genus[which(rowSums(genus1) >=3), ]  

# DOI:10.1186/s42523-023-00278-0

dim(genus1)

 
genus2 <- merge(taxonomy_filter[,c("ASV_ID","genus")], ASV_filter)

genus3 <- aggregate(.~genus2$genus,genus2[,3:ncol(genus2)], sum)

genus4 <- genus3[which(genus3$genus%in%rownames(genus1)),]

rownames(genus4) <- genus4$`genus2$genus`

genus4 <- genus4[,-1]# 
 
#######

genus5 <- genus4

genus5$rowsum <- apply(genus5,1,sum)

genus5 <- genus5[order (genus5$rowsum,decreasing=TRUE),]

genus5 = genus5[,-61]
 
genus6 <- apply(genus5,2,function(x) x/sum(x)) 

genus7 <-  genus6[1:50,]

genus8 <- 1-apply(genus7, 2, sum) 

#合并数据
genus9 <- rbind(genus7,genus8)

row.names(genus9)[51]="Others"

genus9 <- as.data.frame(genus9)

genus9$genus <- rownames(genus9)

##################################################
tax <- taxonomy_filter[,c("phylum","genus")]

tax <- tax[!duplicated(tax[,"genus"]),]

phylum <- genus4

phylum$genus <- rownames(genus4)

tax <- tax[which(tax[,"genus"]%in%phylum$genus),]

phylum.df <- merge(tax[,c("phylum","genus")], phylum)

phylum.df <- phylum.df[which(phylum.df$genus%in%genus9$genus),]

phylum.df1 <- aggregate(.~phylum.df$phylum,phylum.df[,3:ncol(phylum.df)], sum)

rownames(phylum.df1) <- phylum.df1$`phylum.df$phylum`

phylum.df1 <- phylum.df1%>%select(-`phylum.df$phylum`)

 
A=t(phylum.df1)  

C=A/rowSums(A)

phylum.top <- t(C)
 
df <- reshape2::melt(as.matrix(phylum.top))

colnames(df) <- c("phylum","Sample_ID","Proportion")

rhg_cols1 <-colorRampPalette(brewer.pal(8, "Set2"))(length(unique(df$phylum)))

P.phylum <- ggplot(df,aes(x=Sample_ID,y=Proportion*100,fill= phylum))+
  geom_bar(stat="identity")+
  theme_classic()+
  scale_fill_manual(values = rhg_cols1)+
  labs(fill = "phylum")+
  ylab("Relative abundance (%)")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1,family="serif"),
        legend.text = element_text(size = 6,family="serif"))+
  guides(color = guide_legend(override.aes = list(size = 4))) 

print(P.phylum)


####################

genus9 <- genus9%>%select(-"genus")
 
df <- reshape2::melt(as.matrix(genus9))

colnames(df) <- c("genus","Sample_ID","Proportion")

rhg_cols1 <-colorRampPalette(brewer.pal(8, "Set2"))(length(unique(df$genus)))

P.genus <- ggplot(df,aes(x=Sample_ID,y=Proportion*100,fill= genus))+
  geom_bar(stat="identity")+
  theme_classic()+
  scale_fill_manual(values = rhg_cols1)+
  labs(fill = "genus")+
  ylab("Relative abundance (%)")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1,family="serif"),
        legend.text = element_text(size = 6,family="serif"))+
  guides(color = guide_legend(override.aes = list(size = 4))) 

print(P.genus)

################################################################################

#############04 Linear discriminant analysis

################################################################################

##########
#####LDA
 
group <- substr(colnames(genus4),1,1)
group <- str_replace_all(group,"A","GHPA")
group <- str_replace_all(group,"B","Ctrl")

feature_table <- as.data.frame(genus4)

sample_table <-data.frame(sample_ID=colnames(genus4),
                          Group=group)
rownames(sample_table) <- sample_table$sample_ID

tax_table <-data.frame("ID"=rownames(genus4),
                       "genus"=rownames(genus4))

rownames(tax_table) <- tax_table$ID

dataset <- microtable$new(sample_table = sample_table,
                          otu_table = feature_table, 
                          tax_table = tax_table)
#The feature abundance table; rownames are features (e.g. OTUs/ASVs/species/genes); column names are samples.
 
#Differential Abundance Analysis）
lefse <- trans_diff$new(dataset = dataset,   
                        method = "lefse",  
                        group = "Group",   
                        alpha = 0.05, # significance threshold to select taxa 
                        add_sig_label = "Significance",
                        taxa_level = "genus", # 
                        lefse_min_subsam = 8,
                        p_adjust_method = "none"
) 
 
View(lefse$res_diff)

rownames(lefse$res_diff) <- sapply(str_split(rownames(lefse$res_diff),"g__"), "[",3)

lefse$res_diff$Taxa <- sapply(str_split(lefse$res_diff$Taxa,"g__"), "[",3)

index <- which(lefse$res_diff$LDA > 2.5 &lefse$res_diff$P.adj<0.05)

cut.les <- as.data.frame(lefse$res_diff)[index,]

#openxlsx::write.xlsx(cut.les,file="lefse_result_cut.xlsx")
#openxlsx::write.xlsx(lefse$res_diff,file="lefse_result_nocut.xlsx")
 
cut.les$P.adj <- p.adjust(cut.les$P.adj,method = "BH")

cut.les$Group <- factor(cut.les$Group,levels = c("GHPA", "Ctrl"))

cut.les<- arrange(cut.les,desc(P.adj))

cut.les <- cut.les %>%  
  arrange(Group,P.adj)

cut.les$Taxa <- factor(cut.les$Taxa,levels = rev(cut.les$Taxa),ordered = T)

cut.lesGH <- cut.les[which(cut.les$Group%in%"GHPA"),] 

cut.lesnon <- cut.les[which(cut.les$Group%in%"Ctrl"),]

mytheme1 <- theme_classic() +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(family = "serif",size = 9),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        legend.title.position = "left",
        legend.key.height=unit(0.5, "cm"),
        legend.key.width=unit(0.4, "cm"),
        legend.title = element_text(family = "serif",size = 9,hjust = 0.5, angle = 90),
        plot.margin=unit(x=c(0.5,0.2,0,0.2),units="inches"))
 

f1 <- ggplot(cut.lesGH, aes(x = LDA,y = Taxa,fill = -log10(P.adj))) +
  geom_col(color = "white",width=0.85,linewidth=0.7) +
  labs(y="",fill="-log10 P.adj")+
  scale_color_manual(values =  "#EF9703")+
  scale_fill_gradient(low = "#ffffff",high = "#FA9405",
                      limits = c(0,3),breaks = c(0,1,2,3))+
  scale_x_continuous(breaks = seq(0, 5, by = 2),
                     limits = c(0, 5),expand = c(0,0)) +
  
  mytheme1

mytheme2 <- theme_classic() +
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

f2 <- ggplot(cut.lesnon, aes(x = LDA,y = Taxa,fill = -log10(P.adj))) +
  geom_col(color = "white",width=0.15,linewidth=0.1) +
  labs(x="LDA score",y="",fill="-log10 P.adj")+
  scale_color_manual(values =  "#5E7BBB")+
  scale_fill_gradient(low = "#ffffff",high = "#5E7BBB",
                      limits = c(0,3),breaks = c(0,1,2,3))+
  scale_x_continuous(breaks = seq(0, 5, by = 2),
                     limits = c(0, 5),expand = c(0,0)) +
  scale_y_discrete(expand = expansion(mult = c(0.1, 0.1)))+
  
  mytheme2
 
plot_grid(f1, f2, ncol = 1, align = "v")

###################

taxonomy_expanded1 <- taxonomy_expanded[which(taxonomy_expanded[,6]%in%rownames(genus4)),]

taxonomy_expanded2 <- taxonomy_expanded1 %>% 
  as.data.frame()%>% 
  distinct(Genus, .keep_all = TRUE)

cut.les1 <- cut.les

cut.les1$Taxa <- paste0("g__",cut.les1$Taxa)

cut.les1 <- cut.les1[c(2,4,7)]

colnames(cut.les1) <- c("f","DIAGNOSIS","pvalue")

cut.taxonomy <- taxonomy_expanded2[which(taxonomy_expanded2$Genus%in%cut.les1$f),]

Phylum <- data.frame(
  f=unique(cut.taxonomy$Phylum),
  DIAGNOSIS=rep("GHPA",length(unique(cut.taxonomy$Phylum))),
  pvalue =rep(0.00001,length(unique(cut.taxonomy$Phylum)))
)

Class <- data.frame(
  f=unique(cut.taxonomy$Class),
  DIAGNOSIS=rep("GHPA",length(unique(cut.taxonomy$Class))),
  pvalue =rep(0.0001,length(unique(cut.taxonomy$Class)))
)
Order <- data.frame(
  f=unique(cut.taxonomy$Order),
  DIAGNOSIS=rep("GHPA",length(unique(cut.taxonomy$Order))),
  pvalue =rep(0.001,length(unique(cut.taxonomy$Order)))
)

Family <- data.frame(
  f=unique(cut.taxonomy$Family),
  DIAGNOSIS=rep("GHPA",length(unique(cut.taxonomy$Family))),
  pvalue =rep(0.01,length(unique(cut.taxonomy$Family)))
)


dt <- rbind(cut.les1,Phylum,Class,Order,Family)
#data(df_alltax_info) demo data
#data(df_difftax)

taxa <- taxonomy_expanded2

dt <- cut.les1

ggdiffclade(obj=taxa,
            nodedf=dt,
            factorName="DIAGNOSIS",
            layout="radial",
            skpointsize=0.9,
            cladetext=2,
            linewd=0.3,
            taxlevel=5,
            # This argument is to remove the branch of unknown taxonomy.
            reduce=TRUE) + 
  scale_fill_manual(values=c("#5E7BBB","#EF9703"))+
  guides(color = guide_legend(keywidth = 0.1, keyheight = 0.6,
                              order = 2,ncol=1)) +
  theme(panel.background=element_rect(fill=NA),
        legend.position="right",
        plot.margin=margin(0,0,0,0),
        legend.spacing.y=unit(0.02, "cm"), 
        legend.title=element_text(size=7.5), 
        legend.text=element_text(size=5.5), 
        legend.box.spacing=unit(0.02,"cm")
  )

###########################################################################

###########05 differential metabolites between GHPA and Ctrl

###########################################################################

###############
POS <- openxlsx::read.xlsx("Data/file1_resultPosAndNeg.xlsx",sheet=1)
NEG <- openxlsx::read.xlsx("Data/file1_resultPosAndNeg.xlsx",sheet=2)

Metabolites.m <- rbind(POS,NEG)

delet.col <- c( "m/z", "rt(s)","Name","adduct", 
                "score","HMDB","KEGG","SuperClass","Class","SubClass",
                "QC-1","QC-2","QC-3","QC-4","QC-5","QC-6" ,"QC-7" )
#########

#########
Metabolites.m <- Metabolites.m%>%select(-delet.col)

group <- substr(colnames(Metabolites.m[,-1]),1,1)

group <- str_replace_all(group,"A","GHPA")

group <- str_replace_all(group,"B","Ctrl")

##########

initial_value<-40

theta_value <- 0.9 #[0,1]

plex_value <- 15 # < 1/3 samples


#######

set.seed(5)

tsne_out<-Rtsne(log2(t(Metabolites.m[,-1])+1),
                dims = 2,initial_dims = initial_value,
                pca = FALSE,
                perplexity = plex_value,
                theta = theta_value, 
                max_iter = 1000
)

tsne <- data.frame(tSNE1=tsne_out$Y[,1],
                   tSNE2=tsne_out$Y[,2],
                   Group=group,
                   ID=colnames(Metabolites.m[,-1]))
 
#############
tsne_df <- data.frame(tsne_out$Y, group = group)

colnames(tsne_df) <- c("Dim1", "Dim2", "group")

# Calculate p-value using t-test for one of the dimensions  

group1 <- tsne_df %>% filter(group == "GHPA") %>% pull(Dim1)

group2 <- tsne_df %>% filter(group == "Ctrl") %>% pull(Dim1)

# Perform t-test
t_test_result <- t.test(group1, group2)

# Extract p-value
p_value <- t_test_result$p.value

p_value
#0.0003985451
##########
ggplot(tsne,aes(tSNE1,tSNE2))+
  geom_point(aes(color=group),size=1.5)+
  scale_color_lancet()+
  scale_color_manual(values=c("#5E7BBB","#EF9703"))+
  geom_text_repel(aes(tSNE1,tSNE2, label=ID),size=4, family="Times",arrow = arrow(length=unit(0.01, "npc")),
                 force = 1, max.iter = 3e3)+
  theme_bw()+
  xlim(-100,100)+#300
  ylim(-150,150)+#150
  theme(plot.margin = unit(rep(1.5,4),"lines"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.9,0.2), # legend???????????????
        legend.background = element_rect(size = 1, colour = "white"))

tsne_out$Y[1,]

 
############################

POS.cleaned <- POS[complete.cases(POS$Name), ]

NEG.cleaned <- NEG[complete.cases(NEG$Name), ]

vip.df <- c("POS.cleaned","NEG.cleaned")

for(ion in 1:2){
df <- get(vip.df[ion])

rownames(df) <- df$ID

df <- df[,12:71]

group <- substr(colnames(df),1,1)

group <- str_replace_all(group,"A","GHPA")

group <- str_replace_all(group,"B","Ctrl")


set.seed(1)

oplsda = opls(log2(t(df)), group,predI = 1, orthoI = NA, crossvalI=7)

vip <- getVipVn(oplsda)

vip <- as.data.frame(vip)

vip$FeatureName<-rownames(vip)

#
all.metabolits <- rownames(df)

input <-as.data.frame(t(df))

gene.count <- length(all.metabolits)

group.factor <- as.factor(group)

deg.count <- NULL

i <- 0

pb <- txtProgressBar(min = 0, max = gene.count, style = 3, char = "+")

for (g in all.metabolits) {
  
  i <- i + 1
  
  setTxtProgressBar(pb, i)
  
  # Display the progress bar!
  
  t.testp <- try(t.test(get(g) ~ group.factor, 
                        data = input), 
                 silent = FALSE)
  
  
  deg.p <- t.testp$p.value
  
  flush.console()
  
  # This variable, deg.count, stores all the differently expressed genes.
  
  if (!is.na(deg.p) & deg.p < 1) deg.count <- c(deg.count, g,deg.p) else next
  
}

close(pb)

deg.count1<-matrix(deg.count,ncol=2,byrow = T)

colnames(deg.count1)<-c("FeatureName","Pvalue")

deg.count1<-as.data.frame(deg.count1)

deg.count1$Pvalue<-as.vector(deg.count1$Pvalue)

deg.count1$Pvalue<-as.numeric(deg.count1$Pvalue)

mm1 <- merge(vip,deg.count1,by="FeatureName")

fdrp.pos=p.adjust(mm1$Pvalue, "BH")

mm1$AdjP <- fdrp.pos

print("pvalue")
#########FC#

e.p <- which(group%in%"GHPA")  

c.p <- which(group%in%"Ctrl")  

FC <-apply(df, 1, function(x) {mean(x[e.p])/mean(x[c.p])})

FC <-as.data.frame(FC)

FC$"log2(FC)"<-log2(FC$FC)

FC$FeatureName<-rownames(FC)

print("FC")

mm1 <- merge(mm1,FC,by="FeatureName")

assign(paste0(vip.df[ion],"DM"),mm1)

}

POSDM_fileter <- POS.cleanedDM[abs(POS.cleanedDM$`log2(FC)`)>1 &POS.cleanedDM$Pvalue<0.05&POS.cleanedDM$AdjP<0.1 & POS.cleanedDM$vip>1, ]

POSDM_fileter <- merge(POSDM_fileter,POS,by.x = "FeatureName",by.y="ID")

NEGDM_fileter <- NEG.cleanedDM[abs(NEG.cleanedDM$`log2(FC)`)>1 &NEG.cleanedDM$AdjP<0.1& NEG.cleanedDM$Pvalue<0.05& NEG.cleanedDM$vip>1, ]

NEGDM_fileter <- merge(NEGDM_fileter,NEG,by.x = "FeatureName",by.y="ID")

DM_fileter <- rbind(POSDM_fileter,NEGDM_fileter)
############

######## 

DM <- DM_fileter %>% 
  mutate(
    diff=case_when(
      log2(FC) > 0 ~ "up",
      
      log2(FC) <0 ~ "down")
  )


DM_cut <- DM[,c("Name","vip","Pvalue","AdjP","log2(FC)","diff","adduct",    
                "HMDB","KEGG","SuperClass","Class","SubClass" )]


colnames(DM_cut) <- c("Name","VIP","Pvalue","AdjP","log2FC","diff","adduct",    
                      "HMDB","KEGG","SuperClass","Class","SubClass" )



ggplot(DM_cut, aes(x = "Pvalue", y = "VIP", z = "log2FC"))+geom_point(alpha=.5)


P <- plot_ly(DM_cut, 
             x = ~Pvalue, 
             y = ~VIP, 
             z = ~log2FC, 
             color = ~DM_cut$diff,
             text=~DM_cut$Name,
             mode="markers",
             colors = c("#5E7BBB","#EF9703"),
             marker = list(size = 5))%>% 
  add_markers(alpha=0.6)%>%
  layout(annotations = list(x = DM_cut$Pvalue,
                            y = DM_cut$VIP,
                            z = DM_cut$log2FC,
                            text =DM_cut$Name,
                            showarrow = FALSE,
                            xanchor = "center",
                            yanchor = "bottom",
                            font = list(size = 6)))

#openxlsx::write.xlsx(DM_cut,file="DM.xlsx")

DM_cut_abandence <-DM[,c(9,17:76)]

rownames(DM_cut_abandence) <- DM_cut_abandence$Name

DM_cut_abandence<- DM_cut_abandence%>% select(-"Name")

colnames(DM_cut_abandence) #check the colnames

#openxlsx::write.xlsx(DM_cut_abandence,file="DM_cut_abandence.xlsx",rowNames = TRUE)
#####################

anno<- DM[,c("SuperClass","Class")]

anno$SuperClass <- as.factor(anno$SuperClass)

anno$Class <- as.factor(anno$Class)

anno<-as.data.frame(anno)

rownames(anno)<-rownames(DM_cut_abandence)

annocol<-as.data.frame(group)###sam.lab1

rownames(annocol)<-colnames(DM_cut_abandence)

colnames(annocol)<-"Group"

col_colors<-get_anno_for_heatmap2(anno,annocol)

col_colors$Group <- c("GHPA"="#EF9703","Ctrl"="#5E7BBB")
 

a<-hclust(dist_n(t(DM_cut_abandence),mtd="manhattan"))  

df2<-dendro_data(a,type="rectangle") #"rectangle", "triangle"

tree.manhattan <- df2$labels$label

Group<- substr(df2$labels$label,1,1)

table(Group)

ggplot(segment(df2))+
  geom_text(data=df2$labels,aes(x=x,y=y,label=label,color=Group),
            angle=90,hjust=1,vjust=0.3,size=3)+
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend))+
  theme(panel.background = element_rect(fill = "white"))

p.heatmap <- pheatmap(as.matrix(log2(DM_cut_abandence[,as.character(tree.manhattan)])),
                      show_rownames = T,
                      show_colnames = T,
                      cluster_cols = F,
                      cluster_rows=T,
                      show_row_dend=F,
                      fontsize_row=5, # 
                      border_color = "NA", 
                      scale = "row",  
                      #angle_col=90, 
                      annotation_row = anno,
                      annotation_col = annocol,
                      annotation_colors = col_colors,
                      color =colorRampPalette(c("#0225A2", "white","#FF2400"))(30)
                      #color =colorRampPalette(c("#5BB93B","#030303", "#C22D1D"))(100)
)
font_family <- "serif"

p.heatmap@row_names_param$gp$fontfamily = font_family

p.heatmap@matrix_legend_param = gpar(title_gp = gpar(fontfamily = font_family), labels_gp =  gpar(fontfamily = font_family))

p.heatmap@top_annotation@anno_list[[1]]@name_param$gp$fontfamily = font_family  

p.heatmap@left_annotation@anno_list[[1]]@name_param$gp$fontfamily = font_family  

p.heatmap@left_annotation@anno_list[[2]]@legend_param = gpar(title_gp = gpar(fontfamily = font_family), labels_gp =  gpar(fontfamily = font_family))

p.heatmap
 
##########################################################

##########06  SparseCCA

#########################################################
 
################# Input genes and taxa matrix ###########

rownames(DM_cut_abandence)[c(10,17,21,27,28,31,40,47,48,52)] <- c("Carboxylic acid#",
                                                                  "Quinoline class#",
                                                                  "Carboxamide#",
                                                                  "Methlester#",
                                                                  "Glycylglycine#",
                                                                  "Brugierol",
                                                                  "Benzo#",
                                                                  "beta zearalenol" ,
                                                                  "Dihomoprostaglandin#",
                                                                  "Phosphoserine#")
X <- as.matrix(t(DM_cut_abandence))
Y <- as.matrix(t(genus4))


## Ensure same sampleIDs in both genes and microbes matrices
rownames(X) <- str_replace_all(rownames(X),
                               "-",".")
 

ASV.diff <- paste0("g__",cut.les$Taxa)

Y <- Y[rownames(X),ASV.diff]#ASV.diff

dim(Y)

stopifnot(all(rownames(X) == rownames(Y)))

#colnames(Y) <- sapply(str_split(lefse$res_diff$Taxa[index],"\\|"), "[",3)


## select tuning parameters using grid-search

X1 <- scale(X, scale = TRUE)

rownames(X1) <- str_replace_all(rownames(X1),
                                "-",".")

Y1 <- scale(Y, scale = TRUE)

X1 <- X1[rownames(Y1),]

stopifnot(all(rownames(X1) == rownames(Y1)))


scoreXcv <- c()
scoreYcv <- c()
penaltyX <- seq(0.1,0.4,length=10)
penaltyY <- seq(0.15,0.4,length=10)

corr_all <- matrix(nrow = length(penaltyX), ncol =  length(penaltyY))

num_samples <- nrow(X1)

start_time <- Sys.time()
for( i in 1:length(penaltyX)){
  for(j in 1:length(penaltyY)){
    
    for(k in 1:num_samples){
      
      print(paste0("Index: i = ",i,", j =", j," k = ",k)); flush.console()
      
      res <- CCA(X1[-k,],Y1[-k,], penaltyx = penaltyX[i], penaltyz = penaltyY[j], K=1, niter = 5, trace = F)
      
      scoreXcv[k] <- X1[k,]%*%res$u 
      scoreYcv[k] <- Y1[k,]%*%res$v
    }
    ## correlation between scores for X and Y for all held out samples.
    corr_all[i,j] = cor(scoreXcv,scoreYcv) 
  }
}
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste0("Time elapsed for all = ", time_elapsed))


row.names(corr_all) <- as.character(penaltyX)
colnames(corr_all) <- as.character(penaltyY)

corr_all_df <- as.data.frame(corr_all)
rownames(corr_all_df)
colnames(corr_all_df)


# find index with max absolute corr
bestpenalty <- which(abs(corr_all) == max(abs(corr_all)), arr.ind = TRUE)
bestpenalty

bestpenaltyX <- penaltyX[bestpenalty[1]]
bestpenaltyX 
bestpenaltyY <- penaltyY[bestpenalty[2]]
bestpenaltyY

## order abs corr to get top 10 corr
index <- order(abs(corr_all), decreasing = T)
abs(corr_all)[index][1:10] ## top 5 absolute corr

cca.k = 10 ## number of desired components

## Run sparse CCA using selected tuning param using permutation search
cca <- run_sparseCCA(as.matrix(Y1), as.matrix(X1), cca.k, bestpenaltyX, bestpenaltyY,
                     outputFile=paste0("/Users/huachunyin/Desktop/Xinqiao/微生物与代谢/CCA.output.",bestpenaltyX,"_",bestpenaltyY,".txt"))

## canonical correlation for each component:
cca[[1]]$cors

## average number of genes and microbes in resulting components
avg_meta <- get_avg_features(cca[[1]]$u, cca.k)
avg_meta

avg.microbes <- get_avg_features(cca[[1]]$v, cca.k)
avg.microbes


## Test significance of correlation using LOOCV

cca.k = 10
scoresXcv <- matrix(nrow = nrow(X), ncol = cca.k)
scoresYcv <-  matrix(nrow = nrow(Y), ncol = cca.k)
corr_pval <- c()
corr_r <- c()
for(i in 1:nrow(Y1)){ 
  res <- CCA(X[-i,],Y[-i,], penaltyx=bestpenaltyX, penaltyz=bestpenaltyY, K=cca.k, trace = F) 
  for(j in 1:cca.k){
    print(paste0("i = ", i," K = ", j)); flush.console()
    scoresXcv[i,j] <- X[i,]%*%res$u[,j]
    scoresYcv[i,j] <- Y[i,]%*%res$v[,j]
  }
}
for(j in 1:cca.k){
  # plot(scoresXcv,scoresYcv)
  corr <- cor.test(scoresXcv[,j],scoresYcv[,j])
  corr_pval[j] <- corr$p.value
  corr_r[j] <- corr$estimate
}
corr_pval

length(which(corr_pval < 0.1)) 
which(corr_pval < 0.1)


corr_padj <- p.adjust(corr_pval, method = "BH")
corr_padj
which(corr_padj < 0.1)

length(which(corr_padj < 0.05))
#8
## LOOCV corr
corr_r

outputFile <- paste0("/Users/huachunyin/Desktop/Xinqiao/微生物与代谢/crc_sparseCCA_summary_",bestpenaltyX,"_",bestpenaltyY,".txt")
sink(outputFile)
cat(paste0(" bestpenaltyX = ", bestpenaltyX, ", bestpenaltyY = ", bestpenaltyY))
cat(paste0("\n cor(Xu,Yv): \n"))
cat(paste0(signif(cca[[1]]$cors, digits = 4)))
cat(paste0("\n Avg. no. of metabolites across components = ",avg_meta))
cat(paste0("\n Avg. no. of microbes across components= ", avg.microbes))
cat(paste0("\n P-value for components (LOOCV): \n"))
cat(paste0(signif(corr_pval, digits = 4)))
cat(paste0("\n LOOCV corr: \n"))
cat(paste0(signif(corr_r, digits = 4)))
cat(paste0("\n No. of components with p-value < 0.1 = ", length(which(corr_pval < 0.1))))
cat(paste0("\n No. of components with p-value < 0.05 = ", length(which(corr_pval < 0.05))))
cat(paste0("\n No. of components with FDR < 0.1 = ", length(which(corr_padj < 0.1))))
cat(paste0("\n Significant components: \n" ))
cat(paste0(which(corr_padj < 0.1)))
sink()

write.table(cca[[2]], file = paste0("CCA_var_metabolits.txt"), sep="\t", row.names = T, col.names = NA )
write.table(cca[[3]], file = paste0("CCA_var_microbes.txt"), sep="\t", row.names = T, col.names = NA )

## only spit out significant components
sig <- which(corr_padj < 0.05)
dirname <- paste0("sig_meta_taxa_components_",bestpenaltyX,"_", bestpenaltyY,"_padj/")
ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
save_CCA_components(cca[[1]],sig,dirname)

####################

#network analysis

##################
#dirname <- "sig_meta_taxa_components_0.4_0.344444444444444_padj"
setwd(dirname)
metabolites.list <- c()
taxa.list <- c()
dir.list <- dir()

for(li in dir.list){
  
  a <- read.table(li, sep="\t",header = T)
  
  metabolites <- a$metabolites
  
  taxa <- a$taxa
  
  metabolites.list <- c(metabolites, metabolites.list)
  
  taxa.list <- c(taxa,taxa.list)
}

metabolites.list <- na.omit(metabolites.list)

taxa.list <- na.omit(taxa.list)

metabolites.list <- metabolites.list[!duplicated(metabolites.list)]

taxa.list <- taxa.list[!duplicated(taxa.list)]

#######################

## Lasso tutorial
## Sambhawa Priya

## Goal
## Our aim in this tutorial is to show how to run our lasso analysis on a small set of meabolites (~2-3 meabolites)
## to identify host meabolite-taxa associations for the demo dataset.


## Install and import libraries
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
  sapply(pkg, require, character.only = TRUE)
}

nfold <- 10
## Extract the features in the matrix

for(i in 1:length(metabolites.list)){
  m.name <- metabolites.list[i] 
  
  X2 <- X1[,m.name]
  
  Y2 <- Y1[,taxa.list]
  
  ## Make sure y_i is numeric before model fitting
  stopifnot(class(X2) == "numeric")
  
  ## Fit lasso CV model
  dim(Y2)
  
  fit.model <- fit.cv.lasso(Y2, X2,  kfold = nfold)
  
  bestlambda <- fit.model$bestlambda
  r.sqr <- fit.model$r.sqr ## note this will give us R^2 for the meabolite's final model fit using bestLambda
  ## This R^2 reflects final model R^2 for this meabolite using all the microbes in the model,
  ## and does not correspond to each meabolite-microbe pair.
  
  ## Estimate sigma and betainit using the estimated LOOCV lambda.
  ## Sigma is the standard deviation of the error term or noise.
  sigma.myfun <- estimate.sigma.loocv(as.matrix(Y2), X2, bestlambda, tol=1e-5)
  sigma <- sigma.myfun$sigmahat
  beta <- as.vector(sigma.myfun$betahat)[-1] ## remove intercept term
  sigma.flag <- sigma.myfun$sigmaflag
  
  ## Inference using lasso projection method, also known as the de-sparsified Lasso,
  ## using an asymptotic gaussian approximation to the distribution of the estimator.
  
  lasso.proj.fit <- lasso.proj(as.matrix(Y2), X2, multiplecorr.method = "BH", betainit = beta, sigma = sigma, suppress.grouptesting = T)
  ## A few lines of log messages appear here along with a warning about substituting sigma value (standard deviation of error term or noise)
  ## because we substituted value of sigma using our computation above.
  # Warning message:
  #   Overriding the error variance estimate with your own value.
  
  ## get 95% confidence interval (CI)
  lasso.ci <- as.data.frame(confint(lasso.proj.fit, level = 0.95))
  
  ## prep lasso output dataframe
  lasso.df <- data.frame(Metabolites = rep(m.name, length(lasso.proj.fit$pval)),
                         taxa = names(lasso.proj.fit$pval.corr),
                         r.sqr = r.sqr,
                         pval = lasso.proj.fit$pval,
                         ci.lower = lasso.ci$lower, ci.upper = lasso.ci$upper,
                         se=lasso.proj.fit$se,
                         adjP=lasso.proj.fit$pval.corr,
                         beta=lasso.proj.fit$betahat,
                         row.names=NULL)
  
  
  ## sort by p-value
  lasso.df <- lasso.df[order(lasso.df$pval),]
  lasso.df <- lasso.df[lasso.df$pval<0.05,]
  if(nrow(lasso.df)>0){
    write.table(lasso.df, file = "lasso_results1.txt", quote = FALSE, sep = "\t", 
                row.names = FALSE, append = TRUE)
  }else{
    next
  }
}

##################
cca_lass.re <- read.table("lasso_results1.txt",sep="\t",head=T,  stringsAsFactors = F)

cca_lass.re$adjP=p.adjust(cca_lass.re$pval, "BH")

cca_lass.re.adj <- cca_lass.re[cca_lass.re$adjP<0.05&abs(cca_lass.re$r.sqr)>0.3,]

unique(cca_lass.re.adj$Metabolites)
unique(cca_lass.re.adj$taxa)
openxlsx::write.xlsx(cca_lass.re.adj, file = "cca_lass.re.adj.xlsx")

##################
cca_lass.re.adj <- openxlsx::read.xlsx( "cca_lass.re.adj1.xlsx")

matrix <- xtabs(r.sqr ~ Metabolites + taxa, data = cca_lass.re.adj)
 