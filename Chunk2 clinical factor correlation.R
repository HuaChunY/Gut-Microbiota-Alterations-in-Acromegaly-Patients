####################################################################################################################################################################
#      ☆        % Project: Gut Microbiota Alterations in Acromegaly Patients Are Associated with Neutrophil Depletion-Induced Inflammation #
#   ☆ \|/ ☆    % Author: HuaChun Yin                                         
#  ☆  \|/  ☆   % Date: Apr. 4th, 2025                                  
# ☆   \|/   ☆  %                                                          
#  ☆  \|/  ☆   % Environment:   R version 4.4.2           
#  ☆ __|__ ☆   % EPlatform: Mac-IOS(64-bit)                                  
#                % CHUNK2:this script includes clinical information                        
################################################################################################################################################################### 
 
packages <- c("dplyr","doParallel", 
              "randomForest","stringr","ggrepel",
              "pROC","caret","tidyr",
              "ggplot2","ggpubr","scales",
              "psych","tidyr","forcats",
              "stringr","corrplot",
              "ComplexHeatmap","RColorBrewer",
              "ggsci","microeco",
              "tidyverse","magrittr",
              "data.table","foreach",
              "lme4","nlme",
              "factoextra","vegan",
              "ggalluvial","PMA",
              "ggrepel","scales","hdi",
              "stabs")  
#check.packages(packages)

lapply(packages, library, character.only = TRUE)
###########################################################################

# 01 clinical factor correlation

###########################################################################
 
all_patients <- openxlsx::read.xlsx("all_patients.xlsx",
                          startRow = 2,  
                          colNames = TRUE, 
                          detectDates = TRUE)

  
##############polyps.pie

polyps.pie <- all_patients[,c("Gender", "Group","intestinal.polyps")]


counts <- table(polyps.pie$Gender, polyps.pie$Group,polyps.pie$intestinal.polyps)

 
counts_df <- as.data.frame(counts)
 
colnames(counts_df) <- c("Gender", "Group", "intestinal.polyps", "Count")

counts_df$Gender <- factor(counts_df$Gender,levels=unique(counts_df$Gender))
counts_df$Group <- factor(counts_df$Group,levels=unique(counts_df$Group))

counts_df$GG <- paste0(counts_df$Gender,"-",counts_df$Group)
counts_df$GGpolyps <-  paste0(counts_df$GG,counts_df$intestinal.polyps)
counts_df$GG<- factor(counts_df$GG)

################
 
df <- counts_df[,1:4] 

summary_table <- df%>%
  group_by(Gender, intestinal.polyps) %>%
  summarise(Count = sum(Count)) %>%
  spread(intestinal.polyps, Count, fill = 0) %>%
  as.data.frame()

 
table <- matrix(c(summary_table$Y[1], summary_table$N[1], 
                  summary_table$Y[2], summary_table$N[2]),
                nrow = 2, byrow = TRUE)
dimnames(table) <- list(Gender = c("Female", "Male"), polyps = c("Y", "N"))

#  
result <- fisher.test(table)
result
result <- chisq.test(table)
result

summary_table <- df %>%
  group_by(Gender, Group) %>%
  summarise(Count = sum(Count)) %>%
  spread(Group, Count, fill = 0) %>%
  as.data.frame()

 
table <- matrix(c(summary_table$GHPA[1], summary_table$Ctrl[1], 
                  summary_table$GHPA[2], summary_table$Ctrl[2]),
                nrow = 2, byrow = TRUE)
dimnames(table) <- list(Gender = c("Female", "Male"), Group = c("GHPA", "Ctrl"))

 
result <- fisher.test(table)
result
result <- chisq.test(table)
result

summary_table <- df %>%
  group_by(intestinal.polyps, Group) %>%
  summarise(Count = sum(Count)) %>%
  spread(Group, Count, fill = 0) %>%
  as.data.frame()

table <- matrix(c(summary_table$GHPA[1], summary_table$Ctrl[1], 
                  summary_table$GHPA[2], summary_table$Ctrl[2]),
                nrow = 2, byrow = TRUE)
dimnames(table) <- list(polyps = c("N", "Y"), Group = c("GHPA", "Ctrl"))
result <- fisher.test(table)
result
result <- chisq.test(table)
result

##################
dat1 <-  counts_df %>%
  select(GG,Count) %>% group_by(GG) %>% 
  summarise(total = sum(Count)) %>%
  ungroup() %>% 
  mutate(perc = total/sum(total),
         y = cumsum(total) - 0.5*total,
         label = paste0(GG,"(",percent(round(perc,3)),")"))


dat2 <-  counts_df %>%
  select(GGpolyps,Count) %>% group_by(GGpolyps) %>% 
  summarise(total = sum(Count)) %>%
  ungroup() %>% 
  mutate(perc = total/sum(total),
         y = cumsum(total) - 0.5*total)

dat2$label2 <-  sapply(str_split(dat2$GGpolyps,"N"), "[",1) 
dat2$label2 <-  sapply(str_split(dat2$label2,"Y"), "[",1) 

dat2 <-  dat2 %>%
  select(GGpolyps,total,label2,y,perc) %>% group_by(label2) %>% 
  mutate(perc2 = total/sum(total),
         label = paste0(GGpolyps,"(",percent(round(perc2,3)),")"))
#############################
table <- as.data.frame(dat2[,c(1:3)])
table$GGpolyps <- sapply(strsplit(table$GGpolyps, "-.{4}"), function(x) if(length(x) > 1) x[2] else NA)

summary_table <- table %>%
  group_by(GGpolyps, label2) %>%
  summarise(Count = sum(total)) %>%
  spread(label2, Count, fill = 0) %>%
  as.data.frame()

table <- matrix(c(summary_table$`Female-GHPA`[1], summary_table$`Female-Ctrl`[1], 
                  summary_table$`Male-GHPA`[1], summary_table$`Male-Ctrl`[1], 
                  summary_table$`Female-GHPA`[2], summary_table$`Female-Ctrl`[2], 
                  summary_table$`Male-GHPA`[2], summary_table$`Male-Ctrl`[2]),
                nrow = 2,byrow = TRUE)
dimnames(table) <- list(polyps = c("N","Y"), label = c("Female-GHPA", "Female-Ctrl",
                                                             "Male-GHPA","Male-Ctrl"))
result <- fisher.test(table)
result
result <- chisq.test(table)
result


dat1$GG <- factor(dat1$GG,levels = c("Female-GHPA",
                                     "Male-GHPA",
                                     "Female-Ctrl",
                                     "Male-Ctrl"))

dat2$GGpolyps <- factor(dat2$GGpolyps,levels = c("Female-GHPAN",
                                                 "Female-GHPAY",
                                                 "Male-GHPAN",
                                                 "Male-GHPAY",
                                                 "Female-CtrlN",
                                                 "Female-CtrlY",
                                                 "Male-CtrlN",
                                                 "Male-CtrlY"))
# 绘制双层饼图
pie <- ggplot()+
  geom_bar(data = dat2,aes(x=2,y=total,fill=fct_reorder(GGpolyps,y,.desc = TRUE)),stat="identity",width = 1,color="white")+  
  geom_text_repel(data =dat2,size=4,direction = "x",point.padding = 0,box.padding = 0,nudge_x = .4,aes(x=2,y=as.numeric(y),family = 'serif',label=label))+
  geom_bar(data = dat1,aes(x=1,y=total,fill=fct_reorder(GG,y,.desc = TRUE)),stat="identity",width = 1,color="white")+
  geom_text_repel(data =dat1,size=4,direction = "x",point.padding = 0,box.padding = 0,nudge_x = .4,aes(x=1,y=as.numeric(y),family = 'serif',label=label))+
  coord_polar(theta = "y")+
  scale_y_continuous(labels = scales::percent) +
  theme_void()+  
  scale_fill_simpsons() +  
  theme(legend.position = 'none')  

print(pie)


#########GH level

clinic.p <- c("GH","IGF1","BMI","Age","Glycated_hemoglobin","Triglycerides","Cholesterol","HDL","LDL")
for(i in 1:9){
  
  name <- clinic.p[i]
  
  
  GH_patients1 <- all_patients[,c("Group",name)]
  colnames(GH_patients1) <- c("Group","Value")
  
  GH_patients1 <- na.omit(GH_patients1)
  GHN <- length(which(GH_patients1$Group%in%"GHPA"))
  nonN <- length(which(GH_patients1$Group%in%"Ctrl"))
  y.min <- min(GH_patients1$Value)-1
  y.max <-max(GH_patients1$Value)+1.3
  lab.max <- y.max-0.2
  sig<-compare_means(Value~Group, 
                     data=GH_patients1,
                     method = "wilcox.test",
                     p.adjust.method = "BH") #t.test
  print(name)
  print(sig$p.signif)
  GH_patients1$Group <- factor(GH_patients1$Group, levels=c("Ctrl","GHPA"))
  
  P <- ggplot(data=GH_patients1, 
              aes(x=Group,y=Value))+
    annotate("text", x=1.5, y=lab.max, label= sig$p.signif)+
    geom_boxplot(aes(color=Group),size=0.8, outliers=F, position=position_dodge(1))+
    geom_jitter(alpha=0.5,aes(color=Group),
                position=position_jitterdodge(jitter.width = 0.3, 
                                              jitter.height = 0, 
                                              dodge.width = 0.4
                ))+
    scale_color_manual(values = c("#5E7BBB","#EF9703"))+
    # facet_wrap(~ Name, scales="free")+
    theme_bw()+
    ylab("secretion levels")+
    labs(title=paste0(name," level"))+
    scale_x_discrete(labels=c(paste0("Ctrl (N=",nonN,")"),paste0("GHPA (N=",GHN,")")))+
    ylim(y.min,y.max)+
    theme(text = element_text(family="serif",size = 8),#"Times New Roman","serif"
          panel.grid = element_blank(),
          axis.ticks.x.bottom = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust=0.5),
          axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=10),
          legend.position = "none"
    )
  
  assign(name,P)
  
}

ggarrange( GH, IGF1, nrow = 1,ncol = 2) 

############
all_patients$Gender <- gsub("Female", "0",all_patients$Gender)

all_patients$Gender <- gsub("Male", "1",all_patients$Gender)

all_patients$Gender <- as.numeric(all_patients$Gender)

all_patients$intestinal.polyps <- gsub("Y", "1",all_patients$intestinal.polyps)
all_patients$intestinal.polyps <- gsub("N", "0",all_patients$intestinal.polyps)
all_patients$intestinal.polyps <- as.numeric(all_patients$intestinal.polyps)


all_patients$Hypertension <- gsub("Y", "1",all_patients$Hypertension)
all_patients$Hypertension <- gsub("N", "0",all_patients$Hypertension)
all_patients$Hypertension <- as.numeric(all_patients$Hypertension)

all_patients$Diabetes <- gsub("Y", "1",all_patients$Diabetes)
all_patients$Diabetes <- gsub("N", "0",all_patients$Diabetes)
all_patients$Diabetes <- as.numeric(all_patients$Diabetes)

############  GHPA patients
 
GH_patients1 <- all_patients[which(all_patients$Group%in%"GHPA"),] 

GH_patients2_1 <- GH_patients1%>%select(-"ID")

for (ii in 2:9) {
  
  name <- clinic.p[ii]
  
  
  GH_patients2 <- GH_patients2_1[,c("GH",name)]
  colnames(GH_patients2) <- c("GH","Value")
  
  GH_patients2 <- na.omit(GH_patients2)
  nrow <- nrow(GH_patients2)
 
  P.cor <- ggplot(GH_patients2, aes(x = log2(Value), y = log2(GH))) +
    geom_point() + 
    geom_smooth(method = "lm",color="#EF9703")+
    stat_cor(method = "spearman")+
    xlab(paste0("Log2(",name,")"))+
    theme_classic()+
    theme(text = element_text(family="serif",size = 9),
          legend.position = "none")
  
  assign(paste0(name,".cor"),P.cor)
  print(paste0(name,".cor"))
}


##################Figure S1 A
P <- ggarrange( 
  Glycated_hemoglobin,
  Triglycerides,
  Cholesterol,
  LDL,
  IGF1.cor, 
  Glycated_hemoglobin.cor,
  Triglycerides.cor,
  Cholesterol.cor, 
  LDL.cor,
  nrow = 3,ncol = 4)
 
########################
# DM and GH clinical factor
######################## 

#DM_cut_abandence  <- openxlsx::read.xlsx("DM_cut_abandence.xlsx",rowNames = TRUE)

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




###################################neutrophil

neutrophil.dt <- all_patients[,c("ID","Neutrophil.(10^9/L)" )]

DM1 <- DM_cut_abandence[which(rownames(DM_cut_abandence)%in%"Taurine"),]
DM1 <- t(DM1)
neutrophil <- cbind(neutrophil.dt,DM1)
neutrophil <- neutrophil[,-1]
colnames(neutrophil) <- c("Neutrophil","Taurine")


corr.test(neutrophil[,1],
          neutrophil[,2],
          method = "spearman",adjust="BH")


ggplot(neutrophil, aes(x = Taurine, y = Neutrophil)) +
  geom_point() +  
  geom_smooth(method = "lm",color="#EF9703")+
  stat_cor(method = "spearman")+
  xlab("Taurine")+
  theme_classic()+
  theme(text = element_text(family="serif",size = 9),
        legend.position = "none")

####################################

all_patients1 <- all_patients[,c("GH", "Age","Gender", "Hypertension","Diabetes","Glycated_hemoglobin","Triglycerides",
                                 "Cholesterol","HDL","LDL","IGF1","BMI")]
rownames(all_patients1) <- all_patients$ID

DM_cut_abandence1 <- DM_cut_abandence 

for(i in 1:ncol(all_patients1)){
 
  a <- colnames(all_patients1)[i]
  
  NA.name <- rownames(all_patients1)[which(is.na(all_patients1[,a]))]
  
  patients <- as.matrix(na.omit(all_patients1[,a]))

  
  NA.pos<- which(colnames(DM_cut_abandence1)%in%NA.name)
 
   print(i)
  
  
  if(length(NA.pos)>0){
    
    P.cor<- corr.test(t(DM_cut_abandence1[,-NA.pos]),
                      patients,
                      method = "spearman",adjust="BH")
    
  }
  else{
    
    P.cor<- corr.test(t(DM_cut_abandence1),
                      patients,
                      method = "spearman",adjust="BH")
    
  }
   
  assign(a,P.cor)
}



cor.matrix <- cbind(GH$r,
                    Age$r,
                    Gender$r, 
                    Hypertension$r,
                    Diabetes$r,
                    Glycated_hemoglobin$r,
                    Triglycerides$r,
                    Cholesterol$r,
                    HDL$r,
                    LDL$r,
                    IGF1$r,
                    BMI$r)

colnames(cor.matrix) <-c("GH", 
                         "Age",
                         "Gender",
                         "Hypertension",
                         "Diabetes",
                         "Glycated hemoglobin",
                         "Triglycerides",
                         "Cholesterol",
                         "HDL",
                         "LDL",
                         "IGF1",
                         "BMI")
 


adjp.matrix <- cbind(GH$p.adj,
                  Age$p.adj,
                  Gender$p.adj,
                  Hypertension$p.adj,
                  Diabetes$p.adj,
                  Glycated_hemoglobin$p.adj,
                  Triglycerides$p.adj,
                  Cholesterol$p.adj,
                  HDL$p.adj,
                  LDL$p.adj,
                  IGF1$p.adj,
                  BMI$p.adj)

 
colnames(adjp.matrix)<-colnames(cor.matrix) 
 
cor.matrix <- t(cor.matrix)
adjp.matrix <- t(adjp.matrix)
 
par(family= "serif")# 
 corrplot(cor.matrix, 
         p.mat = round(adjp.matrix,2), 
         sig.level = 0.05,
         pch.cex=1,
         tl.srt=45,
         addgrid.col="grey90",
         insig = 'label_sig', 
         tl.col = "black",
         method = 'square', 
         #addCoef.col = 'black', 
          col = rev(COL2(n=100)), 
         number.cex = 0.4, 
         addCoef.col = NULL,
         tl.cex = 0.5,
         cl.cex = 0.4)
 # 7.95*4.56
#########################16s and clinic factor
 
 ####################
 
feature_table <- read.csv("extdata/feature_table.csv",row.names = 1)
rownames(feature_table) <- sapply(str_split(rownames(feature_table),"g__"), "[",2)

cut.les <- openxlsx::read.xlsx("16s/lefse_result_cut.xlsx")

#cut.les$Taxa <- sapply(str_split(cut.les$Taxa,"\\|"), "[",2)
#cut.les$Taxa <- sapply(str_split(cut.les$Taxa,"g__"), "[",2)

feature_table <- feature_table[which(rownames(feature_table)%in%cut.les$Taxa,),]

#########
###################################neutrophil

neutrophil.dt <- all_patients[,c("ID","Neutrophil.(10^9/L)" )]

levels.Bilophila <- feature_table[which(rownames(feature_table)%in%"Bilophila"),]
levels.Bilophila <- as.data.frame(t(levels.Bilophila))

levels.Bilophila$ID <- rownames(levels.Bilophila)

levels.Bilophila$ID <- gsub("\\.", "-", levels.Bilophila$ID)

neutrophil <- merge(neutrophil.dt,levels.Bilophila,by="ID")
neutrophil <- neutrophil[,-1]
colnames(neutrophil) <- c("Neutrophil","Bilophila")


corr.test(neutrophil[,1],
          neutrophil[,2],
          method = "spearman",adjust="BH")


ggplot(neutrophil, aes(x = log2(Bilophila+1), y = Neutrophil)) +
  geom_point() +  
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman")+
  xlab("Bilophila")+
  theme_classic()+
  theme(text = element_text(family="serif",size = 9),
        legend.position = "none")

#############
########################

group <- substr(neutrophil.dt$ID,1,1)

group <- str_replace_all(group,
                         "A","GHPA")
group <- str_replace_all(group,
                         "B","Ctrl")

neutrophil.dt$Group <- factor(group, levels=c("Ctrl","GHPA"))

colnames(neutrophil.dt) <- c("ID","neutrophil","Group")

sig<-compare_means(neutrophil~Group, 
                   data=neutrophil.dt,
                   method = "wilcox.test",
                   p.adjust.method = "BH",) #t.test

y.max <-max(neutrophil.dt$neutrophil)+1.3

lab.max <- y.max-0.2

# 计算 SEM
sem <- function(x) { sd(x) / sqrt(length(x)) }
# stat_summary 中的 fun.min / fun.max
summary_min <- function(x) { mean(x) - sem(x) }
summary_max <- function(x) { mean(x) + sem(x) }


ggplot(data = neutrophil.dt, aes(x=Group,y=neutrophil)) + 
  stat_summary(
    aes(color = Group),
    geom = "errorbar",
    width = 0.4,       
    size = 1,          
    fun.min = summary_min,
    fun.max = summary_max ) +
  geom_bar(
    aes(color = Group), 
    stat = "summary", 
    fun = "mean", 
    fill = "white",    
    linewidth = 1, 
    width = 0.5 ) +
  geom_jitter(alpha=0.5,aes(color=Group),
              position=position_jitterdodge(jitter.width = 0.3, 
                                            jitter.height = 0, 
                                            dodge.width = 0.4
              )) + 
  annotate("segment", x = 1, xend = 2, y = 9.4, yend =9.4, size = 1) +
  annotate("text", x = 1.5, y = 9, label = "**", vjust = -1, size = 5)+
  # 坐标轴和标签设置
  ylab("Abundance levels")+
  labs(title="Neutrophil cells")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10)) +
  scale_color_manual(values = c("#5E7BBB","#EF9703")) +
  theme_classic() +
  theme(text = element_text(family="serif",size = 12),
    legend.position = "none",
    axis.title = element_text( color = "black"),
    axis.text = element_text( color = "black"),
    axis.line = element_line(color = "black", linewidth = 1.1),
    axis.ticks = element_line(color = "black", linewidth = 1.1),
    axis.ticks.length = unit(0.07, "in") )

#2*2.8

#######
#########alpha diversity
#alpha_diversity <- openxlsx::read.xlsx("16s/alpha_re.xlsx")
#feature_table <-alpha_diversity[,c(1,4:8)]

#rownames(feature_table) <- feature_table$Sample
#feature_table <- feature_table[,-1]
#feature_table <- feature_table%>%mutate(across(everything(), as.numeric))
 

feature_table <- t(feature_table)

rownames(all_patients1) <- str_replace_all(rownames(all_patients1),"-",".")

for(i in 1:ncol(all_patients1)){
  
  
  a <- colnames(all_patients1)[i]
  
  NA.name <- rownames(all_patients1)[which(is.na(all_patients1[,a]))]
  
  patients <- as.matrix(na.omit(all_patients1[,a]))
 
  NA.pos<- which(rownames(feature_table)%in%NA.name)
  
  print(a)
  
  
  if(length(NA.pos)>0){
    
    P.cor<- corr.test(feature_table[-NA.pos,],
                      as.matrix(patients),
                      method = "spearman",adjust="BH")
    
  }
  else{
    
    P.cor<- corr.test(feature_table,
                      patients,
                      method = "spearman",adjust="BH")
    
  }
  
  assign(a,P.cor)
}

cor.matrix <- cbind(GH$r,
                    Age$r,
                    Gender$r,
                    Hypertension$r,
                    Diabetes$r,
                    Glycated_hemoglobin$r,
                    Triglycerides$r,
                    Cholesterol$r,
                    HDL$r,
                    LDL$r,
                    IGF1$r,
                    BMI$r)

colnames(cor.matrix) <-c("GH", 
                         "Age",
                         "Gender",
                         "Hypertension",
                         "Diabetes",
                         "Glycated hemoglobin",
                         "Triglycerides",
                         "Cholesterol",
                         "HDL",
                         "LDL",
                         "IGF1",
                         "BMI")



adjp.matrix <- cbind(GH$p.adj,
                     Age$p.adj,
                     Gender$p.adj,
                     Hypertension$p.adj,
                     Diabetes$p.adj,
                     Glycated_hemoglobin$p.adj,
                     Triglycerides$p.adj,
                     Cholesterol$p.adj,
                     HDL$p.adj,
                     LDL$p.adj,
                     IGF1$p.adj,
                    BMI$p.adj)

 
 
colnames(adjp.matrix)<-colnames(cor.matrix) 
 
par(family= "serif")# 
corrplot(cor.matrix, 
         p.mat = round(adjp.matrix,2), 
         sig.level = 0.05,
         pch.cex=1,
         tl.srt=90,
         hclust.method = c("ward.D"),
         addgrid.col="grey90",
         insig = 'label_sig', 
         tl.col = "black",
         method = 'square', 
         #addCoef.col = 'black', 
         col = rev(COL2(n=100)), 
         number.cex = 0.2, 
         addCoef.col = NULL,
         tl.cex = 0.6,
         cl.cex = 0.6)

 

########################Taurine

Taurine.plot <- as.data.frame(t(DM_cut_abandence["Taurine",]))

Taurine.plot$Group <- factor(group, levels=c("GHPA","Ctrl"))

sig<-compare_means(Taurine~Group, 
                   data=Taurine.plot,
                   method = "wilcox.test",
                   p.adjust.method = "BH",) #t.test

y.max <-max(Taurine.plot$Taurine)+1.3

lab.max <- y.max-0.2

Taurine.p <- ggplot(data=Taurine.plot, 
                    aes(x=Group,y=Taurine))+
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
  labs(title="Taurine")+
  scale_x_discrete(labels=c( "Ctrl (N=25)" , "GHPA (N=35)"))+
  theme(text = element_text(family="serif",size = 8),#"Times New Roman","serif"
        panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=10),
        legend.position = "none"
  )

###############################################

#################################################

###RF classifier

##############################################
df1 <- all_patients[,c("ID" ,"Neutrophil.(10^9/L)")]

colnames(df1) <- c("ID","neutrophil")

Taurine.plot <- as.data.frame(t(DM_cut_abandence["Taurine",]))

Taurine.plot$ID <- rownames(Taurine.plot)


Bilophila.plot <- as.data.frame(t(feature_table["Bilophila",]))
 
rownames(Bilophila.plot) <- gsub("\\.", "-", rownames(Bilophila.plot))

Bilophila.plot$ID <- rownames(Bilophila.plot)

mergedf <- merge(df1,Taurine.plot,by="ID")

mergedf <- merge(mergedf,Bilophila.plot,by="ID")

rownames(mergedf) <- mergedf$ID

mergedf$ID <- substr(mergedf$ID,1,1)

mergedf$ID <- str_replace_all(mergedf$ID,
                         "A","GHPA")
mergedf$ID <- str_replace_all(mergedf$ID,
                         "B","Ctrl")

mergedf$ID <- factor(mergedf$ID, levels=c("GHPA","Ctrl"))

mergedf$Bilophila <- as.character(mergedf$Bilophila)

mergedf$Bilophila <- as.numeric(mergedf$Bilophila)


set.seed(114)
 

colnames(mergedf)[1] <- "class"

train_idx <- createDataPartition(mergedf$class, p = 0.7, list = FALSE)

train_data <- mergedf[train_idx, ]# c(1,2)

test_data <- mergedf[-train_idx,]

 
ctrl <- trainControl(method = "repeatedcv",  
                     number = 10,           
                     repeats = 5,         
                     classProbs = TRUE, 
                     summaryFunction = twoClassSummary,   
                     sampling = "up",    
                     savePredictions = "final")

model <- train(class ~ ., data = train_data, 
               method = "rf", 
               trControl = ctrl,
               tuneGrid = expand.grid(mtry = c(2, 3, 5, 7)),   
               ntree = 1000,   
               metric = "ROC")   

test_prob <- predict(model, newdata = test_data, type = "prob")[, 2]

 
 
roc <- roc(test_data$class, test_prob)

pROC::ci(roc)

best <- ci.coords(roc, "best",ret = c("specificity", "sensitivity"),best.policy = "random") 

best  

#windowsFonts(A = windowsFont("Times New Roman"))

plot.roc(roc, 
         col = "#377EB8", 
         legacy.axes = TRUE, 
         percent = TRUE, 
         xlab = "False Positive Percentage", 
         ylab = "TRUE Positive Percentage", 
         lwd = 4, 
         print.auc = TRUE, 
         print.auc.x = 0.45, 
         print.auc.y = 0.45,
         auc.polygon = TRUE, 
         #auc.polygon.col = "#377EB822"
         ) 

plot.roc(roc, 
         col = "#84CC8C", 
         legacy.axes = TRUE, 
         percent = TRUE, 
         xlab = "False Positive Percentage", 
         ylab = "TRUE Positive Percentage", 
         lwd = 4, 
         print.auc = TRUE, 
         print.auc.x = 0.85, 
         print.auc.y = 0.85,
         auc.polygon = TRUE, 
         #auc.polygon.col = "#377EB822", 
         add = TRUE) 
plot.roc(roc, 
         col = "#FFA631",
         legacy.axes = TRUE, 
         percent = TRUE, 
         xlab = "False Positive Percentage", 
         ylab = "TRUE Positive Percentage", 
         lwd = 4, 
         print.auc = TRUE, 
         print.auc.x = 0.85, 
         print.auc.y = 0.85,
         auc.polygon = TRUE, 
         #auc.polygon.col = "#377EB822", 
         add = TRUE) 

plot.roc(roc, 
         col = "#C5A5D2", 
         legacy.axes = TRUE, 
         percent = TRUE, 
         xlab = "False Positive Percentage", 
         ylab = "TRUE Positive Percentage", 
         lwd = 4, 
         print.auc = TRUE, 
         print.auc.x = 0.85, 
         print.auc.y = 0.85,
         auc.polygon = TRUE, 
         #auc.polygon.col = "#377EB822", 
         add = TRUE) 

legend("bottomright", 
       legend = c("Bilophila", "Taurine", "Neutrophil", "Conbined"), 
       col =c("#377EB8", 
              "#84CC8C",
              "#FFA631",
              "#C5A5D2"), 
       lwd = 4)

par(pty = "m")

