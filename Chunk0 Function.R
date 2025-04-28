####################################################################################################################################################################
#      ☆        % Project: Gut Microbiota Alterations in Acromegaly Patients Are Associated with Neutrophil Depletion-Induced Inflammation #
#   ☆ \|/ ☆    % Author: HuaChun Yin                                         
#  ☆  \|/  ☆   % Date: Apr. 4th, 2025                                  
# ☆   \|/   ☆  %                                                          
#  ☆  \|/  ☆   % Environment:   R version 4.4.2           
#  ☆ __|__ ☆   % EPlatform: Mac-IOS(64-bit)                                  
#                %                
################################################################################################################################################################### 


check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


############################################################################################################################

#Chunk 1 function

############################################################################################################################
filter_ASV <- function(ASV, qt){
  ASV.t <- ASV
  ASV.sd <- transform(as.data.frame(ASV.t), SD=apply(as.data.frame(ASV.t),1, sd, na.rm = TRUE))
  
  SD_quantile <- quantile(ASV.sd$SD) 
  SD_cutoff <- SD_quantile[qt] 
  ASV.sd <- ASV.sd[order(ASV.sd$SD, decreasing = T),]
  variable.ASV <- rownames(ASV.sd[ASV.sd$SD > SD_cutoff,])
  
  select <- which(rownames(ASV.t) %in% variable.ASV)
  ASV <- ASV.t[select,]
}

CompositionFigure <- function(x,y,n){ # x is the ASV matrix, y is the taxonomy matrix, change the metaphlan result to a composition table, select top n most abundant features
  #x
  #      A.1 A.2
  #ASV1  956  285
  #ASV2  199  814
  #ASV3    0    0
  #ASV4  232    0
  #ASV5    0  348
  #ASV6 1978    0
  
  #y
  #           phylum
  #ASV1     Proteobacteria
  #ASV2         Firmicutes
  #ASV3   Actinobacteriota
  #ASV4         Firmicutes
  #ASV5         Firmicutes
  #ASV6   Actinobacteriota
  
  
  
  percentages <- prop.table(as.matrix(t(x)), margin = 1) * 100
  
  
  result_df <- as.data.frame(t(percentages))
  
  new_df <- data.frame(matrix(ncol = ncol(y)+2, nrow = n*ncol(x)))
  
  colnames(new_df) <- c("Sample_ID","Proportion",colnames(y))
  
  
  # Sort each column individually and extract the top n rows"
  for (ii in 1:ncol(result_df)) {
    
    print(ii)
    col <-  colnames(result_df)[ii]
    
    sorted_col <- order(result_df[[col]], decreasing = TRUE)
    top_ten <- sorted_col[1:n]
    
    chunk <- seq(1, nrow(new_df)+n, n)
    
    new_df[chunk[ii]:c(chunk[ii+1]-1),"Proportion"] <- result_df[top_ten,col]
    new_df[chunk[ii]:c(chunk[ii+1]-1),"Sample_ID"] <- rep(col, n)
    
    for (i in 1:ncol(y)) {
      
      new_df[chunk[ii]:c(chunk[ii+1]-1),i+2] <- y[top_ten,i]
      
    }
    
    
  }
  
  
  
  #pal_igv("default", alpha = 0.66)(51)
  
  for (p1 in 1:ncol(y)) {
    
    p <- colnames(y)[p1]
    
    new_df1 <- new_df[,c("Sample_ID","Proportion",p)]
    
    colnames(new_df1) <- c("Sample_ID","Proportion", "Type")
    new_df1 <-  na.omit(new_df1)
    
    rhg_cols1 <-colorRampPalette(brewer.pal(8, "Set2"))(length(unique(new_df1$Type)))
    P.class <- ggplot(new_df1,aes(x=Sample_ID,y=Proportion,fill= Type))+
      geom_bar(stat="identity")+
      theme_classic()+
      scale_fill_manual(values = rhg_cols1)+
      labs(fill = p)+
      theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1,family="serif"))
    
    
    assign(paste0(p,"P"),P.class)
    
  }
  
  return(get(paste0(p,"P")))
  return(new_df)
  
  
}

get_anno_for_heatmap2<-function(annocol,annorow=NULL,color=NULL,only.color=F){
  require(plyr)
  require(stringr)
  if(is.null(color)){
    require(RColorBrewer)
    color=c(brewer.pal(12,"Set3"),brewer.pal(12,"Paired"),brewer.pal(8,"Set2"),brewer.pal(9,"Set1"),brewer.pal(8,"Dark2"))
  }
  
  annocolor=do.call(as.list,list(x=annocol))
  annocolor=lapply(annocolor,function(x){if(is.factor(x)){x=levels(x);a=color[1:length(x)];names(a)=x;return(a)}else{x=unique(x);a=color[1:length(x)];names(a)=x;return(a)}})
  if(!is.null(annorow)){
    annocolor.row<-do.call(as.list,list(x=annorow))
    annocolor.row=lapply(annocolor.row,function(x){if(is.factor(x)){x=levels(x);a=color[1:length(x)];names(a)=x;return(a)}else{x=unique(x);a=color[1:length(x)];names(a)=x;return(a)}})
  }else{annocolor.row=NULL}
  annocolor=c(annocolor,annocolor.row)
  annocolor_col<-as.list(annocol)
  annocolor_row<-as.list(annorow)
  annocolor<-c(annocol,annorow)
  annocolor<-lapply(annocolor,function(x){if(is.factor(x)){x=levels(x);return(x)}else{x=unique(x);return(x)}})
  annocolor<-do.call(c,annocolor)
  annocolor<-data.frame(var_name=as.factor(stringr::str_replace(names(annocolor),"[0-9]{1,}$","")),
                        var=annocolor,
                        color=color[1:length(annocolor)])
  annocolor<-split(annocolor,annocolor$var_name)
  annocolor<-lapply(annocolor,function(x){a=x$var;b=as.character(x$color);names(b)=a;return(b)})
  
  
  
  # if(only.color){
  #   anno_res<-annocolor
  # }else{anno_res<-list(annocol=annocol,
  #                      annorow=annorow,
  #                      annocolor=annocolor)}
  # return(anno_res)
  
  
  
}

# Calculating with many distance metrics.
dist_n <- function(x, mtd = "euclidean", p = NULL){
  
  if (!require("philentropy")) install.packages("philentropy")
  
  if (mtd == "maximum") {
    dist(x = x,
         method = mtd,
         diag = FALSE,
         upper = FALSE,
         p = p)
  }else{
    # library(philentropy)
    distance(x = x,
             method = mtd,
             p = p,
             test.na = TRUE,
             unit = "log",
             est.prob = NULL,
             use.row.names = TRUE,
             as.dist.obj = TRUE,
             diag = FALSE,
             upper = FALSE)
  }
}

get_tuning_params <- function(X, Z, outputFile = NULL, pdfFile = NULL, num_perm=25){
  perm.out <- CCA.permute(X,Z,typex="standard",typez="standard", nperms = num_perm)
  ## can tweak num. of permutations if needed (default = 25)
  
  ## get p-value for first left and right canonical covariates resulting from the selected
  ## tuning parameter value
  perm.out.pval <- perm.out$pvals[which(perm.out$penaltyxs == perm.out$bestpenaltyx)]
  
  if(!is.null(outputFile)){
    sink(outputFile)
    print(perm.out)
    print(paste0("P-value for selected tuning params:",perm.out.pval)) #0
    sink()
  }
  
  if(!is.null(pdfFile)){
    pdf(pdfFile)
    plot(perm.out)
    dev.off()
  }
  
  return(perm.out)
}

run_sparseCCA <- function(X, Z, CCA.K, penaltyX, penaltyZ, vInit=NULL, outputFile=NULL){
  CCA.out <-  CCA(X,Z,typex="standard",typez="standard",K=CCA.K,
                  penaltyx=penaltyX,penaltyz=penaltyZ,
                  v=vInit) ## standardize=T by default
  if(!is.null(outputFile)){
    sink(outputFile)
    print(CCA.out)
    sink()
  }
  
  ## add rownames to output factors
  rownames(CCA.out$u) <- colnames(X)
  rownames(CCA.out$v) <- colnames(Z)
  ## compute contribution of selected features to each of the samples.
  CCA_var_Meta <- X %*% CCA.out$u ## canonical variance for metabolites 
  CCA_var_microbes <- Z %*% CCA.out$v ## canonical variance for microbes
  
  return(list(CCA.out, CCA_var_Meta, CCA_var_microbes))
  
}

get_avg_features <- function(cca_cov, CCA.K){
  num_features <- 0
  for(k in 1:CCA.K){
    num_features <- num_features + length(which(cca_cov[,k]!=0))
  }
  avg_features <- num_features/CCA.K
}

save_CCA_components <- function(CCA.out, CCA.K, dirname){
  ## Print canonical covariates in files 
  for(i in CCA.K){
    #i <- 2 ##debug
    print(paste0("Writing significant component = ", i))
    selected_X <- which(CCA.out$u[,i]!=0) 
    selected_X <- rownames(CCA.out$u)[selected_X]
    coeff_X <- unname(CCA.out$u[selected_X,i])
    selected_Z <- which(CCA.out$v[,i]!=0)
    selected_Z <- rownames(CCA.out$v)[selected_Z]
    coeff_Z <- unname(CCA.out$v[selected_Z,i])
    ## Make all vectors of same length to avoid repetition of elements from shorter vectors.
    n <- max(length(selected_X), length(selected_Z))
    length(selected_X) <- n                      
    length(selected_Z) <- n
    length(coeff_X) <- n
    length(coeff_Z) <- n
    
    selected_XZ <- as.data.frame(cbind(metabolites = selected_Z, metabolites_coeff = coeff_Z,
                                       taxa = selected_X, taxa_coeff = coeff_X))
    write.table(selected_XZ, file=paste0(dirname,"metabolites_taxa_component_",i,".txt"), sep = "\t", col.names = NA)
  }
  
}
 
estimate.sigma.loocv <- function(x, y_i, bestlambda, tol) {
  
  
  ## Fit a lasso object
  lasso.fit = glmnet(x,y_i,alpha = 1) ## this is same as cv.fit$glmnet.fit from loocv code below.
  beta <- as.vector(coef(lasso.fit, s = bestlambda)) ## This gives coefficients of fitted model, not predicted coeff.
  # try(if(length(which(abs(beta) > tol)) > n) stop(" selected predictors more than number of samples! Abort function"))
  
  y = as.vector(y_i)
  
  yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
  ## predicted coefficients, same as coefficient of fitted model lasso. Either one is fine.
  # beta = predict(lasso.fit,s=bestlambda, type="coef")
  df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
  n = length(y_i)
  ss_res = sum((y - yhat)^2)
  
  if((n-df-1) >= 1) {
    sigma = sqrt(ss_res / (n-df-1))
    sigma.flag = 0
  } else{
    sigma = 1 ## conservative option
    # sigma = ss_res ## lenient option
    sigma.flag = 2
  }
  
  
  return(list(sigmahat = sigma, sigmaflag = sigma.flag, betahat = beta)) ## we return beta to be used later in hdi function.
  
}

fit.cv.lasso <- function(x, y_i, kfold){
  
  lambdas = NULL
  r.sqr.final <- numeric()
  r.sqr.final.adj <- numeric()
  r.sqr.CV.test <- numeric()
  lasso.cv.list <- list()
  
  ## glmnet CV
  cv.fit <- cv.glmnet(x, y_i, alpha=1, nfolds=kfold, type.measure = "mse", keep =TRUE, grouped=FALSE, standardize = T)  
  lambdas = data.frame(cv.fit$lambda,cv.fit$cvm)
  
  ## get best lambda -- lambda that gives min cvm
  bestlambda <- cv.fit$lambda.min
  bestlambda_index <- which(cv.fit$lambda == bestlambda)
  
  ## Get R^2 of final model
  final_model <- cv.fit$glmnet.fit
  r_sqr_final_model <- cv.fit$glmnet.fit$dev.ratio[bestlambda_index]
  
  ## Get adjusted R^2
  r_sqr_final_adj <- adj_r_squared(r_sqr_final_model, n = nrow(x), 
                                   p = sum(as.vector(coef(cv.fit$glmnet.fit, 
                                                          s = cv.fit$lambda.min)) > 0))
  
  return(list(bestlambda = bestlambda, r.sqr = r_sqr_final_model, 
              r.sqr.adj = r_sqr_final_adj
  ))
}

## functions to compute R2
r_squared <- function(y, yhat) {
  ybar <- mean(y)
  ## Total SS
  ss_tot <- sum((y - ybar)^2)
  ## Residual SS
  ss_res <- sum((y - yhat)^2)
  ## R^2 = 1 - ss_res/ ss_tot
  1 - (ss_res / ss_tot)
}
## Function for Adjusted R^2
## n sample size, p number of prameters
adj_r_squared <- function(r_squared, n, p) {
  1 - (1 - r_squared) * (n - 1) / (n - p - 1)
}

############################################################################################################################

#Chunk 3 function

############################################################################################################################




enrichment.kegg <- function(
    df = df, #gene list
    background.num=9501, #  DOI:10.1038/s41467-020-17916-9
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

if(is.null(class1)&is.null(class2)){
  pathwaymapping <- pathway_gene
}

if(is.null(class1)==FALSE&is.null(class2)){
  pathwaymapping <- pathway_gene[which(pathway_gene$level1_name%in%class1),]
}

if(is.null(class1)&is.null(class2)==FALSE){
  pathwaymapping <- pathway_gene[which(pathway_gene$level2_name%in%class2),]
}
if(is.null(class1)==FALSE&is.null(class2)==FALSE){
  pathwaymapping <- pathway_gene[which(pathway_gene$level1_name%in%class1),]
  
  pathwaymapping <- pathway_gene[which(pathway_gene$level2_name%in%class2),]
}

enrichres <- pathwaymapping
enrichres$exist_k <- sapply(enrichres$KOs, function(kos) {
  found_genes <- df[df %in% unlist(strsplit(kos, ","))]
  length(found_genes)
})

enrichres$gene_hit <- sapply(enrichres$KOs, function(kos) {
  found_genes <- df[df %in% unlist(strsplit(kos, ","))]
  paste(found_genes, collapse = ",")  #  
})

if(is.null(min_exist_KO)){
  min_exist_KO <- 3
  enrichres <- enrichres[enrichres$exist_k>min_exist_KO,]
}else{
  enrichres <- enrichres[enrichres$exist_k>min_exist_KO,]
}


M <-background.num
n <- length(df) #  


enrichres$p.value <- c(rep("NA",nrow(enrichres)))
for(i in 1:nrow(enrichres)){
  
  N <- enrichres[i,"K_num"]#  
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




netVisual_diffInteraction.YHC <- function (object, comparison = c(1, 2), measure = c("count", 
                                                    "weight", "count.merged", "weight.merged"), color.use = NULL, 
          color.edge = c("#b2182b", "#2166ac"), title.name = NULL, 
          sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, 
          top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, 
          vertex.size.max = 15, vertex.label.cex = 1, vertex.label.color = "black", 
          edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6, 
          label.edge = FALSE, edge.label.color = "black", edge.label.cex = 0.8, 
          edge.curved = 0.2, shape = "circle", layout = in_circle(), 
          margin = 0.2, arrow.width = 1, arrow.size = 0.2) 
{
  options(warn = -1)
  measure <- match.arg(measure)
  obj1 <- object@net[[comparison[1]]][[measure]]
  obj2 <- object@net[[comparison[2]]][[measure]]
  net.diff <- obj2 - obj1
  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  }
  else if (measure %in% c("weight", "weight.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                          df.net[["target"]]), sum)
    net[is.na(net)] <- 0
  }
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  net[abs(net) < stats::quantile(abs(net), probs = 1 - top, 
                                 na.rm = T)] <- 0
  g <- graph_from_adjacency_matrix(net, mode = "directed", 
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max + 
    5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, -atan(coords_scale[igraph::V(g), 
                                                                             2]/coords_scale[igraph::V(g), 1]), pi - atan(coords_scale[igraph::V(g), 
                                                                                                                                       2]/coords_scale[igraph::V(g), 1]))
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1], 
                               color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, 
                                               alpha.edge)
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * 
      edge.width.max
  }
  else {
    igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
  }

  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle <- NA
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,
                                                             1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  
  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x = 1:length(igraph::V(g)), 
                               direction = -1, start = 0)
  label.dist <- vertex.weight/max(vertex.weight) + 2
  plot(g, edge.curved = edge.curved, vertex.shape = shape, 
       layout = coords_scale, margin = margin, vertex.label.dist = label.dist, 
       vertex.label.degree = label.locs, vertex.label.family = "Helvetica", 
       edge.label.family = "Helvetica")
  if (!is.null(title.name)) {
    text(0, 1.5, title.name, cex = 1.1)
  }
  gg <- recordPlot()
  return(gg)
}

 