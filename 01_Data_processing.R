# R script for RNA-seq data processing
# Shuichi N Kudo
# Last update 2025. 2. 26
# R version 4.3.1
#
# functions ####################################################################
reconst_GE_matrix = function(dataset,scale="Expected_count"){
  if(scale=="Expected_count"){
    matrix_data = select(dataset,Sample_id,Gene_id,Expected_count) %>% pivot_wider(names_from=Sample_id,id_cols=Gene_id,values_from=Expected_count) %>% data.frame(row.names=1)
    
  } else if (scale=="TPM"){
    matrix_data = select(dataset,Sample_id,Gene_id,TPM) %>% pivot_wider(names_from=Sample_id,id_cols=Gene_id,values_from=TPM) %>% data.frame(row.names=1)
    
  } else if (scale=="RPK"){
    matrix_data = select(dataset,Sample_id,Gene_id,RPK) %>% pivot_wider(names_from=Sample_id,id_cols=Gene_id,values_from=RPK) %>% data.frame(row.names=1)
  }
  return(matrix_data)
}

sort_mat=function(GE_matrix_A){
  Sp_list = substr(colnames(GE_matrix_A),1,2)
  Sp_unq = unique(Sp_list)
  n_Sp = length(Sp_unq)
  for(i in 1:length(Sp_unq)){
    GE_matrix_Sp = GE_matrix_A[,Sp_list==Sp_unq[i]]
    Tissue_list = substr(colnames(GE_matrix_Sp),4,4)
    # separate matrix (tissues)
    mat_L = GE_matrix_Sp[,Tissue_list=="L"]
    mat_B = GE_matrix_Sp[,Tissue_list=="B"]
    mat_F = GE_matrix_Sp[,Tissue_list=="F"]
    Dates_L = substr(colnames(mat_L),1,6)
    Dates_B = substr(colnames(mat_B),1,6)
    Dates_F = substr(colnames(mat_F),1,6)
    # ordering matrix of each tissue
    order_L = order(Dates_L)
    order_B = order(Dates_B)
    order_F = order(Dates_F)
    mat_L_sorted = mat_L[,order_L]
    mat_B_sorted = mat_B[,order_B]
    mat_F_sorted = mat_F[,order_F]
    mat_Sp = cbind(mat_L_sorted, mat_B_sorted)
    mat_Sp= cbind(mat_Sp, mat_F_sorted)
    if(i==1){
      mat_ret = mat_Sp
    } else {
      mat_ret = cbind(mat_ret,mat_Sp)
    }
  }
  return(mat_ret)
}


# (1) Filtering of low expression gene #########################################
Sp_ls = c("Le","Lg","Qg","Qa")

# calculation of RPK
for(i in 1:length(Sp_ls)){
  Sp = Sp_ls[i]
  raw_data = read.csv(paste("data/dataset_raw_",Sp,".csv",sep=""),row.names=1)
  raw_data["length0_count"] = raw_data$Expected_count
  raw_data[raw_data$Effective_length==0,]$length0_count = 0
  raw_data["RPK"] = raw_data$length0_count/(raw_data$Effective_length/1000)
  raw_data[is.na(raw_data$RPK),]$RPK = 0
  rpk_matrix = sort_mat(reconst_GE_matrix(raw_data,scale="RPK"))
  write.csv(rpk_matrix,paste("data/datamatrix_rpk_",Sp,".csv",sep=""))
}

# calculate sum of RPK per gene within all samples for each species
# filtering condition
cutoff = c(0,0.1,0.5,1,2)

for(i in 1:length(Sp_ls)){
  Sp = Sp_ls[i]
  rpk_matrix = read.csv(paste("data/datamatrix_rpk_",Sp,".csv",sep=""),row.names=1)
  n_tot = dim(rpk_matrix)[1]
  filtering_table = data.frame(matrix(NA,nrow=n_tot,ncol=(length(cutoff)+5)))
  colnames(filtering_table) = c("Mean","N_sample","N_leaf","N_bud","N_flower",paste("c",cutoff,sep=""))
  row.names(filtering_table) = row.names(rpk_matrix)
  sample_list = colnames(rpk_matrix)
  tissue_list = substr(sample_list,4,4)
  N_L = length(tissue_list[tissue_list=="L"])
  N_B = length(tissue_list[tissue_list=="B"])
  N_F = length(tissue_list[tissue_list=="F"])
  mean_vec = apply(rpk_matrix,1,function(x){mean(x,na.rm=T)})
  filtering_table["Mean"] = mean_vec
  filtering_table["N_sample"] = dim(rpk_matrix)[2]
  filtering_table["N_leaf"] = N_L
  filtering_table["N_bud"] = N_B
  filtering_table["N_flower"] = N_F
  for(j in 1:length(cutoff)){
    c = cutoff[j]
    filter_c_vec = numeric(n_tot)
    if(c==0){
      filter_c_vec[mean_vec>c] = 1
    } else {
      filter_c_vec[mean_vec>=c] = 1
    }
    filtering_table[,j+5] = filter_c_vec
  }
  write.csv(filtering_table,paste("tables/filtering_result_",Sp,".csv",sep=""))
}

# summarize the filtering result
filtering_summary = data.frame(matrix(NA,ncol=4,nrow=length(cutoff)))
colnames(filtering_summary) = Sp_ls
row.names(filtering_summary) = paste("c",cutoff,sep="")
for(i in 1:length(Sp_ls)){
  Sp = Sp_ls[i]
  filtering_table = read.csv(paste("tables/filtering_result_",Sp,".csv",sep=""),row.names=1)
  N_tot = dim(filtering_table)[1]
  print(paste(Sp,N_tot))
  count_table = filtering_table[,6:dim(filtering_table)[2]]
  sum_vec = apply(count_table,2,sum)
  removed_count = N_tot - sum_vec
  removed_percent = 100*removed_count/N_tot
  filtering_summary[,i] = paste(removed_count," (",round(removed_percent,1)," %)",sep="")
}
write.csv(filtering_summary,"tables/filtering_result_summary.csv")

# final output
Sp_ls = c("Le","Lg","Qg","Qa")
cutoff = 1
for(i in 1:length(Sp_ls)){
  Sp = Sp_ls[i]
  filtering_table = read.csv(paste("tables/filtering_result_",Sp,".csv",sep=""),row.names=1)
  rpk_matrix = read.csv(paste("data/datamatrix_rpk_",Sp,".csv",sep=""),row.names=1)
  head(filtering_table)
  
  filtering_table_in = subset(filtering_table,filtering_table[paste("c",cutoff,sep="")]==1)
  sum_in = apply(filtering_table[,6:11],2,sum)
  tot_gen = dim(filtering_table)[1]
  plot(c(0,0.1,0.5,1,2,5),tot_gen-sum_in,xlab="cutoff",ylab="# removed genes",main=Sp)
  abline(v=cutoff,lty=2)
  
  rpk_matrix_in = rpk_matrix[row.names(rpk_matrix)%in%row.names(filtering_table_in),]
  write.csv(rpk_matrix_in,paste("data/datamatrix_filteredbycutoff",cutoff,"_",Sp,".csv",sep=""))
  print(dim(rpk_matrix_in))
}


# (2) GeTMM normalization ######################################################
# installation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

# read package
library(edgeR)
#edgeRUsersGuide()

source("scripts/my_tools.R")

Sp_ls = c("Le","Lg","Qg","Qa")
for(i in 1:length(Sp_ls)){
  Sp = Sp_ls[i]
  
  # read dataset (RPK)
  rpk_data = read.csv(paste(paste("data/datamatrix_filteredbycutoff1_",Sp,".csv",sep="")),row.names=1)
  rpk_data_withoutNA = na.omit(rpk_data[apply(rpk_data,1,sum)!=Inf,])
  if(dim(rpk_data)[1]!=dim(rpk_data_withoutNA)[1]){
    print(paste("input data has NA values",(dim(rpk_data)[1]-dim(rpk_data_withoutNA)[1])))
  }
  
  # edgeR input
  group <- c(rep("A",ncol(rpk_data_withoutNA)))
  rpk.norm <- DGEList(counts=rpk_data_withoutNA,group=group)
  rpk.norm <- calcNormFactors(rpk.norm)
  norm.counts.rpk_edger <- cpm(rpk.norm)
  
  # check
  head(norm.counts.rpk_edger)
  # save
  write.csv(norm.counts.rpk_edger,paste("data/datamatrix_getmm_",Sp,".csv",sep=""))
}

# (3) Extraction of single-copy ortholog within chromosomes ####################
# read data
one2one_gene_list = read.csv("tables/one2one_OG_list_only_chr.txt",sep="\t")
ann_Lthed = read.csv("data/annotation_Lthed_r1.0.1_HC.csv")
ann_Qgl = read.csv("data/annotation_Qgl_r1.0.1_HC.csv")

print(paste("Number of 1:1 genes =",dim(one2one_gene_list)[1]))

Sp_ls = c("Le","Lg","Qg","Qa")
for(i in 1:length(Sp_ls)){
  Sp = Sp_ls[i]
  datamatrix = read.csv(paste("data/datamatrix_getmm_",Sp,".csv",sep=""),row.names=1)
  if(substr(Sp,1,1)=="L"){
    ann = ann_Lthed
    OG_gene_list = one2one_gene_list[c("OG","L.edulis")]
  } else {
    ann = ann_Qgl
    OG_gene_list = one2one_gene_list[c("OG","Q.glauca")]
  }
  
  # extract submatrix
  datamatrix_1to1 = datamatrix[row.names(datamatrix)%in%OG_gene_list[,2],]
  
  # reorder OGID and add OGID to the matrix
  OG_gene_list_in = OG_gene_list[OG_gene_list[,2]%in%row.names(datamatrix_1to1),]
  OG_gene_list_sorted = OG_gene_list_in[order(factor(OG_gene_list_in[,2],levels=row.names(datamatrix_1to1 ))),]
  datamatrix_1to1_with_OG = data.frame(datamatrix_1to1)
  datamatrix_1to1_with_OG["OG"] = OG_gene_list_sorted["OG"]
  # save separate file
  write.csv(datamatrix_1to1_with_OG,paste("data/datamatrix_getmm_1to1_",Sp,".csv",sep=""))
  
  # other genes
  datamatrix_not1to1 = datamatrix[row.names(datamatrix)%not.in%OG_gene_list[,2],]
  write.csv(datamatrix_not1to1,paste("data/datamatrix_getmm_not1to1_",Sp,".csv",sep=""))
  
  print(paste("Number of 1:1 genes in",Sp,"datamatrix =",dim(datamatrix_1to1)[1],";others =",dim(datamatrix_not1to1)[1]))
}

# extract common genes
data_Le = read.csv("data/datamatrix_getmm_1to1_Le.csv",row.names=1)
data_Lg = read.csv("data/datamatrix_getmm_1to1_Lg.csv",row.names=1)
data_Qg = read.csv("data/datamatrix_getmm_1to1_Qg.csv",row.names=1)
data_Qa = read.csv("data/datamatrix_getmm_1to1_Qa.csv",row.names=1)

OG_Le = data_Le$OG
OG_Lg = data_Lg$OG
OG_Qg = data_Qg$OG
OG_Qa = data_Qa$OG

OG_LeLg = intersect(OG_Le,OG_Lg)
OG_QgQa = intersect(OG_Qg,OG_Qa)
OG_LeLgQgQa = intersect(OG_LeLg,OG_QgQa)

length(OG_LeLgQgQa)

data_common_Le = data_Le[data_Le$OG%in%OG_LeLgQgQa,]
data_common_Lg = data_Lg[data_Lg$OG%in%OG_LeLgQgQa,]
data_common_Qg = data_Qg[data_Qg$OG%in%OG_LeLgQgQa,]
data_common_Qa = data_Qa[data_Qa$OG%in%OG_LeLgQgQa,]

dim(data_common_Le)
dim(data_common_Lg)
dim(data_common_Qg)
dim(data_common_Qa)

# reorder
data_common_sort_Le = data_common_Le[order(data_common_Le$OG),]
data_common_sort_Lg = data_common_Lg[order(data_common_Lg$OG),]
data_common_sort_Qg = data_common_Qg[order(data_common_Qg$OG),]
data_common_sort_Qa = data_common_Qa[order(data_common_Qa$OG),]

head(data_common_sort_Le$OG)
head(data_common_sort_Lg$OG)
head(data_common_sort_Qg$OG)
head(data_common_sort_Qa$OG)

# integrate
n_Le = dim(data_common_sort_Le)[2] - 1
n_Lg = dim(data_common_sort_Lg)[2] - 1
n_Qg = dim(data_common_sort_Qg)[2] - 1
n_Qa = dim(data_common_sort_Qa)[2] - 1

data_common_LeLg = cbind(data_common_sort_Le[,1:n_Le],data_common_sort_Lg[,1:n_Lg])
data_common_QgQa = cbind(data_common_sort_Qg[,1:n_Qg],data_common_sort_Qa[,1:n_Qa])
data_common_LeLgQaQg = cbind(data_common_LeLg,data_common_QgQa)

row.names(data_common_LeLgQaQg) = data_common_sort_Le$OG
gene_Le = row.names(data_common_LeLg)
gene_Qg = row.names(data_common_QgQa)
gene_label_table = data.frame(L.edulis=gene_Le,
                              Q.glauca=gene_Qg)
row.names(gene_label_table) = row.names(data_common_LeLgQaQg)
data_intgrated = cbind(gene_label_table,data_common_LeLgQaQg)
write.csv(data_intgrated,"data/datamatrix_getmm_1to1_common_LeLgQgQa.csv")




