# R script for calculating pairwise Pearson correlation between the time series of gene expression from two species
# Shuichi N Kudo
# Last update 2025. 9. 1
# R version 4.3.1
#
# functions ####################################################################
mean_individual=function(GE_matrix){
  n_genes = dim(GE_matrix)[1]
  sample_names = colnames(GE_matrix)
  common_names = unique(substr(sample_names,1,11))
  df_temp = as.data.frame(GE_matrix)
  ave_mat = matrix(NA,ncol=length(common_names),nrow=n_genes)
  colnames(ave_mat) = common_names
  row.names(ave_mat) = row.names(df_temp)
  for(i in 1:length(common_names)){
    df_temp_i = df_temp[substr(colnames(df_temp),1,11)==common_names[i]]
    ave_mat[,i] = rowMeans(df_temp_i)
  }
  return(ave_mat)
}


# plot distribution of pairwise correlation between orthologues ################
# data
data_all = read.csv("data/datamatrix_getmm_1to1_common_LeLgQgQa.csv",row.names=1)
ann_Le = read.csv("data/annotation_Lthed_r1.0.1_HC.csv")
ann_Qg = read.csv("data/annotation_Qgl_r1.0.1_HC.csv")
data = data_all[,3:dim(data_all)[2]]
log_exp = log2(data+1)
mean_log_exp = mean_individual(log_exp)
x = mean_log_exp

# target
species_pair_ls = c("Le-Lg","Qg-Qa","Le-Qg","Le-Qa","Lg-Qg","Lg-Qa")
Tis = "L"

# perpare correlation matrix
cor_mat = data.frame(matrix(NA,nrow=dim(x)[1],ncol=length(species_pair_ls)))
colnames(cor_mat) = species_pair_ls
row.names(cor_mat) = row.names(x)

# calculation
for(i in 1:length(species_pair_ls)){
  pair = str_split(species_pair_ls[i],"-")[[1]]
  pair_label = paste(pair[1],"_",pair[2],sep="")
  
  # read a sample matching table between pair of species
  sample_table_12 = read.csv(paste("tables/sample_list_aligned_",pair_label,".csv",sep=""))
  sample_vec_align_1 = paste(substr(sample_table_12[,1],8,9),substr(sample_table_12[,1],11,11),substr(sample_table_12[,1],1,6),sep="_")
  sample_vec_align_2 = paste(substr(sample_table_12[,2],8,9),substr(sample_table_12[,2],11,11),substr(sample_table_12[,2],1,6),sep="_")
  
  # extract target data
  x1 = x[,(substr(colnames(x),1,2)==pair[1])&(substr(colnames(x),4,4)==Tis)]
  x2 = x[,(substr(colnames(x),1,2)==pair[2])&(substr(colnames(x),4,4)==Tis)]
  x1_aligned = x1[,colnames(x1)%in%sample_vec_align_1]
  x2_aligned = x2[,colnames(x2)%in%sample_vec_align_2]
  x1_aligned_ordered = x1_aligned[,order(factor(colnames(x1_aligned),levels=sample_vec_align_1))]
  x2_aligned_ordered = x2_aligned[,order(factor(colnames(x2_aligned),levels=sample_vec_align_2))]
  
  # z-scoring (note that 0 in all sample -> 0)
  x1_non0 = x1_aligned_ordered[apply(x1_aligned_ordered,1,sum)>0,]
  x2_non0 = x2_aligned_ordered[apply(x2_aligned_ordered,1,sum)>0,]
  x1_scaled = x1_non0 %>% t() %>% scale() %>% t()
  x2_scaled = x2_non0 %>% t() %>% scale() %>% t()
  x1_scaled_all = x1_aligned_ordered
  x1_scaled_all[apply(x1_aligned_ordered,1,sum)>0,] = x1_scaled
  x2_scaled_all = x2_aligned_ordered
  x2_scaled_all[apply(x2_aligned_ordered,1,sum)>0,] = x2_scaled
  
  y1 = x1_scaled_all
  y2 = x2_scaled_all
  
  # calculate correlation
  for(j in 1:dim(cor_mat)[1]){
    cor_mat[j,i] = cor(y1[j,],y2[j,],method="pearson")
  }
}

# prepare plot data frame
cor_mat$OG <- row.names(cor_mat)
cor_mat_plot <- melt(cor_mat, id.vars = "OG", variable.name = "Species", value.name = "Value")
n_NA = dim(cor_mat_plot[is.na(cor_mat_plot[,3]),])[1]

# tukey test
res <- lm(Value~Species, d=cor_mat_plot)
tukey_res <- glht(res, linfct=mcp(Species="Tukey"))
summary(tukey_res)
mltv = cld(tukey_res, decreasing=F)
annos = mltv[["mcletters"]][["Letters"]]

# plot
ggplot(cor_mat_plot, aes(x = Species, y = Value))+
  geom_violin(aes(fill = Species),trim = FALSE)+
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +
  theme_classic()+
  theme(legend.position = "none") +
  stat_summary(aes(x= Species,
                   y=Value), geom = 'text', label =annos,  vjust = -7, hjust=2)
ggsave(paste("figures/pairwise_distance_correlation_NA",n_NA ,"_",Tis,".png",sep=""),width=5,height=4)

# save
write.csv(cor_mat,paste("correlation_matrix_",Tis,".csv",sep=""))
