# R script for hierarchical clustering of gene expression profiles
# Shuichi N Kudo
# Last update 2025. 2. 26
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

my_func_hc = function(x,k){
  hc = pheatmap(x,
                cluster_cols=F,cluster_rows=T,scale="none",
                clustering_distance_cols = "euclidean", 
                clustering_distance_rows = "euclidean", 
                clustering_method = "ward.D2",
                silent = T)$tree_row
  hdendro_id = cutree(hc,k=k)
  result = list(cluster=hdendro_id)
  return(result)
}

# packages #####################################################################
library(pheatmap)
library(NbClust)
library(factoextra)
library(Rmisc)

# clustering samples ###########################################################
# Input #
data_all = read.csv("data/datamatrix_getmm_1to1_common_LeLgQgQa.csv",row.names=1)
ann = read.csv("data/annotation_with_ClusterID.csv")
data = data_all[,3:dim(data_all)[2]]
data = data[,(substr(colnames(data),4,4)=="L")|(substr(colnames(data),4,4)=="B")]
head(data)

# data transformation
log_exp = log2(data+1)
mean_log_exp = mean_individual(log_exp)
x = mean_log_exp
y = t(scale(t(x)))
Nsample=dim(x)[2]
Ngene=dim(x)[1]

# determine optimal number of clusters
wss_test = fviz_nbclust(t(y),my_func_hc,method="wss",k.max=20,nboot=100,verbose=interactive())
Nsample=dim(y)[1]
Ngene=dim(y)[2]
Sp_ls = unique(substr(colnames(data),1,2))
Tis_ls = unique(substr(colnames(data),4,4))
data_label = paste("scale_","row","_dist_","euclidean","_",paste(Sp_ls,collapse = ""),"_",paste(Tis_ls,collapse = ""),sep="")
data_info = paste(Ngene,"genes, ",Nsample,"samples",sep="")
wss_data = wss_test$data
write.csv(wss_data,paste("tables/Fig2_clustering_samples_elbow_plot_",data_label,".csv",sep=""))
g1 = wss_test
g1 = g1 + ggtitle(data_label,subtitle=data_info)
g1 = g1 + theme(aspect.ratio=1,plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))
plot(g1)
ggsave(paste("figures/Fig2_clustering_samples_elbow_plot_",data_label,".png",sep=""),plot=g1,height=5,width=5)


# annotation
sample_ls = colnames(y)
species = substr(sample_ls,1,2)
tissue = substr(sample_ls,4,4)
date = as.Date(substr(sample_ls,6,11),format="%y%m%d")
month = as.numeric(substr(date,6,7))
month_name = factor(month,levels=c(3:12,1:2),labels=c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb"))
annotation_col_table = data.frame(Month=month_name,
                                  Species=species,
                                  Tissue=tissue)
row.names(annotation_col_table) = sample_ls

# phenology label
phenology_data = read.csv("data/18_datpheno_Lit.csv")
phenology_data["date"] = as.Date(paste(phenology_data$year,phenology_data$month,phenology_data$day,sep="-"))
phenology_data["new_leaf_ratio"] = phenology_data$new_leaf/phenology_data$total
annotation_col_table["date"] = as.Date(substr(row.names(annotation_col_table),6,11),format="%y%m%d")
annotation_col_table["new_leaf_ratio"] = NA

for(i in 1:dim(annotation_col_table)[1]){
  Sp_i = annotation_col_table$Species[i]
  date_i = annotation_col_table$date[i]
  annotation_col_table$new_leaf_ratio[i] = phenology_data$new_leaf_ratio[(phenology_data$species==Sp_i)&(phenology_data$date==date_i)]
}

annotation_col_table$new_leaf_ratio[is.na(annotation_col_table$new_leaf_ratio)] = 0
annotation_col_table$new_leaf_ratio[annotation_col_table$Tissue=="B"] = 0

month_palette =  colorRampPalette(c(rep(c("#b5f1c6","#fbfec9","#ca3fb9","#292f8e"),3)))(36)[13:24]
annotation_colors = list(Month=c("Mar"=month_palette[1],"Apr"=month_palette[2],"May"=month_palette[3],
                                 "Jun"=month_palette[4],"Jul"=month_palette[5],"Aug"=month_palette[6],
                                 "Sep"=month_palette[7],"Oct"=month_palette[8],"Nov"=month_palette[9],
                                 "Dec"=month_palette[10],"Jan"=month_palette[11],"Feb"=month_palette[12]),
                         Species=c("Le"="#8CB4D4","Lg"="#446A9D","Qg"="#BC7660","Qa"="#E3B692"),
                         Tissue =c("L"="#4B9A9B","B"="#D29E43"),
                         new_leaf_ratio = colorRampPalette(c("white","#b5f1c6"))(256)
)

col_labels =  paste(substr(colnames(y),6,7),substr(colnames(y),8,9),sep="/") # YY/MM
col_labels[seq(2,length(col_labels),2)] = "" # skip odd number

expression_palette = colorRampPalette(c("#2A1A6D","#3E41A8","#91C6DC","#D0EA58","#DFAD4A","#CE422F","#782721"))

# climate information
climate_info = read.csv("tables/sample_list_with_climate_information.csv",row.names=1)
climate_info_in_data = climate_info[climate_info$Tissue%in%c("L","B"),]
climate_info_in_data_sorted = climate_info_in_data[heatmap$tree_col$order,]
climate_info_in_data_sorted$SampleID = factor(climate_info_in_data_sorted$SampleID,levels=climate_info_in_data_sorted$SampleID)

ggplot(climate_info_in_data_sorted,aes(x=SampleID,y=Temperature,color=Temperature))+
  geom_point(size=0.5)+
  #scale_color_gradient2(low="dodgerblue3",mid="#E2D858",high="indianred3",midpoint=20)+
  scale_color_gradientn(colors=c("#292f8e","dodgerblue3","#E2D858","indianred3"))+
  #scale_color_gradient(low="dodgerblue3",high="indianred3")+
  theme_bw()+
  theme(aspect.ratio=0.08,
        axis.text.x = element_blank())
ggsave("figures/Fig2/temperature_label_gradient.pdf",height=2,width=6)


# heatmap
heatmap = pheatmap(y,
                   cluster_cols=T,cluster_rows=F,
                   scale="none",
                   clustering_distance_cols = "euclidean",
                   clustering_distance_rows = "euclidean",
                   clustering_method = "ward.D2",
                   color=expression_palette(256),
                   annotation_col=annotation_col_table[,c("Species","Tissue","Month","new_leaf_ratio")],
                   annotation_colors = annotation_colors,
                   cellwidth=1.75, cellheight=0.02,fontsize = 4,
                   show_rownames = FALSE, show_colnames = F,
                   breaks = seq(-5,5,length.out=256),
                   main = paste(Ngene,"genes, ",Nsample,"samples, ","scale=row, ","dist=euclidean, ","method=ward.D2",sep=""),
                   filename = "figures/Fig2/sample_clustering_euclid_wardD2_without_floweropening_ver20250121.tiff"
)

# save clustering result ##
clustering_data = heatmap$tree_col
sample_ls = colnames(y)
sample_order = clustering_data$order
cluster_id = cutree(clustering_data,5)

clustering_result = data.frame(SampleID=sample_ls[sample_order],
                               Species=substr(sample_ls[sample_order],1,2),
                               Tissue=substr(sample_ls[sample_order],4,4),
                               Year=as.numeric(paste("20",substr(sample_ls[sample_order],6,7),sep="")),
                               Month=as.numeric(substr(sample_ls[sample_order],8,9)),
                               Date=as.numeric(substr(sample_ls[sample_order],10,11)),
                               ClusterID=cluster_id[sample_order])
row.names(clustering_result) = 1:dim(clustering_result)[1]
clustering_result = clustering_result[rev(1:dim(clustering_result)[1]),]
write.csv(clustering_result,"tables/Fig2/summary_clustering_sorted.csv")


# clustering genes #############################################################
# data
data_all = read.csv("data/datamatrix_getmm_1to1_common_LeLgQgQa.csv",row.names=1)
ann_Le = read.csv("data/annotation_Lthed_r1.0.1_HC.csv")
ann_Qg = read.csv("data/annotation_Qgl_r1.0.1_HC.csv")
gene_table_in_data = data_all[,1:2]
data = data_all[,3:dim(data_all)[2]]
data = data[,(substr(colnames(data),4,4)=="L")|(substr(colnames(data),4,4)=="B")]

# data transformation
log_exp = log2(data+1)
mean_log_exp = mean_individual(log_exp)
x = mean_log_exp
y = t(scale(t(data)))

# determine optimal number of clusters
wss_test = fviz_nbclust(y,my_func_hc,method="wss",k.max=20,nboot=100,verbose=interactive())
Nsample=dim(x)[2]
Ngene=dim(x)[1]
Sp_ls = unique(substr(colnames(data),1,2))
Tis_ls = unique(substr(colnames(data),4,4))
data_label = paste("scale_",scale,"_dist_",distance,"_",paste(Sp_ls,collapse = ""),"_",paste(Tis_ls,collapse = ""),sep="")
data_info = paste(Ngene,"genes, ",Nsample,"samples",sep="")
wss_data = wss_test$data
write.csv(gap_data,paste("tables/clustering_genes_gap_plot_",data_label,".csv",sep=""))
g1 = wss_test
g1 = g1 + ggtitle(data_label,subtitle=data_info)
g1 = g1 + theme(aspect.ratio=1,plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))
plot(g1)
ggsave(paste("figures/Fig2_clustering_genes_elbow_plot_",data_label,".png",sep=""),plot=g1,height=5,width=5)


# sample reorder
sample_order = c("Qg_B","Qa_B","Lg_B","Le_B","Qg_L","Qa_L","Lg_L","Le_L")
gaps_sorted = c()
for(i in 1:length(sample_order)){
  y_i =  y[,substr(colnames(y),1,4)==sample_order[i]]
  x_i = x[,substr(colnames(x),1,4)==sample_order[i]]
  if(i==1){
    y_sorted = y_i
    x_sorted = x_i
    gaps_sorted = c(gaps_sorted,dim(y_sorted)[2])
  } else if(i==length(sample_order)) {
    y_sorted = cbind(y_sorted,y_i)
    x_sorted = cbind(x_sorted,x_i)
  } else {
    y_sorted = cbind(y_sorted,y_i)
    x_sorted = cbind(x_sorted,x_i)
    gaps_sorted = c(gaps_sorted,dim(y_sorted)[2])
  }
}

sample_ls = colnames(y_sorted)
species = substr(sample_ls,1,2)
tissue = substr(sample_ls,4,4)
date = as.Date(substr(sample_ls,6,11),format="%y%m%d")
month = as.numeric(substr(date,6,7))
month_name = factor(month,levels=c(3:12,1:2),labels=c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb"))

annotation_col_table = data.frame(Month=month_name,
                                  Species=species,
                                  Tissue=tissue)
row.names(annotation_col_table) = sample_ls
collabs =  paste(substr(colnames(y_sorted),6,7),substr(colnames(y_sorted),8,9),sep="/")
collabs[seq(2,length(collabs),2)] = ""

# phenology label 
phenology_data = read.csv("data/18_datpheno_Lit.csv")
phenology_data["date"] = as.Date(paste(phenology_data$year,phenology_data$month,phenology_data$day,sep="-"))
phenology_data["floweropening_ratio"] = phenology_data$floweropening/phenology_data$total
phenology_data["new_leaf_ratio"] = phenology_data$new_leaf/phenology_data$total

annotation_col_table["date"] = as.Date(substr(row.names(annotation_col_table),6,11),format="%y%m%d")
annotation_col_table["floweropening_ratio"] = NA
annotation_col_table["new_leaf_ratio"] = NA

for(i in 1:dim(annotation_col_table)[1]){
  Sp_i = annotation_col_table$Species[i]
  date_i = annotation_col_table$date[i]
  annotation_col_table$floweropening_ratio[i] = phenology_data$floweropening_ratio[(phenology_data$species==Sp_i)&(phenology_data$date==date_i)]
  annotation_col_table$new_leaf_ratio[i] = phenology_data$new_leaf_ratio[(phenology_data$species==Sp_i)&(phenology_data$date==date_i)]
}

annotation_col_table$floweropening_ratio[is.na(annotation_col_table$floweropening_ratio)] = 0
annotation_col_table$new_leaf_ratio[is.na(annotation_col_table$new_leaf_ratio)] = 0

month_palette =  colorRampPalette(c(rep(c("#b5f1c6","#fbfec9","#ca3fb9","#292f8e"),3)))(36)[13:24]
annotation_colors = list(Month=c("Mar"=month_palette[1],"Apr"=month_palette[2],"May"=month_palette[3],
                                 "Jun"=month_palette[4],"Jul"=month_palette[5],"Aug"=month_palette[6],
                                 "Sep"=month_palette[7],"Oct"=month_palette[8],"Nov"=month_palette[9],
                                 "Dec"=month_palette[10],"Jan"=month_palette[11],"Feb"=month_palette[12]),
                         Species=c("Le"="#8CB4D4","Lg"="#446A9D","Qg"="#BC7660","Qa"="#E3B692"),
                         Tissue =c("L"="#355792","B"="#D29E43"),
                         floweropening_ratio = colorRampPalette(c("white","#ca3fb9"))(256),
                         new_leaf_ratio = colorRampPalette(c("white","aquamarine3"))(256)
)

col_labels =  paste(substr(colnames(y),6,7),substr(colnames(y),8,9),sep="/") # YY/MM
col_labels[seq(2,length(col_labels),2)] = "" # skip odd number

expression_palette = colorRampPalette(c("#2A1A6D","#3E41A8","#91C6DC","#D0EA58","#DFAD4A","#CE422F","#782721"))


# heatmap
heatmap = pheatmap(y_sorted,
                   cluster_cols=F,cluster_rows=T,
                   scale="none",
                   clustering_distance_cols = "euclidean",
                   clustering_distance_rows = "euclidean",
                   clustering_method = "ward.D2",
                   color=expression_palette(256),
                   annotation_col=annotation_col_table[,c("Species","Tissue","Month","new_leaf_ratio","floweropening_ratio")],
                   annotation_colors = annotation_colors,
                   gaps_col = gaps_sorted,
                   cutree_row = 5,
                   cellwidth=4, cellheight=0.02,fontsize = 4,
                   labels_col = collabs,
                   show_rownames = FALSE, show_colnames = T,
                   breaks = seq(-2,3,length.out=256),
                   main = paste(Ngene,"genes, ",Nsample,"samples, ","scale=row, ","dist=euclidean, ","method=ward.D2",sep=""),
                   filename = "figures/Fig2/gene_clustering_euclid_wardD2.tiff"
)

# calculate the number of species in which the significant periodicity was detected 
Tis_ls = c("L","B")
for(i in 1:length(Tis_ls)){
  
  Tis = Tis_ls[i]
  
  rain_Le = read.csv(paste("tables/wave_analysis/rain_result_mean_month_delta11_processed_","Le","_",Tis,".csv",sep=""),row.names=1)
  rain_Lg = read.csv(paste("tables/wave_analysis/rain_result_mean_month_delta11_processed_","Lg","_",Tis,".csv",sep=""),row.names=1)
  rain_Qg = read.csv(paste("tables/wave_analysis/rain_result_mean_month_delta11_processed_","Qg","_",Tis,".csv",sep=""),row.names=1)
  rain_Qa = read.csv(paste("tables/wave_analysis/rain_result_mean_month_delta11_processed_","Qa","_",Tis,".csv",sep=""),row.names=1)
  
  # long -> aperiodic
  rain_Le$periodicity[rain_Le$period_type_wide=="long"] = "aperiodic"
  rain_Lg$periodicity[rain_Lg$period_type_wide=="long"] = "aperiodic"
  rain_Qg$periodicity[rain_Qg$period_type_wide=="long"] = "aperiodic"
  rain_Qa$periodicity[rain_Qa$period_type_wide=="long"] = "aperiodic"
  
  periodicity_matrix = data.frame(Le=rain_Le$periodicity,
                                  Lg=rain_Lg$periodicity,
                                  Qg=rain_Qg$periodicity,
                                  Qa=rain_Qa$periodicity)
  
  # count
  period_number_ls = numeric(dim(periodicity_matrix)[1])
  for(k in 1:dim(periodicity_matrix)[1]){
    period_number = sum(periodicity_matrix[k,]=="periodic")
    period_number_ls[k] = period_number 
  }
  
  #hist(period_number_ls)
  annotation_row = data.frame(period_number=period_number_ls)
  row.names(annotation_row) = row.names(rain_Le)
  
  annotation_colors = list(Month=c("Mar"=month_palette[1],"Apr"=month_palette[2],"May"=month_palette[3],
                                   "Jun"=month_palette[4],"Jul"=month_palette[5],"Aug"=month_palette[6],
                                   "Sep"=month_palette[7],"Oct"=month_palette[8],"Nov"=month_palette[9],
                                   "Dec"=month_palette[10],"Jan"=month_palette[11],"Feb"=month_palette[12]),
                           Species=c("Le"="#8CB4D4","Lg"="#446A9D","Qg"="#BC7660","Qa"="#E3B692"),
                           Tissue =c("L"="#355792","B"="#D29E43"),
                           floweropening_ratio = colorRampPalette(c("white","#ca3fb9"))(256),
                           new_leaf_ratio = colorRampPalette(c("white","aquamarine3"))(256),
                           period_number =  colorRampPalette(c("white","deepskyblue2"))(256)
  )
  
  # heatmap
  heatmap = pheatmap(y_sorted,
                     cluster_cols=F,cluster_rows=T,
                     scale="none",
                     clustering_distance_cols = "euclidean",
                     clustering_distance_rows = "euclidean",
                     clustering_method = "ward.D2",
                     color=expression_palette(256),
                     annotation_col=annotation_col_table[,c("Species","Tissue","Month","new_leaf_ratio","floweropening_ratio")],
                     annotation_row=annotation_row,
                     annotation_colors = annotation_colors,
                     gaps_col = gaps_sorted,
                     cutree_row = 5,
                     cellwidth=4, cellheight=0.02,fontsize = 4,
                     labels_col = collabs,
                     show_rownames = FALSE, show_colnames = T,
                     breaks = seq(-2,3,length.out=256),
                     main = paste(Ngene,"genes, ",Nsample,"samples, ","scale=row, ","dist=euclidean, ","method=ward.D2",sep=""),
                     filename = paste("figures/Fig2/gene_clustering_euclid_wardD2_with_rain_",Tis,".tiff",sep="")
  )
  
  
}

# save
ClusterID_raw = cutree(heatmap$tree_row,5)
Order_clustering = heatmap$tree_row$order
OrthogroupID = heatmap$tree_row$labels

y_sorted_sorted = y_sorted[Order_clustering,]
ClusterID_sorted = ClusterID_raw[Order_clustering]
Cluster_list = unique(ClusterID_sorted)
Cluster_list_new = c("I","II","III","IV","V")
ClusterID_new = ClusterID_raw
for(i in 1:length(Cluster_list)){
  ClusterID_new[ClusterID_raw==Cluster_list[i]] = Cluster_list_new[i]
}

# gene set
summary_gene_cluster = data.frame(OrthogroupID = row.names(x_sorted),
                                  GeneID.L.edulis = ann$GeneID.L.edulis,
                                  GeneID.Q.glauca = ann$GeneID.Q.glauca,
                                  GeneID.A.thaliana_L.edulis = ann$AtID.L.edulis,
                                  GeneID.A.thaliana_Q.glauca = ann$AtID.Q.glauca,
                                  ClusterID = ClusterID_new)

summary_gene_cluster = cbind(summary_gene_cluster,x_sorted)
write.csv(summary_gene_cluster,"tables/Fig2/gene_clustering_sortbyOrthogroupID.csv",row.names=F)


# Comparison of climatic factors across clusters ###############################
# Input
climate_infomation = read.csv("tables/Fig2/climate_information.csv",row.names=1)
clustering_information = read.csv("tables/Fig2/sample_set_summary.csv")
#all(clustering_information$SampleID==climate_infomation$SampleID) # --> TRUE
climate_infomation["Cluster.Symbol"] = clustering_information$Cluster.Symbol

# Statistical test and plot
test_set = c("A-B","C-D-E")
climate_set = c("Temperature","Cumulative7","Cumulative14","Cumulative21","Photoperiod")

for(i in 1:length(test_set)){
  set = str_split(test_set[i],"-")[[1]]
  
  for(j in 1:length(climate_set)){
    clm =  climate_set[j]
    
    # extract subset
    climate_infomation_sub = climate_infomation[climate_infomation$Cluster.Symbol%in%set,]
    comp_df = data.frame(SampleID=climate_infomation_sub$SampleID,
                         ClusterID=climate_infomation_sub$Cluster.Symbol,
                         ClimateValue=climate_infomation_sub[[clm]])
    N_total = dim(comp_df)[1]
    
    # remove NA
    comp_df_rmNA = na.omit(comp_df)
    N_NA = N_total - dim(comp_df_rmNA)[1]
    
    # statistical test
    if(length(set)==2){
      # Mann-Whitney U test
      res = wilcox.test(x=comp_df_rmNA$ClimateValue[comp_df_rmNA$ClusterID==set[1]],y=comp_df_rmNA$ClimateValue[comp_df_rmNA$ClusterID==set[2]],exact=F)
      res_table = data.frame(ClimateVariable=clm,
                             Clusters=test_set[i],
                             Method = "Mann-Whitney U test",
                             kr_test_chi_sq = NA,
                             kr_test_p = NA,
                             p.value=res$p.value,
                             W.stat = res$statistic
      )
      
    } else {
      # Kruskal Wallis test
      kr = kruskal.test(comp_df_rmNA$ClimateValue ~ as.factor(comp_df_rmNA$ClusterID))
      
      # Steel Dwass test
      res = pSDCFlig(comp_df_rmNA$ClimateValue, as.factor(comp_df_rmNA$ClusterID), method="Asymptotic")
      res_table = data.frame(ClimateVariable=clm,
                             Clusters=res$labels,
                             Method = "Steel Dwass test",
                             kr_test_chi_sq = kr$statistic,
                             kr_test_p = kr$p.value,
                             p.value=res$p.val,
                             W.stat = res$obs.stat
      )
    }
    
    if((i==1)&(j==1)){
      res_table_all = res_table
    } else {
      res_table_all = rbind(res_table_all,res_table)
    }
    
  }
}

write.csv(res_table_all,"tables/Fig2/Statistical_test_for_climate_vs_cluster.csv",row.names=F)

# Statistics
cluster_set = c("A","B","C","D","E")
climate_set = c("Temperature","Cumulative7","Cumulative14","Cumulative21","Photoperiod")

for(i in 1:length(climate_set)){
  for(j in 1:length(cluster_set)){
    
    cluster = cluster_set[j]
    climate = climate_set[i]
    
    X = climate_infomation[[climate]][climate_infomation$Cluster.Symbol==cluster]
    N_total = length(X)
    X_rmNA = na.omit(X)
    N_NA = N_total - length(X_rmNA)
    
    stat_table = data.frame(ClimateVariable=climate,
                            Clusters=cluster,
                            N_total=N_total,
                            N_NA=N_NA,
                            Mean=mean(X_rmNA),
                            SD=sd(X_rmNA),
                            Median=median(X_rmNA))
    
    if((i==1)&(j==1)){
      stat_table_all = stat_table
    } else {
      stat_table_all = rbind(stat_table_all,stat_table)
    }
  }
}

write.csv(stat_table_all,"tables/Fig2/statistics_summary_climate_vs_cluster.csv",row.names=F)

# plot
climate_set = c("Temperature","Cumulative7","Cumulative14","Cumulative21","Photoperiod")

color_species = c("Le"="#8CB4D4","Lg"="#446A9D","Qg"="#BC7660","Qa"="#E3B692")

for(i in 1:length(climate_set)){
  climate = climate_set[i]
  plot_data = data.frame(SampleID=climate_infomation$SampleID,
                         Species=climate_infomation$Species,
                         Tissue=climate_infomation$Tissue,
                         Climate=climate_infomation[[climate]],
                         Cluster=climate_infomation$Cluster.Symbol)
  
  # all
  plot_data["height"] = 0
  plot_data$height[plot_data$Cluster=="A"] = 5+1
  plot_data$height[plot_data$Cluster=="B"] = 4+1
  plot_data$height[plot_data$Cluster=="C"] = 3
  plot_data$height[plot_data$Cluster=="D"] = 2
  plot_data$height[plot_data$Cluster=="E"] = 1
  
  ggplot()+
    geom_boxplot(data=plot_data,aes(x=height,y=Climate,group=Cluster),width=0.25)+
    geom_jitter(data=plot_data,aes(x=height,y=Climate,shape=Species,color=Species,group=Cluster),width=0.2,size=1.75)+
    scale_color_manual(values=c("Le"="#8CB4D4","Lg"="#446A9D","Qg"="#BC7660","Qa"="#E3B692"))+
    theme_classic()+
    theme(aspect.ratio=2,axis.text.y = element_blank())+
    labs(x="Cluster",y=climate)+
    coord_flip()
  ggsave(paste("figures/Fig2/climate_vs_cluster_",climate,"_ver2.pdf",sep=""),height=6,width=3)
}






