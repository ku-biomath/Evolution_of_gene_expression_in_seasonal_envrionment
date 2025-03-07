# R script for species comparison of correlation, difference, and peak month
# Shuichi N Kudo
# Last update 2025. 3. 7
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

my_GO_analysis = function(gsc,gene_ls,all.Gene){
  
  # parameter for GOstats
  p <- GSEAGOHyperGParams(
    name = "Paramaters",
    geneSetCollection = gsc,
    geneIds = gene_ls,
    universeGeneIds = all.Gene,
    ontology = "BP",
    pvalueCutoff = 0.05,
    conditional = TRUE,
    testDirection = "over"
  )
  
  # Fisher's exact test
  result <- hyperGTest(p)
  
  # summarize
  raw_pval = pvalues(result)
  adj_pval = p.adjust(raw_pval, method="BH")
  odds = oddsRatios(result)
  expected_count = expectedCounts(result)
  gene_count = geneCounts(result)
  gene_univ_count = universeCounts(result)
  
  summary_table = data.frame(rank=1:length(raw_pval),
                             go_id=names(raw_pval),
                             raw.pvalues=raw_pval,
                             adjusted.pvalues=adj_pval,
                             oddsRatios=odds,
                             expectedCounts=expected_count,
                             geneCounts=gene_count,
                             universeCounts=gene_univ_count,
                             siginificance="")
  
  
  summary_table$siginificance[summary_table$adjusted.pvalues<0.001] = "***"
  summary_table$siginificance[summary_table$adjusted.pvalues<0.01] = "**"
  summary_table$siginificance[(summary_table$adjusted.pvalues>=0.01)&(summary_table$adjusted.pvalues<0.05)] = "*"
  summary_table$Genes = ""
  
  target_GOdatabase = GOdatabase[GOdatabase$Orthogroup_ID%in%gene_ls,]
  
  for(k in 1:dim(summary_table)[1]){
    GOid_k = summary_table$go_id[k]
    OGs = unique(target_GOdatabase$Orthogroup_ID[target_GOdatabase$go_id==GOid_k])
    summary_table$Genes[k] = paste(OGs, collapse = ", ")
  }
  
  row.names(summary_table) = 1:dim(summary_table)[1]
  summary_table = left_join(summary_table,go_termname,by="go_id")
  
  return(summary_table)
}




# library ######################################################################
library("GOstats")
library("GSEABase")
library("qvalue")
library(dplyr)

# Pearson's correlation of seasonal pattern across species #####################
data_all = read.csv("data/datamatrix_getmm_1to1_common_LeLgQgQa.csv",row.names=1)
data = data_all[,3:dim(data_all)[2]]
log_exp = log2(data+1)
mean_log_exp = mean_individual(log_exp)
x = mean_log_exp

# target
species_pair_ls = c("Le_Lg","Qg_Qa","Le_Qg","Le_Qa","Lg_Qg","Lg_Qa")
Tis_ls = c("L","B")

for(j in 1:length(Tis_ls)){
  Tis = Tis_ls[j]
  # perpare correlation matrix
  cor_mat = data.frame(matrix(NA,nrow=dim(x)[1],ncol=length(species_pair_ls)))
  colnames(cor_mat) = species_pair_ls
  row.names(cor_mat) = row.names(x)
  
  # calculation
  for(i in 1:length(species_pair_ls)){
    pair = str_split(species_pair_ls[i],"_")[[1]]
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
    head(x1_aligned_ordered)
    head(x2_aligned_ordered)
    
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
    for(k in 1:dim(cor_mat)[1]){
      cor_mat[k,i] = cor(y1[k,],y2[k,],method="pearson")
    }
  }
  write.csv(cor_mat,paste("tables/pairwise_correlation_",Tis,".csv",sep=""))
}

# GO analysis on the genes with highest / lowest mean correlation ##############

# gene extraction
Tis_ls = c("L","B")
ann = read.csv("data/annotation_with_ClusterID.csv")

for(i in 1:length(Tis_ls)){
  Tis = Tis_ls[i]
  
  # correlation matrix
  cor_mat = read.csv(paste("tables/Fig3/pairwise_correlation_",Tis,".csv",sep=""),row.names=1)
  
  # gene list
  ann_info = ann[c("OrthogroupID","GeneID.L.edulis","GeneID.Q.glauca","AtID.L.edulis","AtID.Q.glauca","AtSymbol.L.edulis","AtSymbol.Q.glauca","Description.L.edulis","Description.Q.glauca","Cluster")]
  
  # check different AtID in gene list
  ann_info["AtID"] = ann_info$AtID.L.edulis
  ann_info$AtID[((is.na(ann_info$AtID.L.edulis)==F)&(is.na(ann_info$AtID.Q.glauca)==F))&(str_sub(ann_info$AtID.L.edulis,1,9)!=str_sub(ann_info$AtID.Q.glauca,1,9))] = paste(ann_info$AtID.L.edulis[((is.na(ann_info$AtID.L.edulis)==F)&(is.na(ann_info$AtID.Q.glauca)==F))&(str_sub(ann_info$AtID.L.edulis,1,9)!=str_sub(ann_info$AtID.Q.glauca,1,9))],
                                                                                                                                                                            ann_info$AtID.Q.glauca[((is.na(ann_info$AtID.L.edulis)==F)&(is.na(ann_info$AtID.Q.glauca)==F))&(str_sub(ann_info$AtID.L.edulis,1,9)!=str_sub(ann_info$AtID.Q.glauca,1,9))],sep=";")
  ann_info["Description"] = ann_info$Description.L.edulis
  ann_info$Description[((is.na(ann_info$AtID.L.edulis)==F)&(is.na(ann_info$AtID.Q.glauca)==F))&(str_sub(ann_info$AtID.L.edulis,1,9)!=str_sub(ann_info$AtID.Q.glauca,1,9))] = paste(ann_info$Description.L.edulis[((is.na(ann_info$AtID.L.edulis)==F)&(is.na(ann_info$AtID.Q.glauca)==F))&(str_sub(ann_info$AtID.L.edulis,1,9)!=str_sub(ann_info$AtID.Q.glauca,1,9))],
                                                                                                                                                                                   ann_info$Description.Q.glauca[((is.na(ann_info$AtID.L.edulis)==F)&(is.na(ann_info$AtID.Q.glauca)==F))&(str_sub(ann_info$AtID.L.edulis,1,9)!=str_sub(ann_info$AtID.Q.glauca,1,9))],sep=";")
  ann_info["AtSymbol"] = ann_info$AtSymbol.L.edulis
  ann_info$AtSymbol[((is.na(ann_info$AtID.L.edulis)==F)&(is.na(ann_info$AtID.Q.glauca)==F))&(str_sub(ann_info$AtID.L.edulis,1,9)!=str_sub(ann_info$AtID.Q.glauca,1,9))] = paste(ann_info$AtSymbol.L.edulis[((is.na(ann_info$AtID.L.edulis)==F)&(is.na(ann_info$AtID.Q.glauca)==F))&(str_sub(ann_info$AtID.L.edulis,1,9)!=str_sub(ann_info$AtID.Q.glauca,1,9))],
                                                                                                                                                                                ann_info$AtSymbol.Q.glauca[((is.na(ann_info$AtID.L.edulis)==F)&(is.na(ann_info$AtID.Q.glauca)==F))&(str_sub(ann_info$AtID.L.edulis,1,9)!=str_sub(ann_info$AtID.Q.glauca,1,9))],sep=";")
  ann_info$AtID[(is.na(ann_info$AtID.L.edulis)==T)&(is.na(ann_info$AtID.Q.glauca)==F)] = ann_info$AtID.Q.glauca[(is.na(ann_info$AtID.L.edulis)==T)&(is.na(ann_info$AtID.Q.glauca)==F)]
  ann_info$Description[(is.na(ann_info$AtID.L.edulis)==T)&(is.na(ann_info$AtID.Q.glauca)==F)] = ann_info$Description.Q.glauca[(is.na(ann_info$AtID.L.edulis)==T)&(is.na(ann_info$AtID.Q.glauca)==F)]
  ann_info$AtSymbol[(is.na(ann_info$AtID.L.edulis)==T)&(is.na(ann_info$AtID.Q.glauca)==F)] = ann_info$AtSymbol.Q.glauca[(is.na(ann_info$AtID.L.edulis)==T)&(is.na(ann_info$AtID.Q.glauca)==F)]
  
  ann_info$AtID[(is.na(ann_info$AtID.L.edulis)==F)&(is.na(ann_info$AtID.Q.glauca)==T)] = ann_info$AtID.L.edulis[(is.na(ann_info$AtID.L.edulis)==F)&(is.na(ann_info$AtID.Q.glauca)==T)]
  ann_info$Description[(is.na(ann_info$AtID.L.edulis)==F)&(is.na(ann_info$AtID.Q.glauca)==T)] = ann_info$Description.L.edulis[(is.na(ann_info$AtID.L.edulis)==F)&(is.na(ann_info$AtID.Q.glauca)==T)]
  ann_info$AtSymbol[(is.na(ann_info$AtID.L.edulis)==F)&(is.na(ann_info$AtID.Q.glauca)==T)] = ann_info$AtSymbol.L.edulis[(is.na(ann_info$AtID.L.edulis)==F)&(is.na(ann_info$AtID.Q.glauca)==T)]
  
  # complete
  ann_info_final = ann_info[c("OrthogroupID","GeneID.L.edulis","GeneID.Q.glauca","AtID","AtSymbol","Description","Cluster")]
  write.csv(ann_info_final,"tables/Gene_list_AtID_integrated_Le_Qg_20250201.csv")
  head(ann_info_final)
  
  # correlation information
  cor_info = cor_mat
  colnames(cor_info) = paste("cor_",colnames(cor_mat),sep="")
  cor_info["Ave_cor"] = apply(cor_mat,1,function(x){mean(x,na.rm=T)})
  gene_info = cbind(ann_info_final,cor_info)
  row.names(gene_info) = 1:dim(gene_info)[1]
  
  # extract genes with higher or lower correlation
  q_cor_095 = quantile(gene_info$Ave_cor,0.95,na.rm=T)
  high_cor = gene_info[(is.na(gene_info$Ave_cor)==F)&(gene_info$Ave_cor>=q_cor_095),]
  q_cor_005 = quantile(gene_info$Ave_cor,0.05,na.rm=T)
  low_cor = gene_info[(is.na(gene_info$Ave_cor)==F)&(gene_info$Ave_cor<=q_cor_005),]
  
  write.csv(high_cor,paste("tables/Fig4/gene_extract_20250203/gene_list_high_correlation_top5%_",Tis,".csv",sep=""))
  write.csv(low_cor,paste("tables/Fig4/gene_extract_20250203/gene_list_low_correlation_bottom5%_",Tis,".csv",sep=""))
  
}

# GO analysis
# preperation #
GOdatabase = read.csv("data/GOdatabase_AtID_from_PLAZA.csv",row.names=1)

# make gene set collection for GOstats
go_termname = GOdatabase  %>% dplyr::select(go_id, term) %>% unique() %>% dplyr::rename(go_Term=term)
goframeData = data.frame(go_id=GOdatabase$go_id,
                         evidence=GOdatabase$evidence,
                         OGID=GOdatabase$Orthogroup_ID)
goFrame <- GOFrame(goframeData)
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

# set background genes
all.Gene = unique(GOdatabase$OGID)

Tis_ls = c("L","B")
file_list = c("high_correlation_top5%","low_correlation_bottom5%")
for(i in 1:length(Tis_ls)){
  for(j in 1:length(file_list)){
    Tis = Tis_ls[i]
    file = file_list[j]
    
    target_gene_set = read.csv(paste("tables/gene_list_",file,"_",Tis,".csv",sep=""),row.names=1)
    
    # extract OrthogroupID
    gene_ls = target_gene_set$OrthogroupID
    print(length(gene_ls))
    
    # GO analysis
    summary_table = my_GO_analysis(gsc,gene_ls,all.Gene)
    
    # save
    write.csv(summary_table,paste("tables/",file,"_GO_analysis","_",Tis,".csv",sep=""))
    
    # create gene list
    # Gene list
    for(k in 1:dim(summary_table)[1]){
      summary_table_k = summary_table[k,]
      Genes_k = summary_table_k$Genes
      Genes_k_ls = str_split(Genes_k,", ")[[1]]
      if(nchar(Genes_k)>1){
        N_gene = length(Genes_k_ls)
        summary_table_melt_k = data.frame(Tissue=rep(Tis,N_gene),
                                          rank=rep(summary_table_k$rank,N_gene),
                                          go_id=rep(summary_table_k$go_id,N_gene),
                                          go_Term=rep(summary_table_k$go_Term,N_gene),
                                          adjusted.pvalues=rep(summary_table_k$adjusted.pvalues,N_gene),
                                          siginificance=rep(summary_table_k$siginificance,N_gene),
                                          Orthogroup=Genes_k_ls
        )
        
        ann_target_sub = target_gene_set[target_gene_set$OrthogroupID%in%summary_table_melt_k$Orthogroup,]
        if(all(ann_target_sub$OrthogroupID==summary_table_melt_k$Orthogroup)==F){
          print("error!")
        }
        ann_target_sub_required = ann_target_sub
        summary_table_melt_k = cbind(summary_table_melt_k,ann_target_sub_required)
        
        if(k==1){
          summary_table_melt = summary_table_melt_k
        } else {
          summary_table_melt = rbind(summary_table_melt,summary_table_melt_k)
        }
      }
    }
    write.csv(summary_table_melt,paste("tables/",file,"_GO_analysis_gene_list","_",Tis,".csv",sep=""))
    
    
    # create unique gene list
    
    GO_gene_list = summary_table_melt
    
    unique_gene_set = unique(GO_gene_list$Orthogroup)
    
    for(k in 1:length(unique_gene_set)){
      gene_k = unique_gene_set[k]
      GO_gene_list_k = GO_gene_list[GO_gene_list$Orthogroup==gene_k,]
      GO_gene_list_unique = GO_gene_list_k[1,]
      GO_gene_list_unique["N_go_Term"] = dim(GO_gene_list_k)[1]
      
      if(k==1){
        GO_gene_list_unique_all = GO_gene_list_unique
      } else {
        GO_gene_list_unique_all = rbind(GO_gene_list_unique_all, GO_gene_list_unique)
        
      }
    }
    
    write.csv(GO_gene_list_unique_all,paste("tables/",file,"_GO_analysis_unique_gene_list","_",Tis,".csv",sep=""))
    
  }
}


# visualization

Tis_ls = c("L","B")
file_list = c("high_correlation_top5%","low_correlation_bottom5%")
for(j in 1:length(file_list)){
  file = file_list[j]
  for(i in 1:length(Tis_ls)){
    Tis = Tis_ls[i]
    GO_result = read.csv(paste("tables/",file,"_GO_analysis","_",Tis,".csv",sep=""),row.names=1)
    
    GO_result["Ratio"] = GO_result$geneCounts/GO_result$universeCounts
    
    n_go_plot = 5
    plot_data = data.frame(rank=GO_result$rank[1:n_go_plot],
                           go_id=GO_result$go_id[1:n_go_plot],
                           adjusted.pvalues=GO_result$adjusted.pvalues[1:n_go_plot],
                           raw.pvalues=GO_result$raw.pvalues[1:n_go_plot],
                           Ratio=GO_result$Ratio[1:n_go_plot],
                           Count=GO_result$geneCounts[1:n_go_plot],
                           term=GO_result$go_Term[1:n_go_plot])
    
    # sort by ratio
    plot_data = plot_data[order(plot_data$Ratio),]
    
    plot_data$go_id = factor(plot_data$go_id,levels=(plot_data$go_id))
    
    # file name
    if(file == "high_correlation_top5%"){
      file_name = "high_correlation_top005"
    } else {
      file_name = "low_correlation_bottom005"
    }
    
    
    ggplot(plot_data, aes(x = Ratio, y = go_id, size = Count, color=-log10(adjusted.pvalues))) +
      geom_point() +
      scale_size_continuous(range = c(3, 10),limits=c(0,110)) +
      scale_color_gradient(low = "blue", high = "red",limits=c(0,8))+
      scale_x_continuous(limits=c(0,1))+
      theme_bw() +
      theme(
        aspect.ratio=2,
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)
      )
    ggsave(paste("figures/gene_extract_20250203/GO_analysis_",file_name,"_",Tis,".pdf",sep=""),height=4.5,width=4.5)
    
  }
  
}




# Multiple comparison of pairwise correlation ##################################

# Nemenyi Test
Tis_ls = c("L","B")

stats <- data.frame(
  Tissue = character(),
  SpeciesPair = character(),
  NumNA = numeric(),
  Min = numeric(),
  Mean = numeric(),
  Median = numeric(),
  Max = numeric(),
  SD = numeric(),
  stringsAsFactors = FALSE
)

Fm_results = data.frame(
  Tissue = character(),
  Chi_square = numeric(),
  df = numeric(),
  P_value = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE
)

for(i in 1:length(Tis_ls)){
  Tis = Tis_ls[i]
  cor_mat = read.csv(paste("tables/Fig3/pairwise_correlation_",Tis,".csv",sep=""),row.names=1)
  
  # remove NA
  cor_mat_rmNA = na.omit(cor_mat)
  
  # convert into long format
  cor_mat_rmNA["OrthogroupID"] = row.names(cor_mat_rmNA)
  cor_mat_long_rmNA = melt(cor_mat_rmNA,measure=species_pair_ls,id.vars = "OrthogroupID")
  head(cor_mat_long_rmNA)
  
  # statistics
  cor_mat_long_rmNA["Rank"] = rank(cor_mat_long_rmNA$value) # Rank
  
  for(j in 1:length(species_pair_ls)){
    pair = species_pair_ls[j]
    cor_mat_j = cor_mat_long_rmNA[cor_mat_long_rmNA$variable==pair,]
    cor_vec = cor_mat_j$value
    stats = rbind(stats,data.frame(
      Tissue = Tis,
      SpeciesPair=pair,
      NumNA = length(cor_mat[[pair]][is.na(cor_mat[[pair]])]),
      Min = min(cor_vec,na.rm=T),
      Mean = mean(cor_vec,na.rm=T),
      Median = median(cor_vec,na.rm=T),
      Max = max(cor_vec,na.rm=T),
      SD = sd(cor_vec,na.rm=T),
      RankSum = sum(cor_mat_j$Rank),
      RankMean = sum(cor_mat_j$Rank)/dim(cor_mat_long_rmNA)[1]
    ))
  }
  
  # Statistical test # - - - -
  # Friedman test
  fm = friedman.test(as.matrix(cor_mat_rmNA))
  chi_sq = fm$statistic
  df = fm$parameter
  p_value = fm$p.value
  signif = ""
  if(p_value<0.1){
    signif = "."
    if(p_value<0.05){
      signif = "*"
      if(p_value<0.01){
        signif = "**"
        if(p_value<0.001){
          signif = "***"
        }
      }
    }
  }
  Fm_results = rbind(Fm_results, data.frame(
    Tissue = Tis,
    Chi_square = chi_sq,
    df = df,
    P_value = p_value,
    Significance = signif,
    stringsAsFactors = FALSE
  ))
  
  # Multiple comparison (Nemenyi Test)
  res = frdAllPairsNemenyiTest(value~variable|OrthogroupID,cor_mat_long_rmNA)
  pval_mat = as.data.frame(res$p.value)
  stat_mat = as.data.frame(res$statistic)
  pval_mat["Tissue"] = Tis
  stat_mat["Tissue"] = Tis
  write.csv(pval_mat,paste("tables/Fig4/pairwise_correlation_multcomp_",Tis,"_pvalue.csv",sep=""))
  write.csv(stat_mat,paste("tables/Fig4/pairwise_correlation_multcomp_",Tis,"_stats.csv",sep=""))
}

write.csv(stats,"tables/Fig4/pairwise_correlation_statistics.csv",row.names=F)
write.csv(Fm_results,"tables/Fig4/pairwise_correlation_friedman.csv",row.names=F)

Tis_ls = c("L","B")
Sp_ls = c("Le","Lg","Qg","Qa")

for(i in 1:length(Tis_ls)){
  # Input
  Tis = Tis_ls[i]
  cor_mat = read.csv(paste("tables/pairwise_correlation_",Tis,".csv",sep=""),row.names=1)
  
  # periodicity matrix
  period_type_matrix = data.frame(matrix(NA,nrow=dim(cor_mat)[1],ncol=length(Sp_ls)))
  row.names(period_type_matrix) = row.names(cor_mat)
  colnames(period_type_matrix) = Sp_ls
  for(j in 1:length(Sp_ls)){
    Sp = Sp_ls[j]
    rain_result = read.csv(paste("tables/wave_analysis/rain_result_mean_month_delta11_processed_",Sp,"_",Tis,".csv",sep=""),row.names=1)
    period_type_matrix[(rain_result$phase_month%in%c(3,4,5,6,7,8,9,10,11))&(rain_result$period_type_wide=="annual"),j] = "growing"
    period_type_matrix[(rain_result$phase_month%in%c(12,1,2))&(rain_result$period_type_wide=="annual"),j] = "winter"
    period_type_matrix[(rain_result$period_type_wide=="short"),j] = "sixmonth"
    period_type_matrix[(rain_result$period_type_wide%in%c("aperiodic","long")),j] = "aperiodic"
  }
  
  # remove NA
  cor_mat = na.omit(cor_mat)
  period_type_matrix = period_type_matrix[row.names(period_type_matrix)%in%row.names(cor_mat),]
  
  # mean
  mean_cor = apply(cor_mat,1,mean)
  
  # sort
  period_type_matrix_sorted = period_type_matrix[rev(order(mean_cor)),]
  cor_mat_sorted = cor_mat[rev(order(mean_cor)),]
  
  # count
  period_type_count_mat = data.frame(Order=1:dim(period_type_matrix_sorted)[1],
                                     OrthogroupID=row.names(period_type_matrix_sorted),
                                     growing=0,
                                     winter=0,
                                     sixmonth=0,
                                     aperiodic=0)
  
  for(k in 1:dim(period_type_count_mat)[1]){
    period_type_k = period_type_matrix_sorted[k,]
    period_type_count_mat$growing[k] = length(period_type_k[period_type_k=="growing"])
    period_type_count_mat$winter[k] = length(period_type_k[period_type_k=="winter"])
    period_type_count_mat$sixmonth[k] = length(period_type_k[period_type_k=="sixmonth"])
    period_type_count_mat$aperiodic[k] = length(period_type_k[period_type_k=="aperiodic"])
  }
  
  # window average
  window_size = 100
  period_type_ave_mat = data.frame(Order=1:(dim(period_type_matrix_sorted)[1]-window_size),
                                   OrthogroupID=row.names(period_type_matrix_sorted)[1:(dim(period_type_matrix_sorted)[1]-window_size)],
                                   growing=0,
                                   winter=0,
                                   sixmonth=0,
                                   aperiodic=0,
                                   window=window_size)
  
  for(k in 1:dim(period_type_ave_mat)[1]){
    period_type_ave_mat[k,3:6] = apply(period_type_count_mat[k:(k+window_size-1),3:6],2,mean)
  }
  
  # covert matrix into long format to plot
  N = dim(period_type_ave_mat)[1]
  plot_data = data.frame(Order=rep(period_type_ave_mat$Order,4),
                         OrthogroupID=rep(period_type_ave_mat$OrthogroupID,4),
                         Count=c(period_type_ave_mat$growing,period_type_ave_mat$winter,period_type_ave_mat$sixmonth,period_type_ave_mat$aperiodic),
                         Type=c(rep("growing",N),rep("winter",N),rep("sixmonth",N),rep("aperiodic",N)))
  plot_data$Type = factor(plot_data$Type,levels=c("growing","winter","sixmonth","aperiodic"))
  
  my_palette_discr = c("darkolivegreen3","royalblue","gray30","gray80")
  
  ggplot(plot_data,aes(x=Order,y=Count,color=Type))+
    geom_line(lwd=0.7)+
    scale_color_manual(values=my_palette_discr)+
    scale_y_continuous(limits=c(0,4))+
    scale_x_reverse()+
    coord_flip()+
    theme_bw()+
    theme(aspect.ratio=5)
  ggsave(paste("figures/correlation_and_portion_of_growing_winter_genes_window",window_size,"_",Tis,"_ver2.pdf",sep=""),height=8,width=3)
  
  
}

# absolute difference ##########################################################
# data
data_all = read.csv("data/datamatrix_getmm_1to1_common_LeLgQgQa.csv",row.names=1)
data = data_all[,3:dim(data_all)[2]]
log_exp = log2(data+1)
mean_log_exp = mean_individual(log_exp)
x = mean_log_exp

# target
species_pair_ls = c("Le-Lg","Qg-Qa","Le-Qg","Le-Qa","Lg-Qg","Lg-Qa")
Tis_ls = c("L","B")

# calculation
for(j in 1:length(Tis_ls)){
  Tis = Tis_ls[j]
  # perpare difference matrix
  diff_mat = data.frame(matrix(NA,nrow=dim(x)[1],ncol=length(species_pair_ls)))
  colnames(diff_mat) = species_pair_ls
  row.names(diff_mat) = row.names(x)
  
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
    
    y1 = x1_aligned_ordered
    y2 = x2_aligned_ordered
    
    # calculate correlation
    for(k in 1:dim(diff_mat)[1]){
      diff_mat[k,i] = sum(abs(y1[k,]-y2[k,]))/length(y1[k,])
    }
  }
  write.csv(diff_mat,paste("tables/absolute_difference_matrix_",Tis,".csv",sep=""))
}

# Peak difference ##############################################################

data_all = read.csv("data/datamatrix_getmm_1to1_common_LeLgQgQa.csv",row.names=1)

# target
species_pair_ls = c("Le-Lg","Qg-Qa","Le-Qg","Le-Qa","Lg-Qg","Lg-Qa")
Tis_ls = c("L","B")

# calculation
for(j in 1:length(Tis_ls)){
  Tis = Tis_ls[j]
  # perpare difference matrix
  peak_diff_mat = data.frame(matrix(NA,nrow=dim(data_all)[1],ncol=length(species_pair_ls)))
  colnames(diff_mat) = species_pair_ls
  row.names(diff_mat) = row.names(data_all)
  
  for(i in 1:length(species_pair_ls)){
    pair = str_split(species_pair_ls[i],"-")[[1]]
    pair_label = paste(pair[1],"_",pair[2],sep="")
    
    # read result files of wave pattern analysis
    rain_sp1 = read.csv(paste("tables/wave_analysis/rain_result_mean_month_delta11_processed_",pair[1],"_",Tis,".csv",sep=""))
    rain_sp2 = read.csv(paste("tables/wave_analysis/rain_result_mean_month_delta11_processed_",pair[2],"_",Tis,".csv",sep=""))
    
    # peak difference
    peak_diff = circular_distance_month_set(rain_sp1$phase_month,rain_sp2$phase_month)
    peak_diff[((rain_sp1$period_type_wide!="annual")|(rain_sp2$period_type_wide!="annual"))] = NA
    
    peak_diff_mat[,i] = peak_diff
    
  }
  write.csv(peak_diff_mat,paste("tables/peak_difference_matrix_",Tis,".csv",sep=""))
}



