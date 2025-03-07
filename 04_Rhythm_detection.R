# R script for detection of genes with seasonal rhythmic expression pattern
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

# Rhythmm detection (RAIN) #####################################################
library(rain)

# Targets
Sp_ls = c("Le","Lg","Qg","Qa")
Tis_ls = c("L","B")

# expression data
data_all = read.csv("data/datamatrix_getmm_1to1_common_LeLgQgQa.csv",row.names=1)
ann_Le = read.csv("data/annotation_Lthed_r1.0.1_HC.csv")
ann_Qg = read.csv("data/annotation_Qgl_r1.0.1_HC.csv")
data = data_all[,3:dim(data_all)[2]]
log_exp = log2(data+1)
mean_log_exp = mean_individual(log_exp)
x = mean_log_exp

# RAIN
for(i in 1:length(Sp_ls)){
  for(j in 1:length(Tis_ls)){
    x_target = x[,(substr(colnames(x),1,2)==Sp)&(substr(colnames(x),4,4)==Tis)]
    starttime = Sys.time()
    rainresult = rain(t(x_target),deltat=1,period=12,period.delta=11,peak.border=c(0,1),
                      nr.series=1,method="independent",na.rm=T,adjp.method="ABH")
    
    write.csv(rainresult,paste("tables/wave_analysis/rain_result_mean_month_delta11_peak0010_",Sp,"_",Tis,".csv",sep=""))
    endtime = Sys.time()
    time_needed = endtime - starttime
    print(time_needed)
  }
}

# Classification ###############################################################
# Threshold
q_threshold = 0.01

# Targets
Sp_ls = c("Le","Lg","Qg","Qa")
Tis_ls = c("L","B")

# expression data
data_all = read.csv("data/datamatrix_getmm_1to1_common_LeLgQgQa.csv",row.names=1)
ann_Le = read.csv("data/annotation_Lthed_r1.0.1_HC.csv")
ann_Qg = read.csv("data/annotation_Qgl_r1.0.1_HC.csv")
data = data_all[,3:dim(data_all)[2]]
log_exp = log2(data+1)
mean_log_exp = mean_individual(log_exp)
x = mean_log_exp

for(i in 1:length(Sp_ls)){
  for(j in 1:length(Tis_ls)){
    Sp = Sp_ls[i]
    Tis = Tis_ls[j]
    rainresult = read.csv(paste("tables/wave_analysis/rain_result_mean_month_delta11_",Sp,"_",Tis,".csv",sep=""),row.names=1)
    rainresult["periodicity"] = "aperiodic"
    rainresult$periodicity[rainresult$pVal<q_threshold] = "periodic"
    rainresult["period_month"] = round(rainresult$period * 28 * (12/365),0) # convert timepoint n into month
    gene_names = data_all[,1:2]
    rainresult = cbind(gene_names,rainresult)
    rainresult["period_type_wide"] = "aperiodic"
    rainresult$period_type_wide[(rainresult$periodicity=="periodic")&(rainresult$period_month<8)] = "short"
    rainresult$period_type_wide[(rainresult$periodicity=="periodic")&(rainresult$period_month>=8)] = "annual"
    rainresult$period_type_wide[(rainresult$periodicity=="periodic")&(rainresult$period_month>16)] = "long"
    rainresult["period_type_sharp"] = "aperiodic"
    rainresult$period_type_sharp[(rainresult$periodicity=="periodic")&(rainresult$period_month<=4)] = "short"
    rainresult$period_type_sharp[(rainresult$periodicity=="periodic")&(rainresult$period_month>=5)&(rainresult$period_month<=7)] = "halfannual"
    rainresult$period_type_sharp[(rainresult$periodicity=="periodic")&(rainresult$period_month>=7)&(rainresult$period_month<=10)] = "intermediate1"
    rainresult$period_type_sharp[(rainresult$periodicity=="periodic")&(rainresult$period_month>=11)&(rainresult$period_month<=13)] = "annual"
    rainresult$period_type_sharp[(rainresult$periodicity=="periodic")&(rainresult$period_month>=14)&(rainresult$period_month<=22)] = "intermediate2"
    rainresult$period_type_sharp[(rainresult$periodicity=="periodic")&(rainresult$period_month>=23)] = "biannual"
    x_target = x[,(substr(colnames(x),1,2)==Sp)&(substr(colnames(x),4,4)==Tis)]
    date_ls = as.Date(substr(colnames(x_target),6,11),format="%y%m%d") 
    rainresult["phase_date"] = date_ls[1] + (rainresult$phase - 1)*28 # convert timepoint n into month
    rainresult["phase_month"] = as.numeric(substr(rainresult$phase_date,6,7))
    write.csv(rainresult,paste("tables/wave_analysis/rain_result_mean_month_delta11_processed_",Sp,"_",Tis,".csv",sep=""))
    
  }
}

# Distribution of period #######################################################
library(ggplot2)

Sp_ls = c("Le","Lg","Qg","Qa")
Tis_ls = c("L","B")

for(i in 1:length(Sp_ls)){
  for(j in 1:length(Tis_ls)){
    Sp = Sp_ls[i]
    Tis = Tis_ls[j]
    
    rain_res = read.csv(paste("tables/wave_analysis/rain_result_mean_month_delta11_processed_",Sp,"_",Tis,".csv",sep=""),row.names=1)
    
    # extract genes with rhythmic expression
    rain_res_per = rain_res[rain_res$period_type_wide%in%c("annual","short"),]
    
    # histogram
    ggplot(rain_res_per,aes(x=period_month))+
      geom_histogram(binwidth=1,color="black",fill="paleturquoise3")+
      scale_x_continuous(breaks=3:17,limits=c(3,17))+
      theme_bw()+
      theme(aspect.ratio=0.7)+
      ggtitle(paste("Species: ",Sp,", Tissue: ",Tis,", Num of periodic genes: ",dim(rain_res_per)[1],sep=""))
    ggsave(paste("figures/period_distribution_",Sp,"_",Tis,".pdf"),height=4,width=6)
    
  }
}

# Distribution of peak month ###################################################
library(ggplot2)
Sp_ls = c("Le","Lg","Qg","Qa")
Tis_ls = c("L","B")

for(i in 1:length(Sp_ls)){
  for(j in 1:length(Tis_ls)){
    Sp = Sp_ls[i]
    Tis = Tis_ls[j]
    
    rain_res = read.csv(paste("tables/wave_analysis/rain_result_mean_month_delta11_processed_",Sp,"_",Tis,".csv",sep=""),row.names=1)
    
    annual = rain_res[rain_res$period_type_wide=="annual",]
    
    plot_data = data.frame(OrthogroupID=row.names(annual),
                           Peak_month=annual$phase_month)
    
    if(Tis=="L"){
      ymax = 1000
    } else {
      ymax = 1600
    }
    
    # phenology data = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    phenology_data = read.csv("data/18_datpheno_Lit.csv")
    
    # convert into monthly average
    phenology_data["FloweropeningRatio"] = phenology_data$floweropening/phenology_data$total
    phenology_data["NewLeafRatio"] = phenology_data$new_leaf/phenology_data$total
    
    phenology_data_target = phenology_data[(phenology_data$species==Sp),]
    mean_pheno_table = data.frame(Month=1:12,
                                  FlowerOpening=NA,
                                  NewLeaf=NA)
    for(t in 1:12){
      mean_pheno_table$FlowerOpening[t] = mean(phenology_data_target$FloweropeningRatio[phenology_data_target$month==t],na.rm=T)
      mean_pheno_table$NewLeaf[t] = mean(phenology_data_target$NewLeaf[phenology_data_target$month==t],na.rm=T)
    }
    
    
    color_leaf = colorRampPalette(c("white","#ca3fb9"))(256)[round(mean_pheno_table$FlowerOpening * 255) + 1]
    color_flower = colorRampPalette(c("white","aquamarine3"))(256)[round(mean_pheno_table$NewLeaf * 255) + 1]
    
    ggplot(plot_data,aes(x=factor(Peak_month)))+
      geom_histogram(stat="count",fill="#76A5CA")+
      scale_y_continuous(limits=c(0,ymax))+
      xlab("Peak [Months]")+
      ggtitle(paste(Sp,", ",Tis,",", dim(annual)[1],"genes"))+
      theme_minimal()+
      coord_polar()
    ggsave(paste("figures/rose_plot_annual_",Sp,"_",Tis,".pdf",sep=""),
           width = 4, height = 4)
    
    
    ggplot(plot_data,aes(x=factor(Peak_month)))+
      geom_histogram(stat="count",fill="#76A5CA")+
      scale_y_continuous(limits=c(0,ymax+200))+
      geom_rect(aes(xmin=(1-0.5), xmax=(1+0.5), ymin=1000, ymax=1100), fill=color_leaf[1], color="white")+
      geom_rect(aes(xmin=(2-0.5), xmax=(2+0.5), ymin=1000, ymax=1100), fill=color_leaf[2], color="white")+
      geom_rect(aes(xmin=(3-0.5), xmax=(3+0.5), ymin=1000, ymax=1100), fill=color_leaf[3], color="white")+
      geom_rect(aes(xmin=(4-0.5), xmax=(4+0.5), ymin=1000, ymax=1100), fill=color_leaf[4], color="white")+
      geom_rect(aes(xmin=(5-0.5), xmax=(5+0.5), ymin=1000, ymax=1100), fill=color_leaf[5], color="white")+
      geom_rect(aes(xmin=(6-0.5), xmax=(6+0.5), ymin=1000, ymax=1100), fill=color_leaf[6], color="white")+
      geom_rect(aes(xmin=(7-0.5), xmax=(7+0.5), ymin=1000, ymax=1100), fill=color_leaf[7], color="white")+
      geom_rect(aes(xmin=(8-0.5), xmax=(8+0.5), ymin=1000, ymax=1100), fill=color_leaf[8], color="white")+
      geom_rect(aes(xmin=(9-0.5), xmax=(9+0.5), ymin=1000, ymax=1100), fill=color_leaf[9], color="white")+
      geom_rect(aes(xmin=(10-0.5), xmax=(10+0.5), ymin=1000, ymax=1100), fill=color_leaf[10], color="white")+
      geom_rect(aes(xmin=(11-0.5), xmax=(11+0.5), ymin=1000, ymax=1100), fill=color_leaf[11], color="white")+
      geom_rect(aes(xmin=(12-0.5), xmax=(12+0.5), ymin=1000, ymax=1100), fill=color_leaf[12], color="white")+
      geom_rect(aes(xmin=(1-0.5), xmax=(1+0.5), ymin=1100, ymax=1200), fill=color_flower[1], color="white")+
      geom_rect(aes(xmin=(2-0.5), xmax=(2+0.5), ymin=1100, ymax=1200), fill=color_flower[2], color="white")+
      geom_rect(aes(xmin=(3-0.5), xmax=(3+0.5), ymin=1100, ymax=1200), fill=color_flower[3], color="white")+
      geom_rect(aes(xmin=(4-0.5), xmax=(4+0.5), ymin=1100, ymax=1200), fill=color_flower[4], color="white")+
      geom_rect(aes(xmin=(5-0.5), xmax=(5+0.5), ymin=1100, ymax=1200), fill=color_flower[5], color="white")+
      geom_rect(aes(xmin=(6-0.5), xmax=(6+0.5), ymin=1100, ymax=1200), fill=color_flower[6], color="white")+
      geom_rect(aes(xmin=(7-0.5), xmax=(7+0.5), ymin=1100, ymax=1200), fill=color_flower[7], color="white")+
      geom_rect(aes(xmin=(8-0.5), xmax=(8+0.5), ymin=1100, ymax=1200), fill=color_flower[8], color="white")+
      geom_rect(aes(xmin=(9-0.5), xmax=(9+0.5), ymin=1100, ymax=1200), fill=color_flower[9], color="white")+
      geom_rect(aes(xmin=(10-0.5), xmax=(10+0.5), ymin=1100, ymax=1200), fill=color_flower[10], color="white")+
      geom_rect(aes(xmin=(11-0.5), xmax=(11+0.5), ymin=1100, ymax=1200), fill=color_flower[11], color="white")+
      geom_rect(aes(xmin=(12-0.5), xmax=(12+0.5), ymin=1100, ymax=1200), fill=color_flower[12], color="white")+
      xlab("Peak [Months]")+
      ggtitle(paste(Sp,", ",Tis,",", dim(annual)[1],"genes"))+
      theme_minimal()+
      coord_polar()
    ggsave(paste("figures/rose_plot_annual_with_phenology_",Sp,"_",Tis,".png",sep=""),
           width = 5, height = 5)
    
    
  }
}


# half-annual # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
Sp_ls = c("Le","Lg","Qg","Qa")
Tis_ls = c("L","B")

for(i in 1:length(Sp_ls)){
  for(j in 1:length(Tis_ls)){
    Sp = Sp_ls[i]
    Tis = Tis_ls[j]
    
    rain_res = read.csv(paste("tables/wave_analysis/rain_result_mean_month_delta11_processed_",Sp,"_",Tis,".csv",sep=""),row.names=1)
    
    sixmonth = rain_res[rain_res$period_type_wide=="short",] # 1,2,(3),9,10,11,12 for Lithocarpus, 5,6,7,8,9,10,11,12 for Quercus
    
    # convert 1-6 (if phase>=6 -> phase = phase-6)
    sixmonth["phase_six"] = sixmonth$phase_month
    sixmonth$phase_six[sixmonth$phase_month>6] = sixmonth$phase_six[sixmonth$phase_month>6]-6
    table(sixmonth$phase_six)
    
    plot_data = data.frame(OrthogroupID=row.names(sixmonth),
                           Peak_month=sixmonth$phase_six)
    
    if(Tis=="L"){
      ymax = 70
    } else {
      ymax = 200
    }
    
    
    ggplot(plot_data,aes(x=factor(Peak_month,levels=c(3,4,5,6,1,2))))+
      geom_histogram(stat="count",fill="#CFA3A3",width=0.5)+
      scale_y_continuous(limits=c(0,ymax))+
      xlab("Peak [Months]")+
      ggtitle(paste(Sp,", ",Tis,",", dim(plot_data)[1],"genes"))+
      theme_classic()+
      theme(aspect.ratio=1.25)
    ggsave(paste("figures/phase_sixmonth_",Sp,"_",Tis,".pdf",sep=""),
           width = 3, height = 3)
    
  }
}

# Comparison of peak month ####################################################
species_pair_ls = c("Qg-Qa","Qg-Le","Qg-Lg","Qa-Lg","Qa-Le","Lg-Le")
Tis_ls = c("L","B")

# two demenssional distribution
for(j in 1:length(Tis_ls)){
  Tis = Tis_ls[j]
  # correlation table
  cor_mat = read.csv(paste("data/pairwise_correlation_",Tis,".csv",sep=""),row.names=1)
  
  for(i in 1:length(species_pair_ls)){
    pair = str_split(species_pair_ls[i],"-")[[1]]
    rain1 = read.csv(paste("tables/wave_analysis/rain_result_mean_month_delta11_processed_",pair[1],"_",Tis,".csv",sep=""),row.names=1)
    rain2 = read.csv(paste("tables/wave_analysis/rain_result_mean_month_delta11_processed_",pair[2],"_",Tis,".csv",sep=""),row.names=1)
    
    # extract OGs that is "annual" in either sp1 or sp2
    target_OGs = row.names(rain1)#[(rain1$period_type_wide=="annual")|(rain2$period_type_wide=="annual")]
    rain1_target = rain1[row.names(rain1)%in%target_OGs,]
    rain2_target = rain2[row.names(rain2)%in%target_OGs,]
    
    count_table = matrix(0,nrow=14,ncol=14)
    r_table = matrix(0,nrow=14,ncol=14)
    row.names(count_table) = -1:12
    colnames(count_table) = -1:12
    row.names(r_table) = -1:12
    colnames(r_table) = -1:12
    
    # count genes with a specific phase
    
    rain1_target$phase_id = rain1_target$phase_month
    rain2_target$phase_id = rain2_target$phase_month
    rain1_target$phase_id[rain1_target$period_type_wide%in%c("aperiodic","long")] = -1
    rain1_target$phase_id[rain1_target$period_type_wide%in%c("short")] = 0
    rain2_target$phase_id[rain2_target$period_type_wide%in%c("aperiodic","long")] = -1
    rain2_target$phase_id[rain2_target$period_type_wide%in%c("short")] = 0
    
    
    for(l in -1:12){
      for(m in -1:12){
        OGs_lm = target_OGs[(rain1_target$phase_id==l)&(rain2_target$phase_id==m)]
        count_table[l+2,m+2] = length(OGs_lm)
        r_table[l+2,m+2] = mean(cor_mat[row.names(cor_mat)%in%OGs_lm, colnames(cor_mat)%in%c(paste(pair[1],"_",pair[2],sep=""),paste(pair[2],"_",pair[1],sep=""))],na.rm=T)
      }
    }
    
    if(Tis=="L"){
      hmax = 500
    } else {
      hmax = 600
    }
    
    count_table_annual = count_table[3:14,3:14]
    r_table_annual = r_table[3:14,3:14]
    
    library(reshape2)
    
    count_table_melt = melt(count_table_annual) # Var1 = pair[1], Var2 = pair[2]
    r_table_melt = melt(r_table_annual) # Var1 = pair[1], Var2 = pair[2]
    count_table_melt["correlation"] = r_table_melt$value
    
    N_minus03 = sum(count_table_melt[(is.na(count_table_melt$correlation)==F)&(count_table_melt$correlation< (-0.3)),]$value)
    M_minus03 = length(count_table_melt[(is.na(count_table_melt$correlation)==F)&(count_table_melt$correlation< (-0.3)),]$value)
    
    count_table_melt$Var1[count_table_melt$Var1==1] = 13
    count_table_melt$Var1[count_table_melt$Var1==2] = 14
    count_table_melt$Var2[count_table_melt$Var2==1] = 13
    count_table_melt$Var2[count_table_melt$Var2==2] = 14
    
    library(ggExtra)
    library(patchwork)
    
    count_table_melt_masked = count_table_melt
    count_table_melt_masked$correlation[(count_table_melt$correlation>(-0.3))] = 0
    
    expression_palette = colorRampPalette(c("#2A1A6D","#3E41A8","dodgerblue3","white","khaki","#DFAD4A","#CE422F"))
    
    if(Tis == "L"){
      max_num = 250
    } else {
      max_num = 450
    }
    
    
    ggplot(count_table_melt[count_table_melt$value>0,], aes(x = Var1, y = Var2,fill=correlation)) +
      #geom_tile(data=count_table_melt_masked[count_table_melt$value>0,], aes(x=Var1,y=Var2,fill=correlation))+
      geom_abline(slope = 1, intercept = 0,color="black")+
      geom_point(aes(size = value), shape = 21, color = "black", alpha = 1) + 
      scale_size_continuous(limits = c(0,max_num), range = c(2, 15), name = "Gene Count") +
      scale_x_continuous(breaks = 3:14, limits = c(3, 14),labels=c(3:12,1:2)) + 
      scale_y_continuous(breaks = 3:14, limits = c(3, 14),labels=c(3:12,1:2)) +
      scale_fill_gradientn(colors=expression_palette(256),limits=c(-1,1))　+
      labs(
        x = paste("Phase in ",pair[1],sep=""),
        y = paste("Phase in ",pair[2],sep=""),
        title = "Gene Count Distribution by Phase Differences"
      ) +
      theme_minimal() +
      theme_bw()+
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      )
    ggsave(paste("figures/phase_and_correlation_distribution/",pair[1],"_",Tis,"_vs_",pair[2],"_",Tis,"ver20250214_ver3.pdf",sep=""),height=6,width=7)
    
  }
}


# Molecular Phenology Divergence Index D ##################################
species_pair_ls = c("Qg-Qa","Qg-Le","Qg-Lg","Qa-Lg","Qa-Le","Lg-Le")
Tis_ls = c("L","B")
library(reshape2)

# two demenssional distribution
for(j in 1:length(Tis_ls)){
  Tis = Tis_ls[j]
  
  for(i in 1:length(species_pair_ls)){
    pair = str_split(species_pair_ls[i],"-")[[1]]
    rain1 = read.csv(paste("tables/wave_analysis/rain_result_mean_month_delta11_processed_",pair[1],"_",Tis,".csv",sep=""),row.names=1)
    rain2 = read.csv(paste("tables/wave_analysis/rain_result_mean_month_delta11_processed_",pair[2],"_",Tis,".csv",sep=""),row.names=1)
    
    # extract OGs that is "annual" in either sp1 or sp2
    target_OGs = row.names(rain1)#[(rain1$period_type_wide=="annual")|(rain2$period_type_wide=="annual")]
    rain1_target = rain1[row.names(rain1)%in%target_OGs,]
    rain2_target = rain2[row.names(rain2)%in%target_OGs,]
    
    count_table = matrix(0,nrow=14,ncol=14)
    row.names(count_table) = -1:12
    colnames(count_table) = -1:12
    
    # count genes with a specific phase
    
    rain1_target$phase_id = rain1_target$phase_month
    rain2_target$phase_id = rain2_target$phase_month
    rain1_target$phase_id[rain1_target$period_type_wide%in%c("aperiodic","long")] = -1
    rain1_target$phase_id[rain1_target$period_type_wide%in%c("short")] = 0
    rain2_target$phase_id[rain2_target$period_type_wide%in%c("aperiodic","long")] = -1
    rain2_target$phase_id[rain2_target$period_type_wide%in%c("short")] = 0
    
    
    for(l in -1:12){
      for(m in -1:12){
        OGs_lm = target_OGs[(rain1_target$phase_id==l)&(rain2_target$phase_id==m)]
        count_table[l+2,m+2] = length(OGs_lm)
      }
    }
    
    count_table_annual = count_table[3:14,3:14]
    
    niche_part_ls1 = numeric(12)
    for(t in 1:12){
      if(t == 1){
        m1m2m3 = c(12,t,t+1)
        A = count_table_annual[c(12,t,t+1),t]
        B = count_table_annual[1:12%not.in%m1m2m3,t]
      } else if(t==12){
        m1m2m3 = c(t-1,t,1)
        A = count_table_annual[c(t-1,t,1),t]
        B = count_table_annual[1:12%not.in%m1m2m3,t]
      } else {
        m1m2m3 = c(t-1,t,t+1)
        A = count_table_annual[c(t-1,t,t+1),t]
        B = count_table_annual[1:12%not.in%m1m2m3,t]
      }
      niche_part = sum(B)/(sum(A)+sum(B))
      niche_part_ls1[t] = niche_part
    }
    
    niche_part_ls2 = numeric(12)
    for(t in 1:12){
      if(t == 1){
        m1m2m3 = c(12,t,t+1)
        A = count_table_annual[t,c(12,t,t+1)]
        B = count_table_annual[t,1:12%not.in%m1m2m3]
      } else if(t==12){
        m1m2m3 = c(t-1,t,1)
        A = count_table_annual[t,c(t-1,t,1)]
        B = count_table_annual[t,1:12%not.in%m1m2m3]
      } else {
        m1m2m3 = c(t-1,t,t+1)
        A = count_table_annual[t,c(t-1,t,t+1)]
        B = count_table_annual[t,1:12%not.in%m1m2m3]
      }
      niche_part = sum(B)/(sum(A)+sum(B))
      niche_part_ls2[t] = niche_part
    }
    
    niche_part_df = data.frame(NichePartition1=niche_part_ls1,
                               NichePartition2=niche_part_ls2)
    colnames(niche_part_df) = c(paste(pair[1],"-",pair[2],sep=""),paste(pair[2],"-",pair[1],sep=""))
    
    
    
    
    if(i==1){
      niche_part_df_all = niche_part_df
    } else {
      niche_part_df_all = cbind(niche_part_df_all,niche_part_df)
    }
    
  }
  
  niche_part_df_all["Month"] = 1:12
  write.csv(niche_part_df_all,paste("tables/Fig3/temporal_niche_partionning_",Tis,"_raw.csv",sep=""))
  
  stat_table = data.frame(Month=1:12,
                          Mean=apply(niche_part_df_all[,1:12],1,mean),
                          SD=apply(niche_part_df_all[,1:12],1,sd))
  stat_table["upper"] = stat_table$Mean + stat_table$SD
  stat_table["lower"] = stat_table$Mean - stat_table$SD
  
  niche_part_df_all_melt = melt(niche_part_df_all,id.vars="Month")
  
  niche_part_df_all_melt$Month[niche_part_df_all_melt$Month==1] = 13
  niche_part_df_all_melt$Month[niche_part_df_all_melt$Month==2] = 14
  stat_table$Month[stat_table$Month==1] = 13
  stat_table$Month[stat_table$Month==2] = 14
  write.csv(stat_table,paste("tables/Fig3/temporal_niche_partionning_",Tis,"_stat.csv",sep=""))
  
  
  niche_part_df_all_melt["genus"] = "different"
  niche_part_df_all_melt$genus[niche_part_df_all_melt$variable%in%c("Le-Lg","Lg-Le","Qg-Qa","Qa-Qg")] = "same"
  
  ggplot()+
    geom_point(data=niche_part_df_all_melt,aes(x=Month,y=value),color="dodgerblue3",alpha=0.7)+
    geom_line(data=stat_table,aes(x=Month,y=Mean),color="dodgerblue3",size=0.8)+
    geom_ribbon(data = stat_table, aes(x = Month, ymin = lower, ymax = upper), alpha = 0.3, fill = "deepskyblue3") +
    scale_x_continuous(breaks=3:14,labels=c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb"))+
    scale_y_continuous(limits=c(0,1))+
    #scale_color_manual(values=c("same"="black","different"="gray50"))+
    labs(x="Month",y="#genes with delta phase<=1 \n /#genes with delta phase>1 ",color="Species pair")+
    ggtitle(paste("Tissue:",Tis))+
    theme_bw()+
    theme(aspect.ratio=0.3)
  ggsave(paste("figures/Fig3/Degree_of_temporal_niche_partitioning_rev_",Tis,"_ver3.pdf",sep=""),height=4,width=8)
  
  ggplot()+
    geom_point(data=niche_part_df_all_melt,aes(x=Month,y=value,color=genus,shape=genus),alpha=0.7)+ #,color="dodgerblue3"
    geom_line(data=stat_table,aes(x=Month,y=Mean),color="dodgerblue3",size=0.8)+
    geom_ribbon(data = stat_table, aes(x = Month, ymin = lower, ymax = upper), alpha = 0.3, fill = "deepskyblue3") +
    scale_x_continuous(breaks=3:14,labels=c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb"))+
    scale_y_continuous(limits=c(0,1))+
    scale_color_manual(values=c("same"="black","different"="gray50"))+
    labs(x="Month",y="#genes with delta phase<=1 \n /#genes with delta phase>1 ",color="Species pair")+
    ggtitle(paste("Tissue:",Tis))+
    theme_bw()+
    theme(aspect.ratio=0.3)
  ggsave(paste("figures/Fig3/Degree_of_temporal_niche_partitioning_rev_",Tis,"_ver4.pdf",sep=""),height=4,width=8)
}



