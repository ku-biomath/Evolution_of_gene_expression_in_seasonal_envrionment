library("ggplot2")
library("NSM3")

# Input ########################################################################
ME_params_table = read.csv("data/molecular_evolution_params_single.csv")
ME_params_table$OGID = substr(ME_params_table$OGID,1,9)
dim(ME_params_table)

ann = read.csv("data/annotation_with_ClusterID.csv")

N_gene_in_data = dim(ann)[1]
ME_params_table_in_data = data.frame(OrthogroupID=ann$OrthogroupID,
                                     LthedID = ann$GeneID.L.edulis,
                                     QglID = ann$GeneID.Q.glauca,
                                     Nsites=NA,
                                     Ssites=NA,
                                     dN=NA,
                                     dS=NA,
                                     omega=NA)
for(i in 1:N_gene_in_data){
  OG_i = ME_params_table_in_data$OrthogroupID[i]
  ME_params_table_i = ME_params_table[ME_params_table$OGID==OG_i,]
  if(dim(ME_params_table_i)[1]>0){
    ME_params_table_in_data$Nsites[i] = ME_params_table_i$Nsites
    ME_params_table_in_data$Ssites[i] = ME_params_table_i$Ssites
    ME_params_table_in_data$dN[i] = ME_params_table_i$dN
    ME_params_table_in_data$dS[i] = ME_params_table_i$dS
    ME_params_table_in_data$omega[i] = ME_params_table_i$omega
  }
}

# output
head(ME_params_table_in_data)
write.csv(ME_params_table_in_data,"tables/Fig5/molecular_evolution_params_in_data_raw.csv")


# Add gene information #########################################################

# periodicity 
phase_table = data.frame(matrix(NA,nrow=dim(ME_params_table_in_data)[1],ncol=4))
label_ls = c("Le_L","Qg_L","Le_B","Qg_B")

for(i in 1:length(label_ls)){
  rain_res = read.csv(paste("tables/wave_analysis/rain_result_mean_month_delta11_processed_",label_ls[i],".csv",sep=""),row.names=1)
  period_type = rain_res$period_type_wide
  phase = rain_res$phase_month
  phase_id = phase
  phase_id[period_type%in%c("aperiodic","long")] = "aperiodic"
  phase_id[period_type=="short"] = "short"
  phase_table[,i] = phase_id
}

colnames(phase_table) = paste("phase_",label_ls,sep="")
row.names(phase_table) = ME_params_table_in_data$OrthogroupID

# output
head(phase_table,10)

# phase difference 
period_comp_table = data.frame(matrix(NA,nrow=dim(ME_params_table_in_data)[1],ncol=0))
row.names(period_comp_table) = ME_params_table_in_data$OrthogroupID

# annual-annual or not
period_comp_table["annual_annual_L"] = FALSE
period_comp_table$annual_annual_L[(phase_table$phase_Le_L%not.in%c("aperiodic","short"))&(phase_table$phase_Qg_L%not.in%c("aperiodic","short"))] = TRUE
period_comp_table["annual_annual_B"] = FALSE
period_comp_table$annual_annual_B[(phase_table$phase_Le_B%not.in%c("aperiodic","short"))&(phase_table$phase_Qg_B%not.in%c("aperiodic","short"))] = TRUE
period_comp_table["annual_annual_Le"] = FALSE
period_comp_table$annual_annual_Le[(phase_table$phase_Le_L%not.in%c("aperiodic","short"))&(phase_table$phase_Le_B%not.in%c("aperiodic","short"))] = TRUE
period_comp_table["annual_annual_Qg"] = FALSE
period_comp_table$annual_annual_Qg[(phase_table$phase_Qg_L%not.in%c("aperiodic","short"))&(phase_table$phase_Qg_B%not.in%c("aperiodic","short"))] = TRUE

# check
table(period_comp_table$annual_annual_L)
table(period_comp_table$annual_annual_B)
table(period_comp_table$annual_annual_Le)
table(period_comp_table$annual_annual_Qg)

# Peak difference
peak_diff_table_L = read.csv(paste("tables/peak_difference_matrix_","L",".csv",sep=""))
peak_diff_table_B = read.csv(paste("tables/peak_difference_matrix_","B",".csv",sep=""))

period_comp_table["DeltaPhase_species_L"] = peak_diff_table_L["Le-Qg"]
period_comp_table["DeltaPhase_species_B"] = peak_diff_table_B["Le-Qg"]

# average
period_comp_table["DeltaPhase_species_MeanByTissue"] = (period_comp_table$DeltaPhase_species_L + period_comp_table$DeltaPhase_species_B)/2

head(period_comp_table)


# Winter or Not 
season_table = data.frame(matrix(NA,nrow=dim(ME_params_table_in_data)[1],ncol=0))
row.names(season_table) = ME_params_table_in_data$OrthogroupID

season_table["Winter_or_Not_L"] = NA
season_table$Winter_or_Not_L[period_comp_table$annual_annual_L==TRUE] = "Others"
season_table$Winter_or_Not_L[(phase_table$phase_Le_L%in%c(12,1,2))&(phase_table$phase_Qg_L%in%c(12,1,2))] = "WW"
table(season_table$Winter_or_Not_L)

season_table["Winter_or_Not_B"] = NA
season_table$Winter_or_Not_B[period_comp_table$annual_annual_B==TRUE] = "Others"
season_table$Winter_or_Not_B[(phase_table$phase_Le_B%in%c(12,1,2))&(phase_table$phase_Qg_B%in%c(12,1,2))] = "WW"
table(season_table$Winter_or_Not_B)

season_table["Winter_or_Not_LB"] = paste(season_table$Winter_or_Not_L,season_table$Winter_or_Not_B,sep="_")
season_table$Winter_or_Not_LB[(period_comp_table$annual_annual_L==FALSE)|(period_comp_table$annual_annual_B==FALSE)] = NA
season_table$Winter_or_Not_LB[season_table$Winter_or_Not_LB=="WW_WW"] = "Winter_in_both"
season_table$Winter_or_Not_LB[season_table$Winter_or_Not_LB=="Others_WW"] = "Winter_in_only_bud"
season_table$Winter_or_Not_LB[season_table$Winter_or_Not_LB=="WW_Others"] = "Winter_in_only_leaf"
season_table$Winter_or_Not_LB[season_table$Winter_or_Not_LB=="Others_Others"] = "Not_winter_in_both"
table(season_table$Winter_or_Not_LB)

# expression level 
# input
data_all = read.csv("data/datamatrix_getmm_1to1_common_LeLgQgQa.csv",row.names=1)
data = data_all[,3:dim(data_all)[2]]
log_exp = log2(data+1)
mean_log_exp = mean_individual(log_exp)
x = mean_log_exp

# median
median_table = data.frame(matrix(NA,nrow=dim(ME_params_table_in_data)[1],ncol=4))
label_ls = c("Le_L","Qg_L","Le_B","Qg_B")
for(i in 1:length(label_ls)){
  x_target = x[,substr(colnames(x),1,4)==label_ls[i]]
  median_ls = apply(x_target,1,median)
  median_table[,i] = median_ls
}
colnames(median_table) = paste("median_",label_ls,sep="")
row.names(median_table) = ME_params_table_in_data$OrthogroupID

median_table["mean_median_Le"] = (median_table$median_Le_L + median_table$median_Le_B)/2
median_table["mean_median_Qg"] = (median_table$median_Qg_L + median_table$median_Qg_B)/2
median_table["mean_median_L"] = (median_table$median_Le_L + median_table$median_Qg_L)/2
median_table["mean_median_B"] = (median_table$median_Le_B + median_table$median_Qg_B)/2
median_table["mean_median_all"] = (median_table$median_Le_L + median_table$median_Le_B + median_table$median_Qg_L + median_table$median_Qg_B)/4

# mean
mean_table = data.frame(matrix(NA,nrow=dim(ME_params_table_in_data)[1],ncol=4))
label_ls = c("Le_L","Qg_L","Le_B","Qg_B")
for(i in 1:length(label_ls)){
  x_target = x[,substr(colnames(x),1,4)==label_ls[i]]
  mean_ls = apply(x_target,1,mean)
  mean_table[,i] = mean_ls
}
colnames(mean_table) = paste("mean_",label_ls,sep="")
row.names(mean_table) = ME_params_table_in_data$OrthogroupID

# mean Le
x_target = x[,substr(colnames(x),1,4)%in%c("Le_L","Le_B")]
mean_ls = apply(x_target,1,mean)
mean_table["mean_Le"] = mean_ls
# mean Qg
x_target = x[,substr(colnames(x),1,4)%in%c("Qg_L","Qg_B")]
mean_ls = apply(x_target,1,mean)
mean_table["mean_Qg"] = mean_ls
# mean L
x_target = x[,substr(colnames(x),1,4)%in%c("Le_L","Qg_L")]
mean_ls = apply(x_target,1,mean)
mean_table["mean_L"] = mean_ls
# mean B
x_target = x[,substr(colnames(x),1,4)%in%c("Le_B","Qg_B")]
mean_ls = apply(x_target,1,mean)
mean_table["mean_B"] = mean_ls
# mean all
x_target = x[,substr(colnames(x),1,4)%in%c("Le_L","Qg_L","Le_B","Qg_B")]
mean_ls = apply(x_target,1,mean)
mean_table["mean_all"] = mean_ls

plot(mean_table$mean_Le_L,mean_table$mean_Le_B)
plot(mean_table$mean_Le_L,mean_table$mean_Qg_L)

# output
head(median_table)
head(mean_table)


# Difference of expression 
diff_mat_L = read.csv(paste("tables/absolute_difference/absolute_difference_matrix_","L",".csv",sep=""),row.names=1)
diff_mat_B = read.csv(paste("tables/absolute_difference/absolute_difference_matrix_","B",".csv",sep=""),row.names=1)

diff_table = data.frame(diff_Le_Qg_L=diff_mat_L$Le_Qg,
                        diff_Le_Qg_B=diff_mat_B$Le_Qg,
                        mean_difference_between_species=(diff_mat_L$Le_Qg+diff_mat_B$Le_Qg)/2
                        )

row.names(diff_table) = row.names(diff_mat_L)

# combination 
param_table = cbind(cbind(cbind(cbind(cbind(cbind(ME_params_table_in_data,phase_table),period_comp_table),median_table),mean_table),season_table),diff_table)
head(param_table)


# remove outlier ###############################################################
# remove NA
param_table_rmNA = param_table[is.na(param_table$omega)==F,]
param_table_NA = param_table[is.na(param_table$omega)==T,]
N_NA = dim(param_table_NA)[1]
param_table_w99 = param_table_rmNA[param_table_rmNA$omega==99,]
param_table_rmw99 = param_table_rmNA[param_table_rmNA$omega!=99,]
param_table_rm_outlier = param_table_rmNA[(param_table_rmNA$dN>0.0003)&(param_table_rmNA$dS>0.001)&(param_table_rmNA$dS<10),]
N_outlier = dim(param_table_rmw99[(param_table_rmw99$dN<0.0003)|(param_table_rmw99$dS<0.001)|(param_table_rmw99$dS>10),])[1]

# 11223 (final) + 483 (outlier+NA) + 43 (w=99) = 11749 (all)

# plot before filtering
ggplot(data=param_table_rmNA,aes(x=log2(dN),y=log2(dS),color=omega))+
  geom_point()+
  geom_vline(xintercept=c(log2(0.0003)),linetype=2)+ # dN
  geom_hline(yintercept=c(log2(0.001),log2(10)),linetype=2)+ # dS
  theme_bw()+
  theme(aspect.ratio=1)
ggsave("figures/scatter_plot_of_dSdN.pdf",height=6,width=6.5)

# remove outlier
param_table_rm_outlier = param_table_rmNA[(param_table_rmNA$dN>0.0003)&(param_table_rmNA$dS>0.001)&(param_table_rmNA$dS<10),]

# plot after filtering
ggplot(data=param_table_rm_outlier,aes(x=log2(dN),y=log2(dS),color=omega))+
  geom_point()+
  theme_bw()+
  theme(aspect.ratio=1)
ggsave("figures/scatter_plot_of_dSdN_after_filtering.pdf",height=6,width=6.5)

write.csv(param_table_rm_outlier,"tables/molecular_evolution_params_in_data_after_filtering.csv",row.names=F)

# regression: mean expression level v.s. dN/dS #################################
param_set = read.csv("tables/molecular_evolution_params_in_data_after_filtering.csv") 

# Variables
Y_vars = c("omega","dN","dS")  # Objective 
X_vars = c("mean_L","mean_B","mean_Le","mean_Qg","mean_all")

results <- data.frame(
  Objective = character(),
  Explanatory = character(),
  Slope = numeric(),
  Slope_sd = numeric(),
  Intercept = numeric(),
  Intercept_sd = numeric(),
  R_squared = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (Y in Y_vars) {
  for (X in X_vars) {
    # regression model
    Y_ls = param_set[[Y]]
    X_ls = param_set[[X]]
    model <- lm(log2(Y_ls)~X_ls)
    model_summary <- summary(model)
    
    # extract parameters
    intercept <- model_summary$coefficients[1]
    intercept_sd <- model_summary$coefficients[3]
    slope <- model_summary$coefficients[2]
    slope_sd <- model_summary$coefficients[4]
    r_squared <- model_summary$r.squared
    p_value <- model_summary$coefficients[8]
    
    results <- rbind(results, data.frame(
      Objective = Y,
      Explanatory = X,
      Slope = slope,
      Slope_sd = slope_sd,
      Intercept = intercept,
      Intercept_sd = intercept_sd,
      R_squared = r_squared,
      P_value = p_value
    ))
    
    # plot
    plot_data = data.frame(X=X_ls,
                           Y=Y_ls,
                           logY=log2(Y_ls))
    ggplot(plot_data, aes(x = X, y = logY)) +
      geom_point(color = "dodgerblue3", size = 1, alpha = 0.3) +
      geom_smooth(method = "lm", color = "gray20", se = TRUE) +
      labs(title = paste("Linear Regression: ","log2(",Y,")"," ~ ",X,sep=""),
           x = X,
           y = paste("log2(",Y,")",sep="")) +
      theme_classic() +
      annotate("text", x = Inf, y = Inf, label = sprintf("Slope: %.3f\nIntercept: %.3f\nR²: %.3f\np: %.3e", 
                                                         slope, intercept, r_squared, p_value),
               hjust = 1.1, vjust = 1.1, size = 5, color = "black", parse = FALSE)+
      theme(aspect.ratio=1)
    ggsave(paste("figures/linear_reguression_expression_level/",X,"_vs_",Y,"_lm.pdf",sep=""),height=5,width=5)

  }
}

write.csv(results, "tables/linear_regression_mean_expression_level.csv",row.names=F)



# regression: peak difference v.s. dN/dS #######################################
param_set = read.csv("tables/molecular_evolution_params_in_data_after_filtering.csv") 


# Variables
Y_vars = c("omega","dN","dS")  # Objective 
X_vars = c("DeltaPhase_species_L","DeltaPhase_species_B","DeltaPhase_species_MeanByTissue")

results <- data.frame(
  Objective = character(),
  Explanatory = character(),
  Slope = numeric(),
  Slope_sd = numeric(),
  Intercept = numeric(),
  Intercept_sd = numeric(),
  R_squared = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (Y in Y_vars) {
  for (X in X_vars) {
    param_set_annual = param_set[(param_set$annual_annual_Le)&(param_set$annual_annual_Qg),]
    
    # regression
    Y_ls = param_set_annual[[Y]]
    X_ls = param_set_annual[[X]]
    model <- lm(log2(Y_ls)~X_ls)
    model_summary <- summary(model)
    
    # extract parameters
    intercept <- model_summary$coefficients[1]
    intercept_sd <- model_summary$coefficients[3]
    slope <- model_summary$coefficients[2]
    slope_sd <- model_summary$coefficients[4]
    r_squared <- model_summary$r.squared
    p_value <- model_summary$coefficients[8]
    
    results <- rbind(results, data.frame(
      Objective = Y,
      Explanatory = X,
      Slope = slope,
      Slope_sd = slope_sd,
      Intercept = intercept,
      Intercept_sd = intercept_sd,
      R_squared = r_squared,
      P_value = p_value
    ))
    
    # plot
    plot_data = data.frame(X=X_ls,
                           Y=Y_ls,
                           logY=log2(Y_ls),
                           OrthogroupID = param_set_annual$OrthogroupID,
                           Highlight = param_set_annual$Highlight)
    
    ggplot(plot_data, aes(x = X, y = logY)) +
      geom_jitter(color="dodgerblue3",size = 1, alpha = 0.3, width=0.1) +
      geom_jitter(data = subset(plot_data, Highlight=="Target"), aes(x=X,y=logY),color="black",size = 1, width=0.1)+
      geom_smooth(method = "lm", color = "gray20", se = TRUE) +
      labs(title = paste("Linear Regression: ","log2(",Y,")"," ~ ",X,sep=""),
           x = X,
           y = paste("log2(",Y,")",sep="")) +
      theme_classic() +
      annotate("text", x = Inf, y = Inf, label = sprintf("Slope: %.3f\nIntercept: %.3f\nR²: %.5f\np: %.3e", 
                                                         slope, intercept, r_squared, p_value),
               hjust = 1.1, vjust = 1.1, size = 5, color = "black", parse = FALSE)+
      theme(aspect.ratio=1)
    ggsave(paste("figures/linear_regression_phase_difference/",X,"_vs_",Y,"_lm.pdf",sep=""),height=5,width=8)
    
  }
}

write.csv(results, "tables/linear_regression_phase_difference.csv",row.names=F)


# Multiple comparison of dN/dS: winter or not in both or either tissue #########

param_set = read.csv("tables/Fig5/molecular_evolution_params_in_data_after_filtering.csv") 

# Variables
Y_vars = c("omega","dN","dS")
X = "Winter_or_Not_LB"

stats <- data.frame(
  EvolutionaryParams = character(),
  CategoryName = character(),
  Category = character(),
  Min = numeric(),
  Mean = numeric(),
  Median = numeric(),
  Max = numeric(),
  SD = numeric(),
  N = numeric(),
  stringsAsFactors = FALSE
)

kw_results = data.frame(
  EvolutionaryParams = character(),
  CategoryName = character(),
  CategoryPair = character(),
  Chi_square = numeric(),
  df = numeric(),
  P_value = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE
)

multcomp_results = data.frame(
  EvolutionaryParams = character(),
  CategoryName = character(),
  CategoryPair = character(),
  W_value = numeric(),
  P_value = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE
)


for(Y in Y_vars){
  # target: genes with conserved peaks (peak difference <= 2)
  param_set_annual = param_set[(param_set$annual_annual_L==TRUE)&(param_set$annual_annual_B==TRUE)&(param_set$DeltaPhase_species_MeanByTissue<=2),]
  table(param_set_annual$Winter_or_Not_LB)
  Y_ls = param_set_annual[[Y]]
  X_ls = param_set_annual[[X]]
  
  # statistics
  for(x in unique(X_ls)){
    stats <- rbind(stats, data.frame(
      EvolutionaryParams = Y,
      CategoryName = X,
      Category = x,
      Min = min(Y_ls[X_ls==x]),
      Mean = mean(Y_ls[X_ls==x]),
      Median = median(Y_ls[X_ls==x]),
      Max = max(Y_ls[X_ls==x]),
      SD = sd(Y_ls[X_ls==x]),
      N = length(Y_ls[X_ls==x])
    ))
  }
  
  # Kruskal Wallis test
  kr = kruskal.test(Y_ls ~ as.factor(X_ls))
  chi_sq = kr$statistic[[1]]
  df = kr$parameter[[1]]
  p_value = kr$p.value[[1]]
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
  
  kw_results = rbind(kw_results, data.frame(
    EvolutionaryParams = Y,
    CategoryName = X,
    Chi_square = chi_sq,
    df = df,
    P_value = p_value,
    Significance = signif,
    stringsAsFactors = FALSE
  ))
  
  # Steel Dwass test
  res = pSDCFlig(Y_ls, as.factor(X_ls), method="Asymptotic")
  p_val = res$p.val
  Dwass_Steel_W_val = res$obs.stat
  Dwass_Steel_Method = res$method
  signif = rep("",length(p_val))
  signif[p_val<0.1] = "."
  signif[p_val<0.05] = "*"
  signif[p_val<0.01] = "**"
  signif[p_val<0.001] = "***" 
  
  multcomp_results = rbind(multcomp_results,data.frame(
    EvolutionaryParams=Y,
    CategoryName=X,
    CategoryPair=res$labels,
    W_value=Dwass_Steel_W_val,
    P_value=p_val,
    Significance=signif))
  
  # plot
  plot_data = data.frame(X=X_ls,
                         Y=Y_ls,
                         logY=log2(Y_ls))
  plot_data$X = factor(plot_data$X,levels=c("Not_winter_in_both","Winter_in_only_bud","Winter_in_only_leaf","Winter_in_both"))
  color_set = c("Not_winter_in_both"="gray60","Winter_in_only_bud"="lightblue3","Winter_in_only_leaf"="darkseagreen","Winter_in_both"="dodgerblue3")
  ggplot(plot_data, aes(x = X, y = logY,fill=X, color=X)) +
    geom_jitter( size = 1, alpha = 0.5, width=0.1) +
    geom_boxplot(alpha=0.3,color="black",outliers = FALSE,width=0.4)+
    scale_color_manual(values=color_set)+
    scale_fill_manual(values=color_set)+
    labs(title = paste("log2(",Y,")"," ~ ",X,sep=""),
         x = X,
         y = paste("log2(",Y,")",sep="")) +
    theme_classic() +
    annotate(
      "text",
      x = 4, y =max(plot_data$logY)+0.5,
      label = paste0(
        "p-value of KW test: ", format.pval(p_value, digits = 3)
      ),
      hjust = 1, color = "black", size = 4
    ) +
    theme(aspect.ratio=1,
          axis.text.x  = element_text(angle = 45),
          legend.position = "none")
  ggsave(paste("figures/winter_or_not/",X,"_vs_",Y,"_multicomp.pdf",sep=""),height=5,width=5)
  
}

write.csv(stats, "tables/winter_or_not_both_tissue_stats.csv",row.names=F)
write.csv(kw_results, "tables/winter_or_not_both_tissue_Kruskal_Wallis.csv",row.names=F)
write.csv(multcomp_results, "tables/winter_or_not_both_tissue_Steel_Dwass.csv",row.names=F)

