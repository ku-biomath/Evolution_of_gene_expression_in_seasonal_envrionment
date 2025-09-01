# R script for Principal Component Analysis (PCA)
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

# packages #####################################################################
library("ggplot2")
library(plotly)
library("scatterplot3d")

# PCA ##########################################################################
# data
data_all = read.csv("data/datamatrix_getmm_1to1_common_LeLgQgQa.csv",row.names=1)
gene_table_in_data = data_all[,1:2]
data = data_all[,3:dim(data_all)[2]]
data = data[,(substr(colnames(data),4,4)=="L")|(substr(colnames(data),4,4)=="B")]

# data transformation
log_exp = log2(data+1)
mean_log_exp = mean_individual(log_exp)
x = mean_log_exp

# climate information
climate_info = read.csv("tables/sample_list_with_climate_information.csv",row.names=1)
climate_info_in_data = climate_info[climate_info$Tissue%in%c("L","B"),]
#all(colnames(x)==climate_info_in_data$SampleID)

# PCA
pca_result = prcomp(t(x),scale=T)
summary(pca_result)

# variance explained
variance_explained = pca_result$sdev^2 / sum(pca_result$sdev^2)
par(pin=c(4,4))
cumulative_variance = numeric(length(variance_explained ))
c_v = 0
for(k in 1:length(variance_explained)){
  c_v = c_v + variance_explained[k]
  cumulative_variance[k] = c_v
}

variance_table = data.frame(PCaxis=1:length(variance_explained),
                            Variance_ratio=variance_explained,
                            Cummulative_variance_ratio=cumulative_variance)
write.csv(variance_table,"tables/Fig2/PCA_variance.csv")

plot_variance = data.frame(PC_axis=seq(1,length(variance_explained ),1),
                           Variance_ratio=variance_explained[1:length(variance_explained )],
                           Cumulative_variance_ratio=cumulative_variance[1:length(variance_explained )])

g = ggplot(plot_variance,aes(x=PC_axis,y=Variance_ratio))
g = g + geom_bar(stat="identity")
g = g + geom_line(aes(x=PC_axis,y=Cumulative_variance_ratio))
g = g + geom_point(aes(x=PC_axis,y=Cumulative_variance_ratio))
g = g + theme_classic()
g = g + theme(aspect.ratio=1)
#g = g + scale_x_continuous(breaks=seq(1,length(variance_explained ),1))
g = g + xlab("PC axis")
g = g + ylab("Variance explaned (ratio)")
plot(g)
ggsave("figures/FigS2/variance_explained_PCA.pdf",width=5,height=5)


g = ggplot(plot_variance[plot_variance$PC_axis<20,],aes(x=PC_axis,y=Variance_ratio))
g = g + geom_bar(stat="identity")
g = g + geom_line(aes(x=PC_axis,y=Cumulative_variance_ratio))
g = g + geom_point(aes(x=PC_axis,y=Cumulative_variance_ratio))
g = g + theme_classic()
g = g + theme(aspect.ratio=1)
#g = g + scale_x_continuous(breaks=seq(1,length(variance_explained ),1))
g = g + xlab("PC axis")
g = g + ylab("Variance explaned (ratio)")
plot(g)
ggsave("figures/FigS2/variance_explained_PCA_20.pdf",width=5,height=5)


# plot PCA
plot_data = data.frame(
  Species = substr(colnames(x),1,2),
  Date = substr(colnames(x),6,11),  # sampling date
  Month = as.numeric(substr(colnames(x),8,9)),
  Temperature = climate_info_in_data$Temperature,
  Tissue = substr(colnames(x),4,4),
  PC1 = pca_result$x[, 1],  # pc1
  PC2 = pca_result$x[, 2],  # pc2
  PC3 = pca_result$x[, 3],  # pc3
  PC4 = pca_result$x[, 4],
  PC5 = pca_result$x[, 5]
)

plot_data$Shape = numeric(dim(plot_data)[1])
plot_data$Shape[plot_data$Species=="Le"] = 16
plot_data$Shape[plot_data$Species=="Lg"] = 17
plot_data$Shape[plot_data$Species=="Qg"] = 18
plot_data$Shape[plot_data$Species=="Qa"] = 15

plot_data$Color = numeric(dim(plot_data)[1])
plot_data$Color[plot_data$Tissue=="L"] = "#355792"
plot_data$Color[plot_data$Tissue=="B"] = "#D29E43"

library(grDevices)

plot_data_L = plot_data[(plot_data$Tissue=="L"),]
palette_function_L <- colorRampPalette(colors = c("cadetblue1","aquamarine2","seagreen3","#4B9A9B"))
plot_data_L$Color <- palette_function_L(100)[cut(plot_data_L$Temperature, breaks = 100, labels = FALSE)]

plot_data_B = plot_data[(plot_data$Tissue=="B"),]
palette_function_B <- colorRampPalette(colors = c("khaki1","#D29E43","lightsalmon2","indianred3"))
plot_data_B$Color <- palette_function_B(100)[cut(plot_data_B$Temperature, breaks = 100, labels = FALSE)]

plot_data_all = rbind(plot_data_L,plot_data_B)

s3d = scatterplot3d(plot_data_all[,6:8],
                    pch = plot_data_all$Shape,
                    color = plot_data_all$Color,
                    angle = 140)
text_position = s3d$xyz.convert(plot_data_all[,6:8])
text_position$x = text_position$x + 0.1
text_position$y = text_position$y + 0.1
text(text_position, labels = plot_data_all$Month,
     cex= 0.8,col = "gray40")
legend(s3d$xyz.convert(-190, 0, -80), legend = c("Le","Lg","Qg","Qa"),
       pch =  c(16,17,18,15))

# colorbar
# Leaf
ggplot(plot_data_L,aes(x=as.Date(Date,format="%y%m%d"),y=Temperature,color=Temperature))+
  geom_point()+
  theme_classic()+
  scale_color_gradientn(colors=palette_function_L(100))
ggsave("figures/ExtendedDataFig5/plot_for_temperature_labels_L.pdf")

# Bud
ggplot(plot_data_B,aes(x=as.Date(Date,format="%y%m%d"),y=Temperature,color=Temperature))+
  geom_point()+
  theme_classic()+
  scale_color_gradientn(colors=palette_function_B(100))
ggsave("figures/ExtendedDataFig5/plot_for_temperature_labels_B.pdf")








