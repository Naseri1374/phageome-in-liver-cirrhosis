library(readxl)
setwd("E:\\dddd\\cirr_group\\phage")
data<-read_excel("phyloseqp.xlsx" , col_names = T , sheet = "pca" )
row.names(data) <- data$vir
data$vir= NULL
ncol(data)
data[,ncol(data)]
label<-data[,ncol(data)]
data<-as.data.frame(sapply(data[,-ncol(data)], as.numeric))
label<-as.factor(label)
data<-cbind(data,label)
dim(data)
model = prcomp(data[,1:ncol(data)-1], scale=TRUE)
library("factoextra")
res.pca<- model
fviz_pca_ind(res.pca)
fviz_pca_ind(res.pca, label="none", habillage=data$label)
p <- fviz_pca_ind(res.pca, label="none", habillage=data$label,
                  addEllipses=TRUE, ellipse.level=0.95)
print(p)
p + scale_color_brewer(palette="Dark2") +
  theme_minimal()
p + scale_color_brewer(palette="Paired") +
  theme_minimal()
################################
model$x[,1]-> pc1
model$x[,2]-> pc2
model$x[,3]-> pc3
colors <- c("orange","#0066cc","#cc0000" ,"darkgreen", "#6633CC","pink2")
colors <- colors[as.numeric(label)]
#2D plot
plot(pc2,pc3, col=colors,xlim = c(0,2))


primary <- 'label'
primary_order_list <- c('HBV', 'NAFLD', 'both_alchol&HBV', 'Alcohol-Russia', 'Alcohol-China')
data$label <- factor(data$label, levels = primary_order_list)
levels(data$label)
