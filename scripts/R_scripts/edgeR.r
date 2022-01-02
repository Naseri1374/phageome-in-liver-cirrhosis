#edgeR
setRepositories()
options(stringsAsFactors = F)
setwd("E:/count/all_count/")
files <- list.files("." , "*.count")
files
cn1 <-lapply(files, read.delim , header = F, comment.char = "*" )
cn <-do.call(cbind , cn1) 
rownames(cn) <- cn[,1]
colnames(cn) <- sub(".count", "",files)

#group<- c(rep("ALC-Russia",27),rep("NAFLD", 9) )
#library(edgeR)
pheno_data = read_excel("E:/count/allcount.xlsx",col_names = T)
pheno_data<-pheno_data[,2]
a<-c()
a<-as.factor(t(pheno_data))
group<-a
y <- DGEList(counts=cn, group=a)
#TMM <- calcNormFactors(y, method="TMM")
#Thr following command compute the optimal CPM threashold
cpm(10 , mean(y$samples$lib.size))
sample_no = 18
keep <- rowSums(cpm(y) > 44.94) > sample_no
y<- y[keep,,keep.lib.size = FALSE]

#Normalization 
y<-calcNormFactors(y)
plotMDS(y) #for coloring the samples plotMDS(y, col=rep(1:2, each = 3))
#######
#group <- factor(c('1','1', '1' , '2', '2', '2'))
design <- model.matrix(~group)
rownames(design) <- colnames(y)

y<- estimateDisp(y, design , robust = T)
fit <- glmQLFit(y , design , robust = T)
qlf <- glmQLFTest(fit , coef=2:6)
topTags(qlf)
summary(decideTests(qlf))
plotMD(qlf)
abline(h=c(-2,2), col = "blue")
finalDegMatrix<- topTags(qlf , n = 479425)
write.csv(finalDegMatrix , "edgeR_other.csv")
