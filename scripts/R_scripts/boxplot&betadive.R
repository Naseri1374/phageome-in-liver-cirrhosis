setwd("E:/dddd/cirr_group/")
library("phyloseq")
library(ggplot2)
library(plyr)
library("ggplot2")
library("dplyr")
library("readxl")
stringsAsFactors = FALSE
otu_mat<- read_excel("phyloseq.xlsx", sheet = "OTU matrix")
tax_mat<- read_excel("phyloseq.xlsx", sheet = "Taxonomy table")
samples_df <- read_excel("phyloseq.xlsx", sheet = "Samples")
row.names(otu_mat) <- otu_mat$OUT
otu_mat <- otu_mat %>% select (-OUT)
row.names(tax_mat) <- tax_mat$OUT
tax_mat <- tax_mat %>% select (-OUT) 
row.names(samples_df) <- samples_df$sample
samples_df <- samples_df %>% select (-sample) 
library("wordspace")
otu_mat <- as.matrix(otu_mat)
otu_mat<-normalize.rows(otu_mat)
otu_mat = round(otu_mat * 1e+04)
otu_mat<-removeBatchEffect(otu_mat, batch=samples_df$grouping)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
row.names(samples)<- a
carbom <- phyloseq(OTU, TAX,samples)
carbom <- transform_sample_counts(carbom, function(x) x / sum(x))
library(ape)
random_tree = rtree(ntaxa(carbom), rooted=TRUE, tip.label=taxa_names(carbom))
carbom <- phyloseq(OTU, TAX, samples,random_tree)
carbom
carbom <- subset_samples(carbom,grouping  =="other")
carbom <- subset_samples(carbom,grouping %in% c( "Chinese ALC" , "Chinese ALC/HBV" , "HBV" , "NAFLD" , "Russian ALC"))
carbom <- subset_taxa(carbom ,Genus %in% c("Actinomyces" , "Alistipes" , "Bacteroides" ,"Megasphaera"
                                           , "Bifidobacterium" , "Blautia","Campylobacter" ,
                                           "Clostridium" , "Erysipelotrichaceae_noname" , "Eubacterium" , "Fusobacterium" , "Haemophilus"
                                           ,"Lachnospiraceae_noname" ,"Lactobacillus" , "Parabacteroides	" ,"Prevotella" , "Roseburia"
                                           ,"Ruminococcus" ,"Siphoviridae_noname" , "Streptococcus" , "Veillonella"))

carbom <- subset_taxa(carbom ,Phylum %in% c("Firmicutes", "Bacteroidetes", "Proteobacteria",
                                            "Fusobacteria", "Actinobacteria", "Verrucomicrobia"))
carbom->phy
glom <- tax_glom(phy, taxrank = 'Phylum')
dat <- psmelt(glom)

# convert Phylum to a character vector from a factor because R
dat$Phylum <- as.character(dat$Phylum)

# group dataframe by Phylum, calculate median rel. abundance
medians <- ddply(dat, ~Phylum, function(x) c(median=median(x$Abundance)))

# find Phyla whose rel. abund. is less than 1%
#remainder <- medians[medians$median <= 0.01,]$Phylum

# change their name to "Remainder"
#dat[dat$Phylum %in% remainder,]$Phylum <- 'Remainder'
#+ ylim(0,0.3)
# boxplot
ff<-ggplot(dat,aes(x=Phylum,y=Abundance, fill=Phylum)) + geom_boxplot() 
ff + theme_classic() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1 ,size = 12, face = "bold"),
          axis.title.x.bottom = element_text( size = 12, face = "bold"),
          axis.title.y.left =  element_text( size = 12, face = "bold"),
          legend.text=element_text(size=12, face = "bold" ),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank()) +  theme(legend.position = "bottom")+scale_fill_manual( values=c("orange","#0066cc","#cc0000" ,"darkgreen", "#6633CC","pink2")) 
##########################################
#simpson & Shannon with *

primary <- 'grouping'
primary_order_list <- c("Chinese ALC" , "Chinese ALC/HBV" , "HBV" , "NAFLD" , "Russian ALC")
samples$grouping <- factor(samples$grouping, levels = primary_order_list)
levels(samples$grouping)

library(RColorBrewer)
ncolors <- length(levels(samples$grouping))
primary_color_list <- brewer.pal(ncolors, 'Dark2')


#plot_richness(carbom, measures = c('Shannon', 'Simpson'))
#ps_sub.rare <- rarefy_even_depth(carbom)
#sample_sums(ps_sub.rare)
#plot_richness(ps_sub.rare, measures = c('Shannon', 'Simpson'))
 
richness.rare <- cbind(estimate_richness(carbom, 
                                         measures = c('Shannon', 'Chao1')),
                       sample_data(carbom)$grouping)
colnames(richness.rare) <- c('Chao1', 'se.chao1','Shannon', 'grouping')
richness.rare$Labels <- rownames(richness.rare)

ad.test.df <- richness.rare[,c('Chao1', 'se.chao1','Shannon', 'grouping')]
colnames(ad.test.df) <- c('Chao1', 'se.chao1','Shannon', 'grouping')
kruskal.test(Shannon ~ grouping, data=ad.test.df)


shannon.plot <- ggplot(ad.test.df, aes(x = grouping, y = Shannon, fill=grouping)) +
  geom_boxplot() + 
  ylab('Shannon Index') +
  xlab(primary) +
  scale_fill_manual(name = primary, values = primary_color_list) +
  labs(fill = primary) +
  theme_cowplot() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank()) 

library(cowplot)
chao1.plot <- ggplot(ad.test.df, aes(x = grouping, y = Chao1, fill=grouping)) +
  geom_boxplot() + 
  ylab('chao1 Index') +
  xlab(primary) +
  scale_fill_manual(name = primary, values = primary_color_list) +
  labs(fill = primary) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank())
diversity.plots <- plot_grid(shannon.plot, chao1.plot, 
                             labels = c('A', 'B'), align = 'h',
                             rel_widths = c(1.5, 2))  
diversity.plots

ad.wilcox.shannon <- pairwise.wilcox.test(ad.test.df$Shannon,
                                          ad.test.df$grouping,
                                          p.adjust.method = 'fdr')
ad.wilcox.simpson <- pairwise.wilcox.test(ad.test.df$Simpson,
                                          ad.test.df$grouping,
                                          p.adjust.method = 'fdr')
ad.wilcox.shannon

library(ggplot2)
library(ggpubr)
library(tidyverse)
shannon.plot.sig <- shannon.plot + 
  stat_compare_means(method = 'wilcox',
                     label = 'p.signif',
                     comparisons = list(c('ALC-China', 'ALC-Russia'),
                                        c('ALC-Russia', 'both_alchol&HBV'),c('ALC-Russia','HBV'),
                                        c('NAFLD', 'HBV'))) 

ad.wilcox.simpson

simpson.plot.sig <- simpson.plot + 
  stat_compare_means(method = 'wilcox',
                     label = 'p.signif',
                     comparisons = list(c('HBV', 'ALC-Russia'),c('ALC-Russia','ALC-China'),
                                        c('NAFLD', 'HBV')))
                            
diversity.plots.sig <- plot_grid(shannon.plot.sig, simpson.plot.sig, 
                                 labels = c('A', 'B'), align = 'h',
                                 rel_widths = c(1.5, 2))  
diversity.plots.sig
#######################################################################
#Fierst Method
library(vegan)
sol <-metaMDS(t(otu_mat))
plot(sol$points)
NMDS = data.frame(MDS1 = sol$points[,1],MDS2=sol$points[,2],group=samples$grouping)
NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$group),mean)
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
  ,group=g))
}

sp<-ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = group)) +
  geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
  theme_classic() +
  xlim(-1, 1.5) 
beta_diversity<-sp + scale_color_manual(values = primary_color_list) +  theme(legend.position = "bottom")
beta_diversity
#####################
#second Method
plot.new()

ord<-ordiellipse(sol, samples$grouping, display = "sites", 
                 kind = "sd", conf = .95, label = T)

df_ell <- data.frame()
for(g in levels(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center)))
                                ,group=g))
}

plot2<-ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = group)) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=1, linetype=2)+
  annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$group)

plot2
#####################
samp <- data.frame(sample_data(carbom))
primary_order_list <- c('HBV', 'NAFLD', 'other', 'both_alchol&HBV', 'ALC-Russia', 'ALC-China')
samp$grouping <- factor(samp$grouping, levels = primary_order_list)
sample_data(carbom) <- samp
suppressPackageStartupMessages(library(cluster))
dist_bray <- distance(carbom, method = "bray")
dist_js <- distance(carbom, method="jsd")
dist_rjs <- sqrt(dist_js)
ord_bray <- ordinate(carbom, method="PCoA", distance=dist_bray)
ord_JS  <- ordinate(carbom, method="PCoA", distance=dist_js)
ord_RJS <- ordinate(carbom, method="PCoA", distance=dist_rjs)
plot_ordination(carbom, ord_RJS, color="grouping") + 
  ggplot2::ggtitle("Bray-Curtis Principal Coordinates Analysis") 
##########################
#PErmANOVA
library(vegan)
primary_order_list <- c('HBV', 'NAFLD', 'other', 'both_alchol&HBV', 'ALC-Russia', 'ALC-China')
samples$grouping <- factor(samples$grouping, levels = primary_order_list)
levels(samples$grouping)
cbn <- combn(x=levels(samples$grouping), m = 2)
s=matrix(nrow = 1, ncol = 3)
for(i in 1:ncol(cbn)){
  physeq_subs <- subset_samples(carbom, grouping %in% cbn[,i])
  metadata_sub <- as(sample_data(physeq_subs), "data.frame")
  permanova_sites <-adonis(distance(physeq_subs, method="bray") ~ grouping,
                           data = metadata_sub)
  a<-(c(cbn[,i],permanova_sites[["aov.tab"]]$`Pr(>F)`[[1]]))
  s<-rbind(s,a)
}
########################
ordu.unwt.uni <- ordinate(carbom, "PCoA", "unifrac", weighted=F)
unwt.unifrac <- plot_ordination(carbom, 
                                ordu.unwt.uni, color="grouping") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac") + geom_point(size = 2)
unwt.unifrac <- unwt.unifrac + theme_classic() + scale_color_brewer("Location", palette = "Set2")
print(unwt.unifrac)

print(unwt.unifrac + stat_ellipse())


##
metadf <- data.frame(sample_data(carbom))

unifrac.dist <- UniFrac(carbom, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

permanova <- adonis(unifrac.dist ~ grouping, data = metadf)

permanova
betadisper(unifrac.dist, metadf$grouping)
View(permutest(ps.disper, pairwise = TRUE))


OTU1 <- read_excel("phyloseqp.xlsx", sheet = "OTU_1")
row.names(OTU1) <- OTU1$grouping
OTU1 <- OTU1 %>% select (-grouping)
N_OTU1<-t(apply(OTU1, 1, function(x)(x-min(x))/(max(x)-min(x))))
write.csv(N_OTU1,"N_OTU1.csv",row.names = T)
######################

#heatmap in with sluster

library(gplots) 
y <- read_excel("phyloseq.xlsx", sheet = "family")
row.names(y) <- y$grouping
y <- y %>% select (-grouping)
y[is.na(y)] <- 0
y<-as.matrix(t(y))
heatmap.2(y)

hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete") 
## Tree cutting
mycl <- cutree(hr, h=max(hr$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
## Plot heatmap 
mycol <- colorpanel(40, "lightpink1", "lightskyblue2", "firebrick2") # or try redgreen(75)

# 1. high resolution of heatmap
tiff("Heatmap.tif", width = 20, height = 10, units = 'in', res = 300)
# 2. Create the plot
par(oma=c(15,20,1,1));heatmap.2(y, Rowv=as.dendrogram(hr), Colv=NA, col=mycol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc,cexCol = 1,keysize=0.5,lwid=c(2,3.2)) 
# 3. Close the file
dev.off()
##########################
otu_mat<- read_excel("phyloseq.xlsx", sheet = "OTU matrix")
tax_mat<- read_excel("phyloseq.xlsx", sheet = "Taxonomy table")
samples_df <- read_excel("phyloseq.xlsx", sheet = "Samples")
row.names(otu_mat) <- otu_mat$OUT
otu_mat <- otu_mat %>% select (-OUT)
row.names(tax_mat) <- tax_mat$OUT
tax_mat <- tax_mat %>% select (-OUT) 
row.names(samples_df) <- samples_df$sample
samples_df <- samples_df %>% select (-sample) 
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
primary_order_list <- c("Chinese ALC" , "Chinese ALC/HBV" , "HBV" , "NAFLD" , "Russian ALC")
samples$grouping <- factor(samples$grouping, levels = primary_order_list)
carbom <- phyloseq(OTU, TAX,samples)
carbom <- subset_samples(carbom,grouping  =="NAFLD")
#relativ abundance
library(devtools)
library(philr)
library(DECIPHER)
library(MicrobeR)
source('C:/Users/Naseri11/Desktop/Naseri/summarize_taxa.R', echo=TRUE)
source('C:/Users/Naseri11/Desktop/Naseri/merge_less_than_top.R', echo=TRUE)
ps <- tax_glom(carbom, "Genus")
ps0 <- merge_less_than_top(ps, top=60)
ps1 <- merge_samples(ps0, "grouping")
#plot_bar(ps2, fill="Phylum")+ scale_fill_manual( values=c("orange","#0066cc","#cc0000" ,"darkgreen", "#CD5C5C","pink2","brown","#800080","#800000","#808000","brown","orange","#0066cc","#cc0000" ,"darkgreen", "#CD5C5C","pink2"))
summarize_taxa(ps1,"Genus")
##################
library(devtools)
library(microbiomeSeq)
ps<- subset_samples(carbom, grouping%in% c("Russian ALC", "HBV"))
ps@sam_data[, 1][["grouping"]]<-factor(ps@sam_data[, 1][["grouping"]])
primary_order_list <- c("Chinese ALC" , "Chinese ALC/HBV" , "HBV" , "NAFLD" , "Russian ALC")
samples$grouping <- factor(samples$grouping, levels = primary_order_list)
physeq <- taxa_level(ps, "Family")
# "log-relative"

deseq_sig <- differential_abundance(physeq, grouping_column  ="grouping" ,output_norm=NULL,pvalue.threshold=0.05,lfc.threshold=0,filename=F)
deseq_sig <- differential_abundance(physeq, "grouping" ,output_norm="log-relative",pvalue.threshold=0.05,lfc.threshold=0,filename=F)
View(deseq_sig[["SignFeaturesTable"]])

p <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p)
p <- plot_MDA(deseq_sig$importance)
print(p)
p <- plot_MA(deseq_sig$SignFeaturesTable)
print(p$maplot)
print(p$lfcplot)

##################333
carbom = subset_samples(carbom,grouping == "HBV" | grouping == "NAFLD"|grouping == "other"|grouping == "ALC-China")
#carbom<- subset_samples(carbom, grouping%in% c('HBV', 'NAFLD', 'other',"ALC-China"))
carbom = filter_taxa(carbom, function(x) sum(x >= 1) == (2), TRUE)
otu_table(carbom)