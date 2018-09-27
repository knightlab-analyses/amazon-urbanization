# #############################################################
#                                                             #
#    Create heatmap of in silico predicted chemical classes   #
#                                                             #
###############################################################

# load libraries
library(dendextend)
library(vegan)
library(gplots)
library(RColorBrewer)
library(plyr)
library(plotrix)

################################
# Superclass Heatmap           #
################################

metadata <- read.csv("full_norm_Apr25metadata_noPeruvianSwabs_housing_metadata.csv",sep=",")

ft <- read.csv("featuretable_superclass.tsv",sep="\t",row.names = 1)
ft <- t(ft)

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 15) 

snames <- as.character(metadata$X.SampleID[match(rownames(ft),as.character(metadata$Sample.name))])
snames[is.na(snames)] <- "blank"

## change according to metadata category
Rowside_species <- as.character(metadata$village_socio[match(rownames(ft),as.character(metadata$Sample.name))])
Rowside_species[is.na(Rowside_species)] <- "blank"
rownames(ft) <- snames

# order featuretable according to village_socio category
Rowside_species[which(Rowside_species == "blank")] <- "1_blank"
Rowside_species[which(Rowside_species == "Checherta")] <- "2_Checherta"
Rowside_species[which(Rowside_species == "Puerto Almendras")] <- "3_Puerto Almendras"
Rowside_species[which(Rowside_species == "Iquitos")] <- "4_Iquitos"
Rowside_species[which(Rowside_species == "Manaus low")] <- "5_Manaus low"
Rowside_species[which(Rowside_species == "Manaus middle")] <- "6_Manaus middle"

ft <- ft[order(Rowside_species),]
Rowside_species <- Rowside_species[order(Rowside_species)]

colpal <- c("#1b9e77","#FF0000", "#008000", "#0000FF","#FFFF00","#f27304")
Rowside_cols <- colpal[as.numeric(as.factor(Rowside_species))]

## import ClassyFire scores and map them on the heatmap (calculate average per chemical class)
cl <- read.table("ClassyFire_InputforCytoscape_Amazon.csv",sep="\t",header=T,comment.char = "",stringsAsFactors = F,quote="")
cl$CF_superclass_score <- as.numeric(cl$CF_superclass_score)

scores <- ddply(cl, .(CF_superclass), summarize,  score=mean(CF_superclass_score))
Colside_cols <-gray(1-scores$score[match(colnames(ft),scores$CF_superclass)])

Nnodes <- table(cl$CF_superclass)
NnodesNam <- names(Nnodes)
Nnodes <- as.vector(Nnodes)
names(Nnodes) <- NnodesNam

colnames(ft) <- paste(Nnodes[match(colnames(ft),names(Nnodes))],colnames(ft), sep = " / ")

pdf(file="Superclass_HeatMapDendrogram.pdf", width=9, height=10)
par(mar=c(0, 0, 0, 0))
heatmap.2(ft,Rowv=FALSE, Colv=TRUE,cexCol = 1,cexRow = 0.3,scale="col",col = my_palette, RowSideColors=Rowside_cols, ColSideColors=Colside_cols, keysize = 0.8, key.ylab = NA, key.title = NA,density.info="none",tracecol=NA,margin=c(5,15)) #margin=c(5,15)
legend(y=0.9, x=0.9, xpd=TRUE,      
       legend = sort(unique(Rowside_species)),
       col = colpal, 
       lty= 1,             
       lwd = 5,           
       cex=.7)
color.legend(0.9,0.93,1,0.96,c(0,0.5,1),gray.colors(100, start = 1, end = 0),cex = 0.7,align="rb")
text(x=0.95, y = 0.98,labels = "ClassyFire Score",cex=0.8, font =2)
dev.off()
