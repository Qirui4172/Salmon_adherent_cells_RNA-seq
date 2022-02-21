#!/usr/bin/env Rscript

#--------------------------------------------------------------------------------------------------
Usage<-function(){
	cat("\n\tUsage: Rscript salmonDI_mRNA_part2.r <readCounts.mx> <sampleInfo> <readCounts/gene> <nonZeroLib/gene> <controlGroup> <adjPvalue> <foldChange> <ssaGeneName>","\n\n",

		"\tParameters (all required)","\n",
		"\t<readCounts.mx>		Read counts matrix file (\"salmonDI_mRNA_readcount.mx\", generated from \"salmonDI_mRNA_part1.sh\")","\n",
		"\t<sampleInfo>		Sample information file (\"sample.info\")","\n",
		"\t<readCounts/gene>	Minimal total read counts in all 12 libraries per gene, genes with read counts less than this number will be filtered out (suggest: 100)","\n",
		"\t<nonZeroLib/gene>	Minimal non-zero libraries per gene, genes with non-zero libraries less than this number will be filtered out (suggest: 5)","\n",
		"\t<controlGroup>		Control group, \"DI\" or \"HK\" (DI: distal intestine, HK: head kidney)","\n",
		"\t<adjPvalue>			BH adjusted p value for running DESeq2 (suggest: 0.05)","\n",
		"\t<foldChange>		Minimal fold change threshold for determining DEGs (suggest: 2)","\n",
		"\t<ssaGeneName>		\"ssa_geneName.list\"","\n\n",

		"\tExample","\n",
		"\tRscript salmonDI_mRNA_part2.r salmonDI_mRNA_readcount.mx sample.info 100 5 DI 0.05 2 ssa_geneName.list","\n\n",

		"\tFunction","\n",
		"\tRun DESeq2, to find differentially expressed genes and generate plots.","\n\n",

		"\tContact: Qirui Zhang (qirui.zhang@med.lu.se)","\n",
		"\tUpdated: 10-06-2020","\n\n"
	)
}

args<-commandArgs(TRUE)
if (length(args)!=8){Usage();quit();}

#--------------------------------------------------------------------------------------------------
# Load libraries
time<-format(Sys.time(), format='%Y-%m-%d %H:%M:%S')
cat(time, "Start analysis", "\n")
cat("Loading libraries ...", "\n")
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
pdf("salmonDI_mRNA_DEGplots.pdf")

#--------------------------------------------------------------------------------------------------
cat("\n", "================================================================================", "\n")
# read data
time<-format(Sys.time(), format='%Y-%m-%d %H:%M:%S')
cat(time, "Reading data ...", "\n")

readCounts<-read.table(args[1], header=T)
sampleInfo<-read.table(args[2], header=T)
sampleInfo$Replicate<-as.factor(sampleInfo$Replicate)
colnames(readCounts)<-sampleInfo$Sample
readCounts<-readCounts[which(rowSums(readCounts) >= as.numeric(args[3])),]
readCounts<-readCounts[which(rowSums(readCounts != 0) >= as.numeric(args[4])),]

#--------------------------------------------------------------------------------------------------
cat("\n", "================================================================================", "\n")
# run DESeq2
time<-format(Sys.time(), format='%Y-%m-%d %H:%M:%S')
cat(time, "Generating metadata, runing DESeq2, and normalizing read counts ...", "\n")

dds<-DESeqDataSetFromMatrix(countData=readCounts,colData=sampleInfo,design=~Tissue)
dds$Tissue<-relevel(dds$Tissue, args[5])
dds<-DESeq(dds)

normalized.counts<-counts(dds, normalized=TRUE)
write.table(as.data.frame(normalized.counts), "salmonDI_mRNA_baseMean.tsv", row.names=T, col.names=T, quote=F, sep="\t")

#--------------------------------------------------------------------------------------------------
cat("\n", "================================================================================", "\n")
# extract comparison results
time<-format(Sys.time(), format='%Y-%m-%d %H:%M:%S')
cat(time, "Extracting comparison results ...", "\n")

adjP<-as.numeric(args[6])
FC<-as.numeric(args[7])
res<-results(dds,alpha=adjP)
summary(res)
deg.num<-sum(res$padj < adjP & abs(res$log2FoldChange) >= log2(FC), na.rm=TRUE)
cat("\t","Total DEG num: ",deg.num,"\n")
up.num<-sum(res$padj < adjP & res$log2FoldChange>=log2(FC), na.rm=TRUE)
cat("\t","Up-regulated DEG num: ",up.num,"\n")
down.num<-sum(res$padj < adjP & -(res$log2FoldChange)>=log2(FC), na.rm=TRUE)
cat("\t","Down-regulated DEG num: ",down.num,"\n\n")
geneName=read.table(args[8], header=F)
geneName$V1=paste(geneName$V1, geneName$V2, sep=":")
rownames(geneName)=geneName$V1
res$gene<-geneName[rownames(res), 3]
deg.all<-res[which(res$padj < adjP & abs(res$log2FoldChange)>=log2(FC)),]
write.table(as.data.frame(deg.all), "salmonDI_mRNA_DEGlist.tsv", row.names=T, col.names=T, quote=F, sep="\t")
#--------------------------------------------------------------------------------------------------
cat("\n", "================================================================================", "\n")
# make sample plots
time<-format(Sys.time(), format='%Y-%m-%d %H:%M:%S')
cat(time, "Generate sample plots:", "\n")

# transform data
cat("Transforming data with DESeq2::rlog ...", "\n")
rld<-rlog(dds, blind=FALSE)

# plot dispersion estimate
cat("Dispersion estimate plot ...", "\n")
plotDispEsts(dds)

# PCA plots
cat("PCA plot ...", "\n")
PCA_Plot<-function(rld, DI.color, HK.color){
	data<-plotPCA(rld, intgroup="Tissue", returnData=TRUE)
	percentVar<-round(100*attr(data, "percentVar"))
	DI.color<-DI.color
	HK.color<-HK.color
	ggplot(data, aes(PC1, PC2, color=Tissue))+geom_point(size=5)+scale_colour_manual(values=c("DI"={DI.color},"HK"={HK.color}))+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(size=1),text=element_text(size=18),axis.text=element_text(size=15),legend.text=element_text(size=15))+xlab(paste0("PC1: ",percentVar[1],"% variance"))+ylab(paste0("PC2: ",percentVar[2],"% variance"))+coord_fixed(ratio=3)
}
PCA_Plot(rld, "indianred1", "turquoise3")
PCA_Plot(rld, "salmon", "firebrick")

# correlation heatmap
cat("Sample correlation heatmap ...", "\n")
sampleDists<-dist(t(assay(rld)))
sampleDistMatrix<-as.matrix(sampleDists)
rownames(sampleDistMatrix)<-rld$Sample
colnames(sampleDistMatrix)<-NULL
colors<-colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors)

#--------------------------------------------------------------------------------------------------
cat("\n", "================================================================================", "\n")
# DEG plots
time<-format(Sys.time(), format='%Y-%m-%d %H:%M:%S')
cat(time, "Generate DEG plots:", "\n")

# MA plot
cat("MA plot ...", "\n")
plotMA(res, main="MA plot", ylim=c(-10,10))

# volcano plot version1
cat("volcano plot version1 ...", "\n")

VolcanoPlotV1<-function(res, down.color, up.color){
	# divide groups
	volcano<-as.data.frame(res)
	volcano$significant<-as.factor(ifelse(!is.na(volcano$padj) & volcano$padj < adjP & abs(volcano$log2FoldChange)>=log2(FC), ifelse(volcano$log2FoldChange>=log2(FC), "Up", "Down"), "No"))
	volcano$padj<-ifelse(is.na(volcano$padj),1,volcano$padj)
	volcano$group<-rep(0)
	volcano[which(volcano$padj >= adjP | (volcano$padj < adjP & abs(volcano$log2FoldChange) < log2(FC))),9]=rep(1)
	volcano[which(volcano$padj < adjP & volcano$padj >= 1e-100 & volcano$log2FoldChange <= -1 & volcano$log2FoldChange >= -12),9]=rep(2)
	volcano[which(volcano$padj < adjP & volcano$padj >= 1e-100 & volcano$log2FoldChange <= 12 & volcano$log2FoldChange >= 1),9]=rep(3)
	volcano[which(volcano$padj < adjP & volcano$padj >= 1e-100 & volcano$log2FoldChange < -12),c("log2FoldChange", "group")]=list(log2FoldChange = -12, group=4)
	volcano[which(volcano$padj < 1e-100 & volcano$log2FoldChange <= -1 & volcano$log2FoldChange >= -12),c("padj","group")]=list(padj = 1e-100, group=4)
	volcano[which(volcano$padj < 1e-100 & volcano$log2FoldChange < -12),c("log2FoldChange", "padj", "group")]=list(log2FoldChange = -12, padj = 1e-100, group=4)
	volcano[which(volcano$padj < adjP & volcano$padj >= 1e-100 & volcano$log2FoldChange > 12),c("log2FoldChange", "group")]=list(log2FoldChange = 12, group=5)
	volcano[which(volcano$padj < 1e-100 & volcano$log2FoldChange <= 12 & volcano$log2FoldChange >= 1),c("padj", "group")]=list(padj = 1e-100, group=5)
	volcano[which(volcano$padj < 1e-100 & volcano$log2FoldChange > 12),c("padj", "log2FoldChange", "group")]=list(log2FoldChange = 12, padj = 1e-100, group=5)

	# plot
	down.color<-down.color
	up.color<-up.color
	p<-ggplot(volcano,aes(log2FoldChange,-log10(padj)))+geom_point(data=volcano[which(volcano$group==1),],color="gray",alpha=0.75)+geom_point(data=volcano[which(volcano$group==2),],color={down.color},alpha=0.75)+geom_point(data=volcano[which(volcano$group==3),],color={up.color},alpha=0.75)+geom_point(data=volcano[which(volcano$group==4),],shape=2,color={down.color},alpha=0.75)+geom_point(data=volcano[which(volcano$group==5),],shape=2,color={up.color},alpha=0.75)
	p+labs(title="Volcano plot",x="log2FoldChange",y="-log10(padj)")+geom_hline(yintercept=-log10(adjP),linetype=2,color="gray80")+geom_vline(xintercept=c(-log2(FC),log2(FC)),linetype=2,color="gray80")+xlim(-12,12)+ylim(0,100)+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(size=1))
}

VolcanoPlotV1(res, "royalblue3", brewer.pal(11,"RdYlBu")[2])
VolcanoPlotV1(res, "firebrick", "salmon")

#--------------------------------------------------------------------------------------------------
# volcano plot version2
cat("volcano plot version2 ...", "\n")
VolcanoPlotV2<-function(res){
	volcano<-as.data.frame(res)
	volcano$significant<-as.factor(ifelse(!is.na(volcano$padj) & volcano$padj < adjP & abs(volcano$log2FoldChange)>=log2(FC), "yes", "no"))
	volcano$padj<-ifelse(is.na(volcano$padj),1,volcano$padj)

	p<-ggplot(volcano, aes(x=log2FoldChange, y=log2(baseMean)))+geom_point(data=subset(volcano, significant=="no"), colour="gray50")+geom_point(data=subset(volcano, significant=="yes"), aes(colour=padj), alpha=0.8)+scale_colour_gradientn("Padj", colours=c(rev(brewer.pal(9, "OrRd"))[c(2,5,6,7,8,9)]), limits=c(0, 0.05), breaks=c(0, 0.01, 0.02, 0.03, 0.04, 0.05), labels=c("0", "0.01", "0.02", "0.03", "0.04", "0.05"))
	p+labs(title="Volcano plot version2",x="log2FoldChange",y="log2baseMean")+geom_vline(xintercept=c(-1,1),linetype=2,color="black")+theme_bw()+theme(axis.line=element_line(colour="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank())+coord_flip()
}
VolcanoPlotV2(res)

# log2FoldChange density curve
cat("log2FoldChange density curve ...", "\n")
density<-as.data.frame(res)
ggplot(density, aes(x=log2FoldChange))+geom_density(fill="gray50")+labs(title="log2FoldChange density",x="log2FoldChange",y="Density")+geom_vline(xintercept=c(-1,1),linetype=2,color="black")+theme_bw()+theme(axis.line=element_line(colour="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank(),aspect.ratio=5)+coord_flip()

# DEG heatmap
cat("DEG heatmap ...", "\n")
HeatmapPlot<-function(rld, DI_anno.color, HK_anno.color, heatmap.color){
	rld.deg<-rld[rownames(deg.all),]
	rownames(rld.deg)<-deg.all$gene
	anno.label<-data.frame(Tissue=c(rep("DI", 6),rep("HK", 6)))
	rownames(anno.label)<-rownames(colData(rld))

	anno.color<-list(Tissue=c(DI=DI_anno.color, HK=HK_anno.color))
	pheatmap(assay(rld.deg),scale="row",main="Heatmap of DEGs",color=heatmap.color,cluster_cols=F,show_rownames=F,fontsize_row=0.5,annotation_col=anno.label,annotation_colors=anno.color,annotation_names_col=F,cellwidth=15,border_color=NA)
}

DI_anno.color<-brewer.pal(8,"Dark2")[1]
HK_anno.color<-brewer.pal(8,"Dark2")[2]
heatmap.color<-colorRampPalette(c(rev(brewer.pal(9,"Blues")[c(2:9)]), brewer.pal(9, "OrRd")[c(2:9)]))(100)
HeatmapPlot(rld, DI_anno.color, HK_anno.color, heatmap.color)

DI_anno.color<-"salmon"
HK_anno.color<-"firebrick"
heatmap.color<-colorRampPalette(brewer.pal(9,"YlGnBu"))(100)
HeatmapPlot(rld, DI_anno.color, HK_anno.color, heatmap.color)

#--------------------------------------------------------------------------------------------------
cat("\n", "================================================================================", "\n")
time<-format(Sys.time(), format='%Y-%m-%d %H:%M:%S')
cat(time, "Done with the whole analysis!", "\n")
cat("朴英鎮 oppa, your \"salmonDI_mRNA_part2.r\" has successfully finished!\nHope you get good results, and enjoy lovely Bodø winter ^^","\n\n")

