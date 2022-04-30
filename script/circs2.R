library(RCircos)
library(RColorBrewer)
library(circlize)
setwd('./')
ref = read.table("chrbed.txt", header = TRUE, sep = "\t")
chrcolor =colorRampPalette(brewer.pal(9, "Set1"))(length(unique(ref$Genome))+1)[1:length(unique(ref$Genome))]#brewer.pal(9,"Set1")
#bed4=RCircos.Link.Data
names(chrcolor) =unique(ref$Genome)
##############读入文件

##要在文件结尾加上一行，否则警告报错：head
pdf('circos.plot.pdf.pdf',w=8,h=6)
library(ComplexHeatmap);library(gridBase) 
f=colorRamp2(breaks=c(7,39,59),
color=c("blue","white","red"))
lgd <- Legend(title = "gene_density", col_fun = f,at=c(7,39,59))
 circle_size = unit(1, "snpc")
lgd_list = packLegend( lgd, max_height = unit(0.8*h, "inch") )
draw(lgd_list, x = circle_size, just = "left")
#incomplete final line found by readTableHeader on 'C:/Users/liufei/Desktop/CIRCOS/chrlen.txt'
##############绘制一个空圈
circos.par("track.height"=0.8,gap.degree=5,start.degree=86,clock.wise=T,cell.padding=c(0,0,0,0))####gap.degree指的是两个track直接的距离
circos.initializeWithIdeogram(ref,plotType = c("axis"))#, "labels"
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0, xlim[2], 1, col =chrcolor[chr] )#
    circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
        facing = "inside", niceFacing = TRUE)
}, track.height = 0.1, bg.border = NA)

#绘制基因密度：
density = read.table("gene.density.xls", header = TRUE, sep = "\t")

f=colorRamp2(breaks=c(7,39,59),
color=c("blue","white","red"))
circos.genomicTrackPlotRegion(density, stack = TRUE,
panel.fun = function(region, value, ...) {
circos.genomicRect(region, value, col = f(value[[1]]), border = NA, ...)
}
,bg.border=NA,track.height = 0.1)



#加GC含量 画散点图：
gc = read.table("gc.xls", header = TRUE, sep = "\t")
circos.genomicTrackPlotRegion(gc,
panel.fun = function(region, value, ...) {
cex = (value[[1]] - min(value[[1]]))/(max(value[[1]]) - min(value[[1]]))#
i = getI(...)
circos.genomicPoints(region, value, cex = cex, pch = 16, 
#col = "red", ...)},
col = "red", ...)},
track.height = 0.13)

circos.yaxis(
at=c(0.4),
             labels.cex=0.5,
             lwd=0,
             tick.length=0,
             labels.col="black",
             col="black")


#深度 绘制柱状图
depth = read.table("coverage.xls", header = T, sep = "\t")
depth[depth$value>=250,]$value=250

# circos.genomicTrackPlotRegion(depth,
# ylim=c(200,250),
# panel.fun = function(region, value, ...){
# circos.genomicLines(region, value, 
# type = "l",
# ytop =250,ybottom = 200,
# col='#FFDAB9',...)},track.height = 0.1)
circos.genomicTrackPlotRegion(depth,panel.fun = function(region, value, ...){
circos.genomicRect(region, value,ytop =value,  ybottom = 110,col = "orange",border=NA)}, bg.border = NA, track.height =0.08)#,track.height = 0.05


circos.yaxis(
at=c(220),
             labels.cex=0.5,
             lwd=0,
             tick.length=0,
             labels.col="black",
             col="black")

pacbio = read.table("pacbio.coverage.xls", header = TRUE, sep = "\t")
pacbio[pacbio$value>=80,]$value=80
circos.genomicTrackPlotRegion(pacbio,
#ylim=c(20,60),
panel.fun = function(region, value, ...){
circos.genomicLines(region, value, 
type = "l",
#ytop =300,ybottom = 200,
col=brewer.pal(9, "Set1")[3],...)},track.height = 0.1)

circos.yaxis(
at=c(40),
             labels.cex=0.5,
             lwd=0,
             tick.length=0,
             labels.col="black",
             col="black")


#dat = read.table("effector.txt", header = TRUE, sep = "\t")
#dat$value=-1
#dat[dat$Class=='Secretome',]$value=1
#dat$genename=''

Secretome = read.table("Secretome.xls", header = TRUE, sep = "\t")
Effectors = read.table("Effectors.xls", header = TRUE, sep = "\t")


circos.genomicLabels(Secretome,
labels.column=6,
cex=0.25,
#col=col_text,
line_lwd=0.5,line_col="blue", side="inside",connection_height=0.05,labels_height=0.001)

circos.genomicLabels(Effectors,
labels.column=6,
cex=0.25,
#col=col_text,
line_lwd=0.5,line_col="red", side="inside",connection_height=0.05,labels_height=0.001)

legend("center",legend=c("Secretome","Effector"),
#fill=NULL,
col=c("blue","red"),
bty="n",
#border=NA,
cex =0.5,
pch=15)

dev.off()


# pdf('yao.pdf',w=8,h=6)
# #incomplete final line found by readTableHeader on 'C:/Users/liufei/Desktop/CIRCOS/chrlen.txt'
# ##############绘制一个空圈
# circos.par("track.height"=0.8,gap.degree=5,start.degree=86,clock.wise=T,cell.padding=c(0,0,0,0))####gap.degree指的是两个track直接的距离
# circos.initializeWithIdeogram(ref,plotType = c("axis"))
# dat = read.table("effector.txt", header = TRUE, sep = "\t")
# dat$value=-1
# dat[dat$Class=='Secretome',]$value=1
# circos.track(ylim = c(-1.1, 1.1), panel.fun = function(x, y) {
    # value = dat$value
    # circos.barplot(value,dat$chromStart,col = ifelse(value > 0, 2, 3), border = ifelse(value > 0, 2, 3),bar_width =2)#, col = ifelse(value > 0, 2, 3)
# })
# circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    # for(pos in dat$chromStart) {
        # value = dat$value
        # circos.boxplot(value, pos)
    # }
# })
# dev.off()


