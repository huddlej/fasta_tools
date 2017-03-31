#plotIdeogram.R



#Last updated: 5/14/13 by BLD

#Purpose: Plot significant genomic regions on an ideogram

#Usage:  Rscript plotIdeogram.R <output_pdf> <ideogram file> <significant regions>

#load libraries
library(RColorBrewer)

##########################################    NOTES    ###################################################
#
#1. Significant regions file must be tab delimited with the following header line: "chr	start	stop"
#2. Ideogram file can be downloaded from UCSC Genome Browser Table and stored locally. Be sure to
#   delete the hash tag in the header line so that it is recognized by R. Also, comment out or delete all the
#   lines corresponding to unplaced contigs, alternative haplotypes, and chrY.
#3. Currently, this does NOT plot the Y chromosome, just chr1-22 and chrX.
#
################################### parse command line args ###############################################

args <- commandArgs(trailingOnly = TRUE)
outPdf=args[1]
ideogramFile=args[2]
sigRegions=args[3]


#******************************* Plotting variables to edit manually *************************************#

   #add figure title here. Otherwise, leave as empty double quotes.
   figName=""

   #colors
   colors <- list("insertions" = "#377eb8",
                  "gaps" = "#000000",
                  "inversions" = "#4daf4a")

   nColors=5 # 3 <= nColors <= 12
   heatmap_colors = brewer.pal(n=nColors, name="YlOrRd")
   #heatmap_colors = sample(heatmap_colors, size=length(heatmap_colors), replace=F)
   print(heatmap_colors)

   #specify how to plot significant regions
   lines = 1  #If regions to display on ideograms are large, set lines=1. Otherwise lines=0 will plot an
              #asterisk next to each significant region for improved visibility.

   fixed = 1  #If there are a lot of regions to show, set fixed=0. This will randomly draw distances from
              #a uniform distribution and plot the lines/asterisks at that distance from the ideogram. This
              #prevents adjacent or overlapping lines from completely obscuring one another.

   base_line_width = 3 # Width of line to display for each region

#********************************************************************************************************#


######################################### plot the ideograms #############################################

pdf (file = outPdf, height = 7, width = 15)

stain = read.table(file = ideogramFile, header = T, sep="\t")

chr_lengths = c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560);

plot(x=c(0.75,length(chr_lengths)), y = c(0,max(chr_lengths)), type = "n", xlab = "", ylab = "", main = figName, axes = F)
for (c in seq(1,23,1)) {
   rect(xleft = (c-0.15), ybottom = max(chr_lengths)-chr_lengths[c], xright = c+0.15, ytop = max(chr_lengths), border = "black", lwd = 1)
}

mtext(text = c(seq(1,22,1),"X"), side = 3, line = 0, at = seq(1,23,1), cex = 1.25, padj = 0.5)

#add the bands
for (i in seq(1,length(stain$chrom),by=1)) {
   chr = sub("chr","",stain$chrom[i])
   if (chr == "X") {chr = 23}
   if (chr != "X") {chr = as.numeric(chr)}

   if (stain$gieStain[i] == "gneg") { rect(xleft = (chr-0.15), ybottom = max(chr_lengths)-stain$chromEnd[i], xright = chr+0.15, ytop = max(chr_lengths)-stain$chromStart[i], border = "black", lwd = 1, col = "gray90")}
   if (stain$gieStain[i] == "gpos25") {rect(xleft = (chr-0.15), ybottom = max(chr_lengths)-stain$chromEnd[i], xright = chr+0.15, ytop = max(chr_lengths)-stain$chromStart[i], border = "black", lwd = 1, col = "gray75")}
   if (stain$gieStain[i] == "gpos50") {rect(xleft = (chr-0.15), ybottom = max(chr_lengths)-stain$chromEnd[i], xright = chr+0.15, ytop = max(chr_lengths)-stain$chromStart[i], border = "black", lwd = 1, col = "gray60")}
   if (stain$gieStain[i] == "gpos75") {rect(xleft = (chr-0.15), ybottom = max(chr_lengths)-stain$chromEnd[i], xright = chr+0.15, ytop = max(chr_lengths)-stain$chromStart[i], border = "black", lwd = 1, col = "gray45")}
   if (stain$gieStain[i] == "gpos100") {rect(xleft = (chr-0.15), ybottom = max(chr_lengths)-stain$chromEnd[i], xright = chr+0.15, ytop = max(chr_lengths)-stain$chromStart[i], border = "black", lwd = 1, col = "black")}
   if (stain$gieStain[i] == "acen") {rect(xleft = (chr-0.15), ybottom = max(chr_lengths)-stain$chromEnd[i], xright = chr+0.15, ytop = max(chr_lengths)-stain$chromStart[i], border = "black", lwd = 1, col = "red")}
   if (stain$gieStain[i] == "gvar") {rect(xleft = (chr-0.15), ybottom = max(chr_lengths)-stain$chromEnd[i], xright = chr+0.15, ytop = max(chr_lengths)-stain$chromStart[i], border = "black", lwd = 1, col = "black")}
   if (stain$gieStain[i] == "stalk") {rect(xleft = (chr-0.15), ybottom = max(chr_lengths)-stain$chromEnd[i], xright = chr+0.15, ytop = max(chr_lengths)-stain$chromStart[i], border = "black", lwd = 1, col = "gray60")}
}

#add physical axis scale bars
max = max(chr_lengths)
par(las = 1)
axis(side = 2, at = c(max-250000000,max-200000000,max-150000000,max-100000000,max-50000000,max), labels = c(250,200,150,100,50,0), lwd = 2)
par(las = 0)
mtext(text = "Mb", side = 2, line = 3, cex = 1.3, padj = 0.5)

############################ plot the significant regions from the input file ###############################

sigRegions = read.table(file=sigRegions, header = TRUE, sep = "\t")

base_offset=0.15
offset_increment=0.1

for (row in seq(1,length(sigRegions$chr))) {
   chr = sub("chr","",sigRegions$chr[row])
   if (chr == "X") {chr = 23}
   if (chr != "X") {chr = as.numeric(chr)}

   region_type <- as.character(sigRegions$type[row])
   if (region_type == "STRs") {
       # STRs are always on the inside (or the right?)
       color <- heatmap_colors[[sigRegions$bin[row]]]
       offset <- 0
       line_end <- 1
       line_width=11
   }
   else {
       # set offset based on region type
       color <- colors[[region_type]]
       offset <- base_offset + as.numeric(match(color, colors))*offset_increment
       line_end <- 2
       line_width=base_line_width
   }

   ymax=max(chr_lengths)-sigRegions$start[row]
   ymin = max(chr_lengths)-sigRegions$stop[row]

   #decide if plotting on left or right
   #r = rbinom(1,1,0.5)

   # force plot on the left
   if (region_type == "inversions") {
       r = 1
       offset = 0.25
   }
   else if (region_type == "gaps") {
       r = 1
       offset = 0.25
   }
   else {
       r = 0
   }

   if (fixed==0) {
      #randomly determine line offset from the chromosome
      offset=runif(1,0.16,0.3)
   }

   if (r == 0) {  #on left
      if (lines==1) {lines(x=c(chr-offset,chr-offset), y=c(ymax,ymin), lwd=line_width, col=color, lend=line_end)}
      if (lines==0) {text(x=(chr-offset), y=(ymax+ymin)/2, labels="*", cex=1.75, col=color)}
   }

   if (r == 1) {  #on right
      lines(x=c(chr+offset,chr+offset), y=c(ymax,ymin), lwd=line_width, col=color, lend=line_end)
      if (lines==0) {text(x=(chr+offset), y=(ymax+ymin)/2, labels="*", cex=1.75, col=color)}
   }
}


#kill the plotting device
dev.off()



