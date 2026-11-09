#!/usr/bin/env Rscript

Sys.setenv(ACEM_BIN_BP="500")
source("RKO_TAZ_v_VAL_project/Final_figure_analysis/outputs/148_acEMseq_VAL_average_compact_waterfalls.R")

# Display VD-favored before Comparable within each H3K27me1 state.
display_order <- c("VD_gain", "Comparable_gain", "VD_maintain",
  "Comparable_maintain")
objects <- objects[display_order]
panel_labels <- panel_labels[display_order]

status_cols <- c("white", "black")

draw_status <- function(obj) {
  nr <- nrow(obj$signal)
  methylated <- as.integer(obj$meta$mean_VAL_5mC >= 0.50)
  z <- matrix(rev(methylated), nrow=1)
  image(1, seq_len(nr), z, col=status_cols,
    breaks=c(-0.5,0.5,1.5), axes=FALSE, xlab="", ylab="",
    useRaster=TRUE)
  box(col="black",lwd=.9)
  title(main=">=50%",cex.main=.69,line=.35)
}

draw_status_legend <- function() {
  plot.new(); plot.window(c(0,1),c(0,1))
  text(.05,.92,"Mean VAL 5mC",adj=c(0,.5),font=2,cex=.83)
  rect(.08,.67,.30,.76,col="black",border="black",lwd=.7)
  text(.38,.715,"Methylated\n(>=50%)",adj=c(0,.5),cex=.70)
  rect(.08,.47,.30,.56,col="white",border="black",lwd=.7)
  text(.38,.515,"Below 50%",adj=c(0,.5),cex=.70)
  yy <- seq(.08,.35,length.out=length(signal_colors)+1L)
  for(i in seq_along(signal_colors))
    rect(.08,yy[i],.30,yy[i+1],col=signal_colors[i],border=NA)
  rect(.08,.08,.30,.35,border="black",lwd=.7)
  text(.38,.08,"0",adj=c(0,.5),cex=.65)
  text(.38,.215,"0.5",adj=c(0,.5),cex=.65)
  text(.38,.35,"1",adj=c(0,.5),cex=.65)
  rect(.08,.005,.24,.045,col=missing_color,border="black",lwd=.5)
  text(.30,.025,"No data",adj=c(0,.5),cex=.62)
}

draw_all <- function() {
  # Each heatmap is followed by its aligned binary methylation-status strip.
  layout(matrix(1:9,nrow=1),widths=c(1,.095,1,.095,1,.095,1,.095,.45))
  par(oma=c(1.0,.5,4.3,.2),bg="white")
  for(i in seq_along(objects)) {
    par(mar=c(2.8,if(i==1)3.2 else .75,2.0,.15),bg="white")
    draw_heat(objects[[i]],panel_labels[i],show_y=i==1)
    par(mar=c(2.8,.05,2.0,.40),bg="white")
    draw_status(objects[[i]])
  }
  par(mar=c(2.8,.1,2.0,.1),bg="white")
  draw_status_legend()
  mtext("TAZ H3K27me1 gain",side=3,outer=TRUE,at=.205,line=1.4,font=2,cex=1.05)
  mtext("TAZ-maintained H3K27me1",side=3,outer=TRUE,at=.615,line=1.4,font=2,cex=1.05)
  mtext("VAL acEM-seq",side=3,outer=TRUE,line=3.15,font=2,cex=1.18)
  mtext("Position relative to ATAC peak center (kb)",side=1,outer=TRUE,
    line=-.10,cex=.86)
}

out_png <- file.path(outdir,
  "VAL_average_acEMseq_5mC_500bp_RANKED_HEATMAPS_WITH_50PCT_STATUS.png")
out_pdf <- file.path(outdir,
  "VAL_average_acEMseq_5mC_500bp_RANKED_HEATMAPS_WITH_50PCT_STATUS.pdf")
png(out_png,width=6000,height=3900,res=600,bg="white")
draw_all(); dev.off()
pdf(out_pdf,width=10,height=6.5,useDingbats=FALSE,bg="white")
draw_all(); dev.off()
cat("Saved:\n",out_png,"\n",out_pdf,"\n")
