#!/usr/bin/env Rscript

Sys.setenv(ACEM_BIN_BP="250")
source("RKO_TAZ_v_VAL_project/Final_figure_analysis/outputs/148_acEMseq_VAL_average_compact_waterfalls.R")
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

root <- "RKO_TAZ_v_VAL_project/Final_figure_analysis/"
outdir <- file.path(root,"outputs/acEMseq_5mC_waterfall_heatmaps",
  "VAL_REPLICATE_AVERAGE_250BP","FRAGMENT_METHYLATION_THRESHOLDS")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

eligible <- rbindlist(lapply(names(objects),function(nm){
  z <- as.data.table(objects[[nm]]$meta)
  z[,region_class:=nm]
  z
}),fill=TRUE)
eligible[,`:=`(
  H3K27me1_state=fifelse(grepl("gain$",region_class),
    "TAZ H3K27me1 gain","TAZ-maintained H3K27me1"),
  accessibility_class=fifelse(grepl("^VD",region_class),"VD-favored","Comparable")
)]
eligible[,H3K27me1_state:=factor(H3K27me1_state,
  c("TAZ H3K27me1 gain","TAZ-maintained H3K27me1"))]
eligible[,accessibility_class:=factor(accessibility_class,c("VD-favored","Comparable"))]

thresholds <- c(50,60,70)
calls <- rbindlist(lapply(thresholds,function(t){
  x <- copy(eligible)
  x[,`:=`(threshold_percent=t,
    methylation_call=fifelse(100*mean_VAL_5mC >= t,"Methylated","Below threshold"))]
  x
}))
calls[,threshold_label:=factor(paste0(threshold_percent,"% threshold"),
  paste0(thresholds,"% threshold"))]
calls[,methylation_call:=factor(methylation_call,c("Below threshold","Methylated"))]

counts <- calls[, .N, by=.(H3K27me1_state,accessibility_class,
  threshold_percent,threshold_label,methylation_call)]
counts[,total:=sum(N),by=.(H3K27me1_state,accessibility_class,threshold_percent)]
counts[,percent:=100*N/total]

tests <- calls[,{
  vd_yes <- sum(methylation_call=="Methylated" & accessibility_class=="VD-favored")
  vd_n <- sum(accessibility_class=="VD-favored")
  co_yes <- sum(methylation_call=="Methylated" & accessibility_class=="Comparable")
  co_n <- sum(accessibility_class=="Comparable")
  ft <- fisher.test(matrix(c(vd_yes,vd_n-vd_yes,co_yes,co_n-co_yes),
    nrow=2,byrow=TRUE))
  .(VD_methylated=vd_yes,VD_total=vd_n,Comparable_methylated=co_yes,
    Comparable_total=co_n,VD_percent=100*vd_yes/vd_n,
    Comparable_percent=100*co_yes/co_n,OR_VD_vs_Comparable=unname(ft$estimate),
    ci_low=ft$conf.int[1],ci_high=ft$conf.int[2],fisher_p=ft$p.value)
},by=.(H3K27me1_state,threshold_percent,threshold_label)]
tests[,BH_p:=p.adjust(fisher_p,"BH")]
tests[,label:=paste0("VD vs Comparable: ",sprintf("%.1f",VD_percent),"% vs ",
  sprintf("%.1f",Comparable_percent),"%\nOR = ",sprintf("%.2f",OR_VD_vs_Comparable),
  "; BH p ",fifelse(BH_p<.0001,"< 0.0001",paste0("= ",sprintf("%.3f",BH_p))))]

counts[,segment_label:=fifelse(percent < 15,
  paste0(format(N,big.mark=",")," (",sprintf("%.1f%%",percent),")"),
  paste0(format(N,big.mark=","),"\n",sprintf("%.1f%%",percent)))]
fills <- c("Below threshold"="white","Methylated"="black")

p <- ggplot(counts,aes(accessibility_class,percent,fill=methylation_call))+
  geom_col(width=.68,color="black",linewidth=.65)+
  geom_text(aes(label=segment_label),position=position_stack(vjust=.5),
    color=ifelse(counts$methylation_call=="Methylated","white","black"),
    fontface="bold",size=2.8,lineheight=.9)+
  geom_text(data=tests,aes(x=1.5,y=Inf,label=label),inherit.aes=FALSE,
    vjust=1.25,size=2.65,fontface="bold")+
  facet_grid(threshold_label~H3K27me1_state)+
  scale_fill_manual(values=fills)+
  scale_y_continuous(limits=c(0,125),breaks=seq(0,100,25),
    labels=function(x) paste0(x,"%"),expand=expansion(mult=c(0,0)))+
  labs(title="Methylation of sufficiently covered H3K27ac-associated fragments",
    subtitle="Mean VAL acEM-seq 5mC; eligible fragments have signal in >=6 of 12 250-bp bins",
    x=NULL,y="Sufficiently covered ATAC-linked fragments",fill=NULL)+
  theme_classic(base_size=11)+
  theme(plot.title=element_text(face="bold",size=14),legend.position="top",
    strip.background=element_rect(fill="white",color="black",linewidth=.8),
    strip.text=element_text(face="bold"),panel.border=element_rect(fill=NA,
      color="black",linewidth=.8),axis.line=element_blank(),
    axis.text.x=element_text(face="bold"))

fwrite(calls,file.path(outdir,"H3K27ac_FRAGMENT_METHYLATION_CALLS_50_60_70.tsv.gz"),sep="\t")
fwrite(counts,file.path(outdir,"H3K27ac_FRAGMENT_METHYLATION_COUNTS_50_60_70.tsv"),sep="\t")
fwrite(tests,file.path(outdir,"H3K27ac_FRAGMENT_METHYLATION_FISHER_TESTS_50_60_70.tsv"),sep="\t")
ggsave(file.path(outdir,"H3K27ac_FRAGMENT_METHYLATION_THRESHOLDS_50_60_70.png"),p,
  width=8.5,height=9,dpi=600,bg="white")
ggsave(file.path(outdir,"H3K27ac_FRAGMENT_METHYLATION_THRESHOLDS_50_60_70.pdf"),p,
  width=8.5,height=9,useDingbats=FALSE,bg="white")
print(counts);print(tests)
