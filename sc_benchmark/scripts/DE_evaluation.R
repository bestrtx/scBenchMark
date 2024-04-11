#差异表达评价模块
load("geneinfo.rds")
scRNA.markers <- FindMarkers(scRNA, ident.1=0) 
#scRNA.markers <- FindAllMarkers(scRNA)

upsample=scRNA.markers[scRNA.markers$avg_log2FC>0,]
downsample=scRNA.markers[scRNA.markers$avg_log2FC<0,]
upgene=geneinfo[geneinfo$DEFacGroup1>1,]$Gene
downgene=geneinfo[geneinfo$DEFacGroup1<1,]$Gene

#计算上调差异基因精度
TP=0
for (tmp in rownames(upsample)) {
  if (tmp %in% upgene){
    TP=TP+1
  }
}
FP=length(rownames(upsample))-TP
FN=length(upgene)-TP
TN=length(rownames(geneinfo))-FN-TP
Recall=TP/(TP+FN)
Accuracy=(TP+TN)/length(rownames(geneinfo))
Precision=TP/(TP+FP)
F1score=2*((Precision*Recall)/(Precision+Recall))

write.table("**DE_evaluation_up**",output_txt, row.names=F, quote=FALSE, col.names = FALSE, append = TRUE)
write.table("TP:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = "")
write.table(TP,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)
write.table("FP:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = "")
write.table(FP,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)
write.table("TN:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = "")
write.table(TN,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)
write.table("FN:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = "")
write.table(FN,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)
write.table("Accuracy:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = "")
write.table(Accuracy,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)
write.table("Recall:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = "")
write.table(Recall,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)
write.table("Precision:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = "")
write.table(Precision,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)
write.table("F1score:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = "")
write.table(F1score,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)


TP=0
for (tmp in rownames(downsample)) {
  if (tmp %in% downgene){
    TP=TP+1
  }
}
FP=length(rownames(downsample))-TP
FN=length(downgene)-TP
TN=length(rownames(geneinfo))-FN-TP
Recall=TP/(TP+FN)
Accuracy=(TP+TN)/length(rownames(geneinfo))
Precision=TP/(TP+FP)
F1score=2*((Precision*Recall)/(Precision+Recall))

write.table("**DE_evaluation_down**",output_txt, row.names=F, quote=FALSE, col.names = FALSE, append = TRUE)
write.table("TP:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = "")
write.table(TP,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)
write.table("FP:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = "")
write.table(FP,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)
write.table("TN:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = "")
write.table(TN,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)
write.table("FN:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = "")
write.table(FN,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)
write.table("Accuracy:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = "")
write.table(Accuracy,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)
write.table("Recall:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = "")
write.table(Recall,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)
write.table("Precision:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = "")
write.table(Precision,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)
write.table("F1score:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = "")
write.table(F1score,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)