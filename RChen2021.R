# RChen for RT-qPCR compute
# Version 2021.6.13
# including <RChen(x,y)>
#执行代码:RChen(主文件名.CSV,扩增效率文件.CSV,内参基因基因名,对照样本样本名)
#多基因(参考详细结果)，多样本，Pfaffl法，单德尔塔法，双德尔塔法，单内参，相对表达计算，结果返回所有重复样品的数据便于后期SPSS分析

RChen=function(X,E,ref,nc){
cat("Welcome to RChen\nRChen★\nCurrent Version: 2021.6.13\n")
cat("Start working\n")
address=getwd()
E=deparse(substitute(E));e=data.frame();e=read.csv(E)
X=deparse(substitute(X));x=data.frame();x=read.csv(X)
cat(paste("成功读取文件",X,"\n"))
ref=deparse(substitute(ref))
colnames(x)=c("samples","name","ct")
colnames(e)=c("enames","effective")
cat("数据表头重置成功\n")
refmean=c()
for(i in 1:nlevels(factor(x$samples))){
refmean=c(refmean,mean(subset(subset(x,x$name==ref),subset(x,x$name==ref)$samples==levels(factor(x$samples))[i])$ct))}
reftab1=data.frame()
reftab1=data.frame(refname=levels(factor(x$samples)),refmeans=refmean)
cat(paste("样本内参基因",ref,"均值计算成功\n"))
x1=data.frame(); x1=subset(x,x$name!=ref);x1$expresslevel=0
eref=c(subset(e,e$enames==ref)$effective)
egene=data.frame(); egene=subset(e,e$enames!=ref)
for(i in 1:nlevels(factor(egene$enames))){
for(j in 1:nlevels(factor(reftab1$refname))){
for (z in 1:nrow(x1)){
if(x1$name[z]==egene$enames[i]){
if(x1$samples[z]==reftab1$refname[j]){
x1$expresslevel[z]=((1+eref)^reftab1$refmeans[j])/((1+egene$effective[i])^(x1$ct[z])) 
# core formula in this program
}}}}}; cat("基因表达数据标准化成功\n")
S=data.frame(); genemean=c()
nc=deparse(substitute(nc))
S=subset(x1,x1$samples==nc)
for(i in 1:nlevels(factor(S$name))){
genemean=c(genemean,mean(subset(S,S$name==levels(factor(S$name))[i])$expresslevel))}
genetab=data.frame()
genetab=data.frame(x=levels(factor(S$name)),y=genemean)
cat(paste("对照样本",nc,"基因均值计算成功\n")); x1$RE=0
for(i in 1:nlevels(factor(genetab$x))){
for(z in 1:nrow(x1)){
if(x1$name[z]==genetab$x[i]){
x1$RE[z]=x1$expresslevel[z]/genetab$y[i]
}}};cat("相对基因表达量计算成功并开始生成表格数据\n")
x2=data.frame()
x2=data.frame(SamplesNames=x1$samples,GeneNames=x1$name,RelativeExpression=x1$RE)
if (file.exists("./R-QPCR")==TRUE){cat("阁下 目标文件夹 R-QPCR 已存在\n")}else{
dir.create("./R-QPCR", recursive=TRUE)
cat("阁下 目标文件夹 R-QPCR 已创建\n")}
setwd( "./R-QPCR")
write.csv(x2,paste("R语言荧光定量PCR 双德尔塔 标准结果.csv",gsub(":","_",Sys.time()),".csv"),row.names=FALSE)
write.csv(x1,paste("R语言荧光定量PCR 详细数据结果.csv",gsub(":","_",Sys.time()),".csv"),row.names=FALSE)
cat(paste("数据已放入",getwd()," 文件夹内 请阁下查收\n",sep=""))
cat("标准数据显示如下：\n")
setwd(address);x2}
