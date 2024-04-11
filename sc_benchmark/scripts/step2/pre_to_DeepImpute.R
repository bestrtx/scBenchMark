source("R_conf.R")
#######################################################
#DeepImpute
load(datapath)
dat=as.matrix(dat)
dat=t(dat)
write.table(dat,"./sc_benchmark/output/step2/deepImpute/DeepImpute_log.csv",row.names=TRUE,col.names=TRUE,sep=",")

system("bash ./sc_benchmark/scripts/start_monitor.sh -n ./sc_benchmark/output/step2/deepImpute/usage -t 5")

#命令行处理

#deepImpute --cores 50 ./output/step2/deepImpute/DeepImpute_log.csv -o ./output/step2/deepImpute/imputed.csv