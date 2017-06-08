args<-commandArgs(trailingOnly = TRUE)
setwd("/opt/workhere/ChromFa")
system(paste("python ",args[1],".chrY.psr.py",sep=""))
setwd("/opt/")
