#!/usr/bin/env Rscript

##measuring runtime
ptm <- proc.time()

##setup, read ab1 file for peaks
library(sangerseqR)
setwd("/home/rf/Desktop")
data1 <- read.abif("2_pnpt1F.ab1")
loclist <- data1@data$PLOC.2
Alist = c()
Clist = c()
Glist = c()
Tlist = c()
for (locnum in loclist){
	Glist <- append(Glist, (data1@data$DATA.9[locnum]))
	Alist <- append(Alist, (data1@data$DATA.10[locnum]))
	Tlist <- append(Tlist, (data1@data$DATA.11[locnum]))
	Clist <- append(Clist, (data1@data$DATA.12[locnum]))}

##detect peaks
calllist <- matrix(ncol = 4)
dimnames(calllist) = list(c(),c("A","C","G","T"))
for (i in 1:length(loclist)){
	max1 = max(Alist[i], Clist[i], Glist[i], Tlist[i])
	calllist <- rbind(calllist, c(Alist[i], Clist[i], Glist[i], Tlist[i]))
	if (max1 == 0){
		break
		}
	}
calllist <- calllist[-1,]
calllist <- calllist[-(length(calllist)/4),]
##trim 15nt from each end
calllist <- calllist[-(1:15),]
calllist <- calllist[-((length(calllist)/4-14):(length(calllist)/4)),]

##getting highest and 2nd highest reads, above threshold
call1 = c()
call2 = c()
thresauto = 0
for (i in 1:(length(calllist)/4)){
	max1 = max(calllist[i,])
	max2 = max(calllist[i,][calllist[i,]!=max1])
	if (max2 == 0){
		max3 = 0}
	else {
		max3 = max(calllist[i,][calllist[i,]<max2])
		if ((max3/max1)>thresauto){
			thresauto = (max3/max1)
			}
		}
	}
##thresh=0.2
thresh=thresauto
for (i in 1:(length(calllist)/4)){
	max1 = max(calllist[i,])
	max2 = max(calllist[i,][calllist[i,]!=max1])
	if (max1 == calllist[i,1]){
		call1 <- append(call1,"A")
		}
	else if (max1 == calllist[i,2]){
		call1 <- append(call1,"C")
		}
	else if (max1 == calllist[i,3]){
		call1 <- append(call1,"G")
		}
	else if (max1 == calllist[i,4]){
		call1 <- append(call1,"T")
		}
	if (max2 <= max1*thresh){
		call2 <- append(call2, call1[i])
		}
	else {
	cat ("secondary call made at ", i, "\n")
	if (max2 == calllist[i,4]){
		call2 <- append(call2,"T")
		}
	else if (max2 == calllist[i,3]){
		call2 <- append(call2,"G")	
		}
	else if (max2 == calllist[i,2]){
		call2 <- append(call2,"C")
		}
	else if (max2 == calllist[i,1]){
		call2 <- append(call2,"A")
		}
	}
	}

##output files to txt
library(stringr)
call1txt <- str_c(call1,collapse="")
call2txt <- str_c(call2,collapse="")

write(call1txt, file = "2_pnpt1F.txt", sep = "\n")
write(call2txt, file = "2_pnpt1F.txt", append = TRUE)

proc.time() - ptm
