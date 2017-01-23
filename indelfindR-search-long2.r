#!/usr/bin/env Rscript

##measuring runtime
ptm <- proc.time()

##setup
library(sangerseqR)
mainpath<-'/home/rf/Desktop/crispr mice sequencing/'
setwd(mainpath)
folderlist<-list.files(mainpath)

##read each file
for (fileab1 in folderlist){
write("-----------------", file = "sequences1.txt", sep = "\n", append = TRUE)
write(fileab1, file = "sequences1.txt", sep = "\n", append = TRUE)
data1 <- read.abif(fileab1)
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
calllist = c()

##the readout for each nt is a 4 digit number
sanger1 <- sangerseq(data1)
refseq <- primarySeq(sanger1, string = TRUE)

##locating a sequence before the cutsite (should be uniform), the actual compared sequence will be 55bp downstream
library(stringr)
site1 <- str_locate(refseq, "GGACCCTATCCT")[2]
site41 <- site1 + 40
site90 <- site1 + 89

##handling long deletions, other errors
if (is.na(site1)) {
	write("notfound", file = "sequences1.txt", sep = "\n", append = TRUE)	
	next
	}
if (length(Alist)<(site90+20)) {
	write("outofrange", file = "sequences1.txt", sep = "\n", append = TRUE)	
	next}
site41 <- site41+55
site90 <- site90+55

##threshold maybe tuned for your needs
##autothreshold can be ran
maxA = max(Alist[site1:(site90+20)])
maxC = max(Clist[site1:(site90+20)])
maxG = max(Glist[site1:(site90+20)])
maxT = max(Tlist[site1:(site90+20)])
signal1 = c()
for (i in site1:(site90+20)){
	signal1 = rbind(signal1, (c(Alist[i]/maxA, Clist[i]/maxC, Glist[i]/maxG, Tlist[i]/maxT)))
	}
signal2 = sort(signal1)
dropmax = 0
cutnum = 0
for (i in 1:(length(signal2)-1)){
	if (signal2[i]>0.01 & signal2[i]<0.25){
		drop1 = signal2[i+1]/signal2[i]
		if (drop1 > dropmax){
			cutnum = signal2[i+1]
			dropmax = drop1
			}
		}
	}
print(cutnum)
##or manually set
##cutnum=0.05
thresh1 <- c(cutnum*maxA,cutnum*maxC,cutnum*maxG,cutnum*maxT)
print(thresh1)
for (i in 1:length(loclist)){
	call1 <- 0
	if (Alist[i] >= thresh1[1]){
		call1 <- call1 + 1000}
	if (Clist[i] >= thresh1[2]){
		call1 <- call1 + 100}
	if (Glist[i] >= thresh1[3]){
		call1 <- call1 + 10}
	if (Tlist[i] >= thresh1[4]){
		call1 <- call1 + 1}
	calllist <- append(calllist, call1)
	}
diglist <- sprintf("%04d", calllist[site41:site90])
refseqfrag <- unlist(diglist)
seqscore <- str_c(refseqfrag,collapse="")

##enter wildtype sequence around the area of interest here
(wtseq <- "GTCTACAGAGTGAGTTCTGGGAGAGCCAGAGCTACACAGAGAAACCCTGTCTTGGAGTGCGGGGAACAGAACTGTTGCAGTTGTGTTAACAGAGAACCCCTCGTGACGCTTGTTTTCTCAGATAACGGTTACTATCTGCCATACTACAAGAGAGAAAGGAATAAGCGGAGCACGCAGATCACAGTCAGGTTCCTGGACAGCCCCCACTACAGCAAGAACATCCGCAAGAAGGACCCTATCCTCCTGCTGCACTGGTGGAAGGAGATATTCGGGACGATCCTGCTTTGAATCGTGGCCACAACGTTTATCGTGCGCAGGCTTTTCCATCCTCAGCCCCACAGGGTAAGATGCTCTGTCAACCTAATGTGCTTCCAAGTGGTTGCTGTGTAGGAAACCTCTGGGAGGGAGAGTTCCGGCACTCAGACTTACCAGAGAATCATCTGCGCTTGTACTGTCTCCGTCACTTAGCCACAGAGTCACTGTACCTGAACAGCCAGTATGACATAGCCTACTGCTGCCAGGTTGCAGTCTCATACAGTGTGTGACTACAACATAGTATTTAGTATCTGTGTATCTGAAATAGATGAGGCACACTTAAAAAAAAAAAAAAAAAAAACACTAGCCACAATTGTTTTCAACTTTATTGCAATGTTACGGGTCGGCTGGAGAGATGGTTCAGTGGTTAGAGGACCTGAGTTCATTTTGCAGCACCCACATCAGACAACCCTTACCTGCCTGTAACTCCAGCTACTGGAAATCTGACA") 
sitewt <- str_locate(wtseq, "GGACCCTATCCT")[2]
sitewt <- sitewt+55

##investigating the -35:+300 range of indels
wtscorelist = c()
for (i in (-35:300)){
	fragname = paste("frag",i,sep="")
	fragseq = str_c(unlist(strsplit(substr(wtseq, sitewt+40+i, sitewt+89+i),split="")),collapse="")
	wtscorelist = rbind(wtscorelist,c(fragname,fragseq))
	for (j in 1:10){
		matchval = 1
		readbase = calllist[site41+j-1]
		if ((substring(fragseq,j,j)=="A") & (readbase%/%1000%%10==1)){
		matchval = 0}
		if ((substring(fragseq,j,j)=="C") & (readbase%/%100%%10==1)){
		matchval = 0}
		if ((substring(fragseq,j,j)=="G") & (readbase%/%10%%10==1)){
		matchval = 0}
		if ((substring(fragseq,j,j)=="T") & (readbase%/%1%%10==1)){
		matchval = 0}
		if (matchval == 1){
		break}
	}
	if (matchval == 0){
	print(fragname)
	write(fragname, file = "sequences1.txt", sep = "\n", append = TRUE)
}}}
proc.time() - ptm

