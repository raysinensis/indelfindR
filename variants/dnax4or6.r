#!/usr/bin/env Rscript

##measuring runtime
ptm <- proc.time()

##setup, read ab1 file for peaks
library(sangerseqR)
setwd("/home/rf/Desktop")
data1 <- read.abif("2-2.ab1")
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

##threshold maybe tuned for your needs
thresh=0.05
for (i in 1:length(loclist)){
	call1 <- 0
	max1 = max(Alist[i], Clist[i], Glist[i], Tlist[i])
	if (Alist[i] >= max1*thresh){
		call1 <- call1 + 1000}
	if (Clist[i] >= max1*thresh){
		call1 <- call1 + 100}
	if (Glist[i] >= max1*thresh){
		call1 <- call1 + 10}
	if (Tlist[i] >= max1*thresh){
		call1 <- call1 + 1}
	calllist <- append(calllist, call1)
	}

##the readout for each nt is a 4 digit number
sanger1 <- sangerseq(data1)
refseq <- primarySeq(sanger1, string = TRUE)

##locating a sequence before the cutsite (should be uniform), the actual compared sequence will be 30bp downstream
library(stringr)
site1 <- str_locate(refseq, "CGATCCTGCTTT")[2]
site31 <- site1 + 30
site40 <- site1 +39
refseqfrag <- unlist(calllist[site31:site40])
seqscore <- str_c(refseqfrag,collapse="")

##enter wildtype sequence around the area of interest here
(wtseq <- "GTCTACAGAGTGAGTTCTGGGAGAGCCAGAGCTACACAGAGAAACCCTGTCTTGGAGTGCGGGGAACAGAACTGTTGCAGTTGTGTTAACAGAGAACCCCTCGTGACGCTTGTTTTCTCAGATAACGGTTACTATCTGCCATACTACAAGAGAGAAAGGAATAAGCGGAGCACGCAGATCACAGTCAGGTTCCTGGACAGCCCCCACTACAGCAAGAACATCCGCAAGAAGGACCCTATCCTCCTGCTGCACTGGTGGAAGGAGATATTCGGGACGATCCTGCTTTGAATCGTGGCCACAACGTTTATCGTGCGCAGGCTTTTCCATCCTCAGCCCCACAGGGTAAGATGCTCTGTCAACCTAATGTGCTTCCAAGTGGTTGCTGTGTAGGAAACCTCTGGGAGGGAGAGTTCCGGCACTCAGACTTACCAGAGAATCATCTGCGCTTGTACTGTCTCCGTCACTTAGCCACAGAGTCACTGTACCTGAACAGCCAGTATGACATAGCCTACTGCTGCCAGGTTGCAGTCTCATACAGTGTGTGACTACAACATAGTATTTAGTATCTGTGTATCTGAAATAGATGAGGCACACTTAAAAAAAAAAAAAAAAAAAACACTAGCCACAATTGTTTTCAACTTTATTGCAATGTTACGGGTCGGCTGGAGAGATGGTTCAGTGGTTAGAGGACCTGAGTTCATTTTGCAGCACCCACATCAGACAACCCTTACCTGCCTGTAACTCCAGCTACTGGAAATCTGACA") 
sitewt <- str_locate(wtseq, "CGATCCTGCTTT")[2]

##investigating the -10-+9 range of indels
for (i in 0:20){
	fragname = paste("frag",i,sep="")
	assign(fragname, unlist(strsplit(substr(wtseq, sitewt+20+i, sitewt+29+i),split="")))}
counter <- 0

##searching for up to 4 seqs to match the sequencing results, usually around 1-2 min
for (i in 0:18){
	fi <- get(paste("frag",i,sep=""))
	for (j in i:19){
		fj <- get(paste("frag",j,sep=""))
		for (k in j:20){
			fk <- get(paste("frag",k,sep=""))
			for (n in k:20){
				fn <- get(paste("frag",n,sep=""))
						temp = c()
						vartempscore = c()
						for (m in 1:10){
							score <- 0
							temp <- append(temp, paste(fi[m],fj[m],fk[m],fn[m],sep=""))		
							if (str_detect(paste(fi[m],fj[m],fk[m],fn[m],sep=""),"A")){
								score <- score + 1000
								}
							if (str_detect(paste(fi[m],fj[m],fk[m],fn[m],sep=""),"C")){
								score <- score + 100
								}
							if (str_detect(paste(fi[m],fj[m],fk[m],fn[m],sep=""),"G")){
								score <- score + 10
								}
							if (str_detect(paste(fi[m],fj[m],fk[m],fn[m],sep=""),"T")){
								score <- score + 1
								}
							vartempscore <- append(vartempscore, score)		
							varname = paste("var",i,j,k,n,sep="a")
							varscore = paste("varscore",i,j,k,n,sep="a")
							assign(varname, temp)
							assign(varscore,vartempscore)
							}
						tempstr <- str_c(vartempscore,collapse="")
						##print(tempstr)
						if (tempstr == seqscore){
						print (varname)
						counter <- counter +1
						}
					}
				}
			}
		}
print("for 4 variables search:")
print(counter)

##if 4 seq failed, a much more intensive search of up to 6 seqs will start, 30min to 1hour
if (counter == 0){
for (i in 0:20){
	fragname = paste("frag",i,sep="")
	assign(fragname, unlist(strsplit(substr(wtseq, sitewt+20+i, sitewt+24+i),split="")))}
site31 <- site1 + 30
site35 <- site1 +34
refseqfrag <- unlist(calllist[site31:site35])
seqscore <- str_c(refseqfrag,collapse="")
for (i in 0:18){
	fi <- get(paste("frag",i,sep=""))
	for (j in i:19){
		fj <- get(paste("frag",j,sep=""))
		for (k in j:20){
			fk <- get(paste("frag",k,sep=""))
			for (n in k:20){
				fn <- get(paste("frag",n,sep=""))
				for (p in n:20){
					fp <- get(paste("frag",p,sep=""))
					for (q in p:20){
						fq <- get(paste("frag",q,sep=""))
						temp = c()
						vartempscore = c()
						for (m in 1:5){
							score <- 0
							temp <- append(temp, paste(fi[m],fj[m],fk[m],fn[m],fp[m],fq[m],sep=""))		
							if (str_detect(paste(fi[m],fj[m],fk[m],fn[m],fp[m],fq[m],sep=""),"A")){
								score <- score + 1000
								}
							if (str_detect(paste(fi[m],fj[m],fk[m],fn[m],fp[m],fq[m],sep=""),"C")){
								score <- score + 100
								}
							if (str_detect(paste(fi[m],fj[m],fk[m],fn[m],fp[m],fq[m],sep=""),"G")){
								score <- score + 10
								}
							if (str_detect(paste(fi[m],fj[m],fk[m],fn[m],fp[m],fq[m],sep=""),"T")){
								score <- score + 1
								}
							vartempscore <- append(vartempscore, score)		
							varname = paste("var",i,j,k,n,p,q,sep="a")
							varscore = paste("varscore",i,j,k,n,p,q,sep="a")
							assign(varname, temp)
							assign(varscore,vartempscore)
							}
						tempstr <- str_c(vartempscore,collapse="")
						##print(tempstr)
						if (tempstr == seqscore){
						print (varname)
						counter <- counter +1
						}
					}
				}
			}
		}
	}
}
print(counter)
}
proc.time() - ptm
