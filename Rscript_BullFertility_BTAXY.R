getwd()
setwd("C:/Users/hendy/Documents/OneDrive - University of Florida/Masters/Research/BullFertility/BTAXY/Autosomal - GenABEL/")


### First: do it in the server:
# module load R/3.2.0
# R
# library(GenABEL)
# data = load.gwaa.data(phenofile = "phenotype.txt", genofile = "outfile")                        #load files to gwaa.data object
# s = summary(data)   
# write.table(s, file = "Data", col.names=T,quote=F)

s <- read.table("Data", header = T)

## Load Results
load(file = "C:/Users/hendy/Documents/OneDrive - University of Florida/Masters/Research/BullFertility/BTAXY/Autosomal - GenABEL/teste/minFP/madd.rda")
load(file = "C:/Users/hendy/Documents/OneDrive - University of Florida/Masters/Research/BullFertility/BTAXY/Autosomal - GenABEL/teste/minFP/mdom.rda")
load(file = "mrec.rda")
load(file = "move.rda")
load(file = "mgen.rda")

## raw Pvalues
pvadd = madd[,4]
pvdom = mdom[,4]
#pvrec = mrec[,4]
#pvove = move[,4]

### DOM and REC - XY, ADD - PAR

#PAR_mdom <- as.data.frame(mdom[306182:309577,]); table(is.na(PAR_mdom[,4]))
#PAR_mrec <- as.data.frame(mrec[306182:309577,]); table(is.na(PAR_mrec[,4]))
PAR_madd <- as.data.frame(madd[306182:309577,]); table(is.na(PAR_madd[,4]))

X_mdom <- as.data.frame(mdom[309578:359331,]); table(is.na(X_mdom[,4])); table(is.na(X_mdom[,2]))
#X_mrec <- as.data.frame(mrec[309578:359331,]); table(is.na(X_mrec[,4])); table(is.na(X_mrec[,2]))

#Y_mdom <- as.data.frame(mdom[359332:360558,]); table(is.na(Y_mdom[,4])); table(is.na(Y_mdom[,2]))
#Y_mrec <- as.data.frame(mrec[359332:360558,]); table(is.na(Y_mrec[,4])); table(is.na(Y_mrec[,2]))




## Filtering Results - Autosomal (all), PAR - 30 (ADD), X and Y - 31 and 32 (DOM and REC):
fadd = PAR_madd[,2] < 0 | is.na(PAR_madd[,2]); table(fadd); table(is.na(PAR_madd[,4]))                     #fadd (filtro): se madd Chisq<0 or Chisq=NA  = TRUE (24103 SNP)
fdom = X_mdom[,2] < 0 | is.na(X_mdom[,2]); table(fdom); table(is.na(X_mdom[,4]))                     #table=summarize file (N_false N_true)
#frec = mrec[,2] < 0 | is.na(mrec[,2]); table(frec); table(is.na(mrec[,4]))                     #if Chisq<0 or NA P-value=NA, compare fadd and madd
#fove = move[,2] < 0 | is.na(move[,2]); table(fove); table(is.na(move[,4]))

#FIL_autosomal = (((fadd | fdom) | frec) | fove); table(FIL_autosomal)         
#FIL_PAR = fadd; table(FIL_PAR)
#FIL_XY = frec; table(FIL_XY)

##



## QQPlot 

par(mfrow=c(1,2))

qqplot(-log10(ppoints(length(PAR_madd[!fadd,4]))), -log10(PAR_madd[!fadd,4]),                             #qqplot(expected,observed)
       xlab = "Expected (-logP)", ylab = "Observed (-logP)", pch = 19, main = "ADD")
abline(0,1, lty = 3)

qqplot(-log10(ppoints(length(X_mdom[!fdom,4]))), -log10(X_mdom[!fdom,4]), 
       xlab = "Expected (-logP)", ylab = "Observed (-logP)", pch = 19, main = "DOM")
abline(0,1, lty = 3)

#qqplot(-log10(ppoints(length(mrec[!FIL,4]))), -log10(mrec[!FIL,4]), 
#xlab = "Expected (-logP)", ylab = "Observed (-logP)", pch = 19, main = "REC")
#abline(0,1, lty = 3)

#qqplot(-log10(ppoints(length(move[!FIL,4]))), -log10(move[!FIL,4]), 
#xlab = "Expected (-logP)", ylab = "Observed (-logP)", pch = 19, main = "OVE")
#abline(0,1, lty = 3)



### Plots 
#install.packages("tidyverse")
#install.packages("ggplot2")
#install.packages("rlang")
#install.packages("gridExtra")
#install.packages("cowplot")
#install.packages("colorspace")
#install.packages("ggsave")

library(ggplot2)
#library(colorspace)

#require(gridExtra)
#require(grid)
#require(cowplot)

##ADD - PAR

addPAR<- s[306182:309577,]; #addPAR$Chromosome[addPAR$Chromosome == 29] <- 30
SNP = data.frame(Name = rownames(addPAR), Chr = addPAR$Chromosome, Pos = addPAR$Position)
SNP$Chr = as.numeric(as.character(SNP$Chr))
#dim(SNP); SNPf = SNP[!fadd,];dim(SNPf)
SNP$Pval = -log10(PAR_madd[,4]) ## Additive 
SNP$within = seq(1:nrow(SNP))

SNP[order(-SNP$Pval),][1:20,]


#log10Pe <- expression(paste("Expected -log"[10], plain(P)))
#log10Po <- expression(paste("Observed -log"[10], plain(P)))
##

colo <- c("blue2")
c <- ggplot(SNP, aes(x = within, y = Pval, col = factor(Chr)))


c + scale_color_manual(values = colo) + 
  geom_point( ) + ylab("-Log10 P\n") + xlab("\nChromosome 30 \n Pseudoautosomal Region\n") +
  theme(axis.title.x = element_text(face = "bold", size = 16)) + 
  theme(axis.title.y = element_text(face = "bold", size = 16)) +
  theme(axis.text.x  = element_text(size = 12)) + 
  theme(axis.text.y  = element_text(size = 14)) +
  ylim(0,6)+ 
  theme_bw() +
  theme(legend.position = "none", panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.title.y = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=16),
        axis.text.x  = element_text(size=12),
        axis.text.y  = element_text(size=14))+
  geom_hline(yintercept = 2, lty = 5, size = 1, colour = "red") +
  tiff("GWAS-BullFertility-PAR.tiff", width = 1000, height = 1000)

dev.off()


##DOM - X specific chromosome

domX<- s[s$Chromosome==31,]
SNP = data.frame(Name = rownames(domX), Chr = domX$Chromosome, Pos = domX$Position)
SNP$Chr = as.numeric(as.character(SNP$Chr))
dim(SNP); SNPf = SNP[!fdom,];dim(SNPf)
SNPf$Pval = -log10(X_mdom[!fdom,4]) ## Dominance
SNPf$within = seq(1:nrow(SNPf))

SNPf[order(-SNPf$Pval),][1:20,]

colo <- c("blue2")
c <- ggplot(SNPf, aes(x = within, y = Pval, col = factor(Chr)))


c + scale_color_manual(values = colo) + 
  geom_point( ) + ylab("-Log10 P\n") + xlab("\nChromosome 31 \n BTAX-specific\n") +
  theme(axis.title.x = element_text(face = "bold", size = 16)) + 
  theme(axis.title.y = element_text(face = "bold", size = 16)) +
  theme(axis.text.x  = element_text(size = 12)) + 
  theme(axis.text.y  = element_text(size = 14)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.title.y = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=16),
        axis.text.x  = element_text(size=12),
        axis.text.y  = element_text(size=14))+
  geom_hline(yintercept = 2, lty = 5, size = 1, colour = "red") +
  tiff("GWAS-BullFertility-X.tiff", width = 1150, height = 550)

dev.off()


