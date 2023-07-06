getwd() #"E:/EPIC data/HASM extended"

# place the below files in the working directory:
# Cross reactive probes on EPIC array.csv
#list_Probes_overlapping_genetic_variants_at_targeted_CpG_sites.csv
#Probes_overlapping_genetic_variants_at_single_base_extension_sites_for_Infinium_Type_I_probes.csv
#

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sva")
require("sva")
BiocManager::install("minfi")
require("minfi")
BiocManager::install("minfiData")
require("minfiData")

#place idats and saample sheet in the minfiData folder
baseDir <- system.file("CE idats",package="minfiData")
targets <- read.metharray.sheet(baseDir)
targets
RGSet <- read.metharray.exp(base = baseDir, targets = targets, extended=TRUE) #dim: 1051943 48 

RGSet@annotation=c(array='IlluminaHumanMethylationEPIC', annotation='ilm10b2.hg19') ## fixes problem with read.methyarray not assiging the correct annotation
RGSet ##1051943 48 
save(RGSet, file=paste0('RGSet.rda'))

BiocManager::install("IlluminaHumanMethylationEPICmanifest")
require("IlluminaHumanMethylationEPICmanifest")


phenoData<-pData(RGSet) ##generate a PhenoData object
phenoData
head(phenoData)

#initial QC at the level of methylated and unmethylated signal ratio
MSet<-preprocessRaw(RGSet) # go from probe based intensity values to methylation level values first
qc<-getQC(MSet)
head(qc)
plotQC(qc) ## saves as PDF to working directory -  all samples good
save(MSet, file=paste0('MSet.rda'))


rownames(RGSet) ## identify column names to label MDS clusters by
sex <- pData(RGSet)$Gender
all_diag<-pData(RGSet)$All_Diag
array_pos<- pData(RGSet)$Array
slide<-pData(RGSet)$Slide


mdsPlot(RGSet, numPositions = 1000, sampNames = names, sampGroups = sex, pch = 16, legendPos = "bottomleft", main = NULL)
mdsPlot(RGSet, numPositions = 1000, sampNames = names, sampGroups = all_diag, pch = 16, legendPos = "bottomleft", main = NULL)
mdsPlot(RGSet, numPositions = 1000, sampNames = names, sampGroups = array_pos, pch = 16, legendPos = "topright", main = NULL)
mdsPlot(RGSet, numPositions = 1000, sampNames = names, sampGroups = slide, pch = 16, legendPos = "topright", main = NULL)


####only thing that clearly seperates is Gender, do again when XY removed

###remove any probes that are not represented by 3 beads in more than 1 sample

###Number of beads
NBeads <- getNBeads(RGSet) ##provides number of beads with each probe on in each sample
head(NBeads )
testfailed<-NBeads<3
head(testfailed) ## TRUE/FALSE statement for whether <3 beads are present in each sample for each probe
sampleNames(RGSet)
test_remove <- names(which(rowMeans(testfailed)>0.025, TRUE)) ##probes with <3 beads in 0.025 (2.5%) = 1 sample 
length(test_remove) #3237 = number of probes to remove

head(rownames(RGSet))
RGSetbB<- RGSet[!rownames(RGSet)%in%test_remove,] ## removes the probes identified above -  1048706 48 - 3237 removed - good
save(RGSetbB, file=paste0('RGSetbB.rda'))


###Need to normalise at this point to keep control probes as loose them once convert to beta values (NB funnorm includes noob as a first step)
##Correlated reps first
sampleNames(RGSetbB)
rownames(RGSetbB)
colnames(RGSetbB)

##reps
#AZAC12       201496850012_R08C01
#AZAC12_rep   201496850197_R08C01
#AZAD01       201496710108_R02C01
#AZAD01_rep   201496850197_R02C01
#AZCC26       201492570151_R02C01
#AZCC26_rep   201496850012_R01C01
#MMP1_H01     201492570151_R07C01
#MMP1_H01_rep 201496850031_R07C01
#MMP1_H04     201496850031_R03C01
#MMP1_H04_rep 201496860001_R03C01


C12<-c("201496850012_R08C01","201496850197_R08C01")
head(getBeta(RGSetbB))
Betas<-getBeta(RGSetbB)


C12b<-(Betas)[,colnames(Betas)%in%C12]
head(C12)
corC12b<- cor(C12b, use="complete") ###  0.997339
plot(C12b)

D01<-c("201496710108_R02C01", "201496850197_R02C01")
D01b<-(Betas)[,colnames(Betas)%in%D01]
corD01b<- cor(D01b, use="complete") ###  0.9969911
corD01b
plot(D01b)

C26<-c("201492570151_R02C01", "201496850012_R01C01")
C26b<-(Betas)[,colnames(Betas)%in%C26]
corC26b<- cor(C26b, use="complete")  ##0.9967819
corC26b
plot(C26b)

H01<-c("201492570151_R07C01", "201496850031_R07C01")
H01b<-(Betas)[,colnames(Betas)%in%H01]
corH01b<- cor(H01b, use="complete") ###  0.9980647
corH01b
plot(H01b)

H04<-c("201496850031_R03C01", "201496860001_R03C01")
H04b<-(Betas)[,colnames(Betas)%in%H04]
corH04b<- cor(H04b, use="complete") ###  0.9972171
corH04b
plot(H04b)


###funnorm
fN<-preprocessFunnorm(RGSetbB)
fNbetas<-getBeta(fN)
C12bfN<-(fNbetas)[,colnames(fNbetas)%in%C12]
head(C12bfN)
corC12bfN<- cor(C12bfN, use="complete") ### 0.9961286, fN =  0.9970829
corC12bfN

D01bfN<-(fNbetas)[,colnames(fNbetas)%in%D01]
corD01bfN<- cor(D01bfN, use="complete") ### 0.9973658, fN =  0.9980659
corD01bfN

C26bfN<-(fNbetas)[,colnames(fNbetas)%in%C26]
corC26bfN<- cor(C26bfN, use="complete") ### 0.9973658, fN =  0.997934
corC26bfN

H01bfN<-(fNbetas)[,colnames(fNbetas)%in%H01]
corH01bfN<- cor(H01bfN, use="complete") ###  0.996143, fN = 0.998445
corH01bfN

H04bfN<-(fNbetas)[,colnames(fNbetas)%in%H04]
corH04bfN<- cor(H04bfN, use="complete") ###  0.996143, fN =  0.9982884
corH04bfN


### All correlations better with funNorm data
save(fN, file=paste0('fN (funNorm preprocessing done.rda')) ##sample names are barcodes
sampleNames(fN)
fN$Sample_ID
sampleNames(fN)<-fN$Sample_ID ## make sample names the sample ID
save(fN, file=paste0(CE_dir, 'fN (funNorm preprocessing done.rda'))
sampleNames(fN)


#### removed XY chromosomes
##XY on fN
fN.xybB <- fN[getAnnotation(fN)$chr%in%c("chrX","chrY"),] ##identify probes annotated to the x and y chromosomes
fN.xybB #19571 48 still, good
save(fN.xybB,file=paste(" fN(XYonly).RData",sep=""))
Betas.xybB<-getBeta(fN.xybB)
head(Betas.xybB)

##  XYheatmap to check samples are correct
by<- maPalette(low="blue", high="yellow") ##marray
xybasecolors <- c("red","blue","green","yellow","purple")


xycol <- xybasecolors[as.numeric(factor(fN$Gender))]
heatmap.2(cor(Betas.xybB),trace="none",Rowv=T,Colv=T,dendrogram="both",col=by,labRow=NA,ColSideColor=xycol)
legend("right",legend=levels(factor(gset.xybB$Gender)),fill=xybasecolors[1:nlevels(factor(gset.xybB$Gender))])
#all good 

### detection p value identification
lumi_dpvalbB<- detectionP(RGSetbB, type = "m+u") ## generates detection P value
head(lumi_dpvalbB)
lumi_failedbB<- lumi_dpvalbB > 0.01 ## identity probes with a detection p >0.01
head(lumi_failedbB)
lumi_dpval_removebB <- names(which(rowMeans(lumi_failedbB)>0.05, TRUE)) # identify probes with a detection p>0.01 in >0.05 (5% = 2 samples) = 2180 probes
head(lumi_dpval_removebB )

##remove failed probes from fN
fNbBbp<- fN[!rownames(fN)%in%lumi_dpval_removebB,] ####861221 48, correct amount removed good
save(fNbBbp,file=paste0("fNbBbp(bad beads and bad detection p removed).rda"))

##remove XY from fN
XYprobes<-rownames(fN.xybB) ## identify probes annotated to the XY probes
fNbBbpnXY<- fNbBbp[!rownames(fNbBbp)%in%XYprobes,] ##842089 48 (19132 removed, total XY should be 19571 but some must have been bad detection p too) ### still same for fN
save(fNbBbpnXY,file=paste0("fNbBbpnXY(bad beads, bad detection p, no XY removed).rda"))

####now just need to remove all the multiple binders etc as provided by https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1

install.packages("readr")
library("readr", lib.loc="~/R/win-library/3.4")


Cross_reactive_probes_on_EPIC_array <- read_csv("Cross reactive probes on EPIC array.csv")
View(list_Probes_overlapping_genetic_variants_at_targeted_CpG_sites)
View(Probes_overlapping_genetic_variants_at_single_base_extension_sites_for_Infinium_Type_I_probes)

head(row.names(list_Probes_overlapping_genetic_variants_at_targeted_CpG_sites))
head(Probes_overlapping_genetic_variants_at_single_base_extension_sites_for_Infinium_Type_I_probes)

target <- list_Probes_overlapping_genetic_variants_at_targeted_CpG_sites$Probe # probes overlapping genetic variants at CpG sites to be profiled  = 12510
sbe<-Probes_overlapping_genetic_variants_at_single_base_extension_sites_for_Infinium_Type_I_probes$PROBE #probes with genetic variants at the extension site required for detection = 414

fNbBbpnXYnt<- fNbBbpnXY[!rownames(fNbBbpnXY)%in%target,] ###remove target, now   830615 48  (11474 removed 12510 in list but some may have alreaDY been removed) same for fN
save(fNbBbpnXYnt,file=paste0("fNbBbpnXYnt(no SNPs at targeted CpG).rda"))

fNbBbpnXYntns<- fNbBbpnXYnt[!rownames(fNbBbpnXYnt)%in%sbe,] ## remove SBE, now 830328 48  (287 removed) 
save(fNbBbpnXYntns,file=paste0("fNbBbpnXYntns(no SNPs at targeted CpG or single base extension for tI).rda"))

head(Cross_reactive_probes_on_EPIC_array)
CRP<- Cross_reactive_probes_on_EPIC_array$X1 #43254
fNbBbpnXYntnsnc<-fNbBbpnXYntns[!rownames(fNbBbpnXYntns)%in%CRP,] ## remove cross reactive probes, 789701 48 , 40627 removed, same for fN

save(fNbBbpnXYntnsnc,file=paste0("fNbBbpnXYntnsnc(no SNPs at targeted CpG or single base extension for tI or cross reactive).rda"))
HEQCfN<-fNbBbpnXYntnsnc ##789701 48
save(HEQCfN,file=paste0("HEQCfN(QC'ed and funNorm normalised).rda")) ##GenomicRatioSet

mdsPlot(test, numPositions = 10000, sampGroups = groups, pch = 16, legendPos = "bottomleft", main = NULL)
#no longer clusters by sex

####combat to remove batch effects due to slide and position

combat<- HEQCfN ### make a new project for combat so don't mess up original


# Combat RUN --------------------------------------------------------------

#first remove cfDNA and fibroblasts samples that were part of a test for a different project
sampleNames(combat)
nonHASM<-c("8_MMP1_HO4","F___AZAD03","L12","L17")   

colnames(combat)
combat<-combat[,!colnames(combat)%in%nonHASM] ##789701 44 
sampleNames(combat)## super

M1<-getM(combat) #gets M values for adjustment (more statistically relavant version of beta value)
head(M1)
D1<-pData(combat) #gets pData
head(D1)

model<- model.matrix(~ All_Diag, data=D1) ##generate a model matrix telling combat what variation we want to MAINTAIN - in this case that is disease
bat<- ComBat(M1, D1$slide, model) ## perform combat on the M values, for slide, using the model to maintain variation of interest

# Error:rror in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
#contrasts can be applied only to factors with 2 or more levels 
#suggests some NAs in data


sum(is.na(M1))#0
sum(is.infinite(M1)) #88, need to remove these as combat won't work with them in (Error in while (change > conv) { : missing value where TRUE/FALSE needed)
sum(is.nan(M1)) #0
M1[!is.finite(M1)] <-NA ## replaces infinate with NA
colMeans(test, na.rm=TRUE)
sum(is.na(M1)) ## 88
M1<-na.omit(M1) ## removes M1
sum(is.na(M1))## 0
dim(rownames(M1))
nrow(M1) ###now 789613, was 789701

bat<- ComBat(M1, D1$Slide, model)

head(bat)
head(M1)


#check reps to see if it has made a positive difference

# uses reps to check correlations before and after each combat run
repC12<-c("AZAC12", "AZAC12_rep")
repD01<-c("AZAD01", "AZAD01_rep")
repC26<-c("AZCC26", "AZCC26_rep")
repH01<-c("MMP1_H01", "MMP1_H01_rep")
repH04<-c("MMP1_H04", "MMP1_H04_rep")

#AZAC12       201496850012_R08C01
#AZAC12_rep   201496850197_R08C01
#AZAD01       201496710108_R02C01
#AZAD01_rep   201496850197_R02C01
#AZCC26       201492570151_R02C01
#AZCC26_rep   201496850012_R01C01
#MMP1_H01     201492570151_R07C01
#MMP1_H01_rep 201496850031_R07C01
#MMP1_H04     201496850031_R03C01
#MMP1_H04_rep 201496860001_R03C01

repC12M1<-(M1)[,colnames(M1)%in%repC12]
corC12M1<- cor(repC12M1) ##0.9956776
repC12C1<-(bat)[,colnames(bat)%in%repC12]
corC12C1<- cor(repC12C1) ##0.9936973

repD01M1<-(M1)[,colnames(M1)%in%repD01]
corD01M1<- cor(repD01M1)### 0.9970861 
corD01M1### 0.9970861
repD01C1<-(bat)[,colnames(bat)%in%repD01]
corD01C1<- cor(repD01C1) ### 0.9946811
corD01C1

repC26M1<-(M1)[,colnames(M1)%in%repC26]
corC26M1<- cor(repC26M1)### 0.9968709 
corC26M1### 0.9968709 
repC26C1<-(bat)[,colnames(bat)%in%repC26]
corC26C1<- cor(repC26C1) ###  0.9941554
corC26C1

repH01M1<-(M1)[,colnames(M1)%in%repH01]
corH01M1<- cor(repH01M1)### 0.997348
corH01M1# 0.997348
repH01C1<-(bat)[,colnames(bat)%in%repH01]
corH01C1<- cor(repH01C1) ###  0.9954983
corH01C1

repH04M1<-(M1)[,colnames(M1)%in%repH04]
corH04M1<- cor(repH04M1)### 0.9972431
corH04M1### 0.9972431
repH04C1<-(bat)[,colnames(bat)%in%repH04]
corH04C1<- cor(repH04C1) ###    0.9920462
corH04C1

##### consistently got worse. Leav slide and mve to position


# Combat by sentrix location ----------------------------------------------

###### now repeat combat on run normalise on chip
## continue to use original combat, so M1 and D1

bat3<- ComBat(M1, D1$Array, model)
head(bat3)

corC12M1 ##0.9956776
repC12C3<-(bat3)[,colnames(bat3)%in%repC12]
corC12C3<- cor(repC12C3)
corC12C3 ##0.9976965


corD01M1### 0.9970861
repD01C3<-(bat3)[,colnames(bat3)%in%repD01]
corD01C3<- cor(repD01C3) 
corD01C3 ## 0.9977471

corC26M1### 0.9968709 
repC26C3<-(bat3)[,colnames(bat3)%in%repC26]
corC26C3<- cor(repC26C3) 
corC26C3 #0.9926786 ## worse

corH01M1# 0.997348
repH01C3<-(bat3)[,colnames(bat3)%in%repH01]
corH01C3<- cor(repH01C3)
corH01C3 ##0.9972513 ## worse

corH04M1### 0.9972431
repH04C3<-(bat3)[,colnames(bat3)%in%repH04]
corH04C3<- cor(repH04C3) 
corH04C3 # 0.9972953 ## better


####3 of 5 improved so go with and make new GRS
pDat<-pData(combat)
sampleNames(combat)
colnames(bat3)
combat2<-makeGenomicRatioSetFromMatrix(bat3, pData=pDat, array='IlluminaHumanMethylationEPIC', annotation='ilm10b2.hg19', what="M")

save(combat2,file=paste0("Combat by Array.rda"))

###### correcting by position improved the mean correlation between most replicates
##### use this corrected data


head(getBeta(combat))
head(getBeta(combat2)) #### betas of combat and combat2 are different so changing the M values is also changing beta values - good

HEpc<- combat2##789613 44
save(HEpc,file=paste0("HEpc (HE post QC, post normalisation and post combat_ ready for analysis.rda"))
###need to remove reps and irrelavant samples before analysis!

pDat<- pData(HEpc)
write.csv(pDat,file="pDat all samples HASMe.txt")


# Gene expression import and normalisation --------------------------------

###Gene Expression Data Analysis Human Gene 2.1 ST chip###
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("oligo")
require(oligo)
BiocManager::install("hugene21sttranscriptcluster.db")
require(hugene21sttranscriptcluster.db)
BiocManager::install("pd.hugene.2.1.st")
require(pd.hugene.2.1.st)

#library(pdInfoBuilder)

### place .cel files and sample sheet in the working directory
celList <- list.celfiles(getwd())
rawData <- read.celfiles(celList, pkgname='pd.hugene.2.1.st')

samplesheet<-read.csv("Affy samaple sheet for Tillie 240817 2.csv",header=TRUE,sep=",")
head(samplesheet)
colnames(samplesheet)
head(samplesheet)
rownames(samplesheet) ### this are numbers, make them the Cel.file.ID
rownames(samplesheet)<-samplesheet$Cel.file.ID
head(samplesheet)
sampleNames(rawData) ### in a different order to sample sheet
samplesheet_order<-samplesheet[sampleNames(rawData),]
head(samplesheet_order)
rownames(samplesheet_order)
rownames(samplesheet_order)==sampleNames(rawData) ### check all sampleNames are in the same order


sampleNames(rawData)<-samplesheet_order$Cel.file.ID.2
rownames(samplesheet_order)<-sampleNames(rawData)

hist(rawData)

normData <- rma(rawData)
#Background correcting
#Normalizing
#Calculating Expression

save(normData,file=paste("normData.RData",sep=""))
head(normData)

write.exprs(normData,file="log2 transformed and normalized data 240817.txt")

head(exprs(normData))

Annot <- data.frame(ACCNUM=sapply(contents(hugene21sttranscriptclusterACCNUM), paste, collapse=", "),
                    SYMBOL=sapply(contents(hugene21sttranscriptclusterSYMBOL), paste, collapse=", "), 
                    DESC=sapply(contents(hugene21sttranscriptclusterGENENAME), paste, collapse=", "),
                    CHR=sapply(contents(hugene21sttranscriptclusterCHR), paste, collapse=", "),
                    CHRLOC=sapply(contents(hugene21sttranscriptclusterCHRLOC), paste, collapse=", "),
                    ENTREZID=sapply(contents(hugene21sttranscriptclusterENTREZID), paste, collapse=", "),
                    GO=sapply(contents(hugene21sttranscriptclusterGO), paste, collapse=", "),
                    PATH=sapply(contents(hugene21sttranscriptclusterPATH), paste, collapse=", "),
                    REFSEQ=sapply(contents(hugene21sttranscriptclusterREFSEQ), paste, collapse=", "))

Annot <- data.frame(ACCNUM=sapply(as.list(hugene21sttranscriptclusterACCNUM), paste, collapse=", "),
                    SYMBOL=sapply(as.list(hugene21sttranscriptclusterSYMBOL), paste, collapse=", "), 
                    DESC=sapply(as.list(hugene21sttranscriptclusterGENENAME), paste, collapse=", "),
                    CHR=sapply(as.list(hugene21sttranscriptclusterCHR), paste, collapse=", "),
                    CHRLOC=sapply(as.list(hugene21sttranscriptclusterCHRLOC), paste, collapse=", "),
                    ENTREZID=sapply(as.list(hugene21sttranscriptclusterENTREZID), paste, collapse=", "),
                    GO=sapply(as.list(hugene21sttranscriptclusterGO), paste, collapse=", "),
                    PATH=sapply(as.list(hugene21sttranscriptclusterPATH), paste, collapse=", "),
                    REFSEQ=sapply(as.list(hugene21sttranscriptclusterREFSEQ), paste, collapse=", "))
head(Annot)

my_frame <- data.frame(exprs(normData))#### so this is the expression data 
head(my_frame)
save(my_frame,file=paste("my_frame.RData",sep=""))


featureData(normData) <- getNetAffx(normData, "transcript")
head(Annot) # looks like all NA but actually it just happens that the first chunk of data are NA and those furtehr down are annotated!
write.csv(Annot,file="anno 100523.txt")
write.csv(my_frame,file="exprs 100523.txt")
all <- merge(Annot, my_frame, by.x=0, by.y=0, all=T)#53617
write.csv(all,file="log2 transformed and normalized data plus some anno 100523.txt")
save(all,file=paste("all.RData",sep=""))



# Analysis ----------------------------------------------------------------


##DNA methylation

####Analysis of HASM extended
###Start with simple asthma vs non-asthma
####Demographics
### Age
# asthma Mean	42.13 SEM	3.32
# non-asthma Mean	44.71 SEM	3.80
#Unpaired t test	
#P value	0.6137

###Sex
##Data analyzed	Asthmatic 	Non-asthmatic	Total
#Male	13	13	26
#Female	10	4	14
#Total	23	17	40
#P value and statistical significance	
#Test	Fisher's exact test
# P value	0.3152

#####Smoking
#Data analyzed	Asthmatic 	Non-asthmatic	Total
#Ex	5	4	9
#Never	17	13	30
#Total	22	17	39
##P value and statistical significance	
#Test	Fisher's exact test
#P value	>0.9999

###OK cool no need to adjust for anything

## need to remove reps!
sampleNames(HEpc) ##HE post combat

reps<-c("AZAC12_rep","MMP1_H01_rep", "MMP1_H04_rep", "AZAD01_rep","AZCC26_rep")
HEpcNR<-HEpc[,!colnames(HEpc)%in%reps] #789613 39 
sampleNames(HEpcNR)
save(HEpcNR,file=paste(HE_dir, "HEpcNR(HE post combat and reps removed).RData",sep=""))  
pDatNR<-pData(HEpcNR)
write.csv(pDatNR,file="pDat all samples HASMe no reps.txt")
BetasNR<-getBeta(HEpcNR)
write.csv(BetasNR,file="Betas all samples HASMe no reps.txt")

### limma on M values

ex<-getM(HEpcNR)
head(ex)
head(pData(HEpcNR))
pDat<-pData(HEpcNR)

design<-model.matrix(~Asthmatic, pDat) ###AsthmaticNon
all(row.names(design)==colnames(ex))  

fit <- lmFit(getM(HEpcNR), design, na.action=na.omit)
fit <-eBayes(fit)
topt <- topTable(fit, coef=c("AsthmaticNon"), adjust="BH", number=Inf, p.value=0.05) # all significant hits
dim(topt) ##36
toptt <- topTable(fit, coef=c("AsthmaticNon"), adjust="BH", number=Inf)
head(toptt)

plot(hist(toptt$P.Value))  ## quite a nice distribution
n<-hist(toptt[,"P.Value"], main=NULL, xlab="Raw p values", ylab = "# of CpGs", freq=TRUE)
box()
abline(h=nrow(toptt)/length(n$breaks),lty=2,col="black")


### lets have a look at what they are 
sigcgs<-row.names(topt)
anno<-getAnnotation(HEpcNR)
head(anno)
siganno<- anno[rownames(HEpcNR)%in%sigcgs,]
write.csv(siganno, file = paste0('siganno annotation for the 36 individually sig hit limma p0.05.csv')) 

stripchart(getBeta(HEpcNR)["cg08064403",]~ pData(HEpcNR)[,14], vertical=TRUE, ylab = "cg08064403 Beta Value", xlab= NULL,) 

## volcano...

# delta beta of these 36
asthmabetas <- getBeta(HEpcNR)[,pData(HEpcNR)$"Asthmatic"%in%"Asthma"]
head(asthmabetas)
colnames(asthmabetas)
Nonasbetas <- getBeta(HEpcNR)[,pData(HEpcNR)$"Asthmatic"%in%"Non"]
colnames(Nonasbetas)

Asthmamean<- as.matrix(rowMeans(asthmabetas))
NAmean<- as.matrix(rowMeans(Nonasbetas))
deltabeta<-NAmean-Asthmamean ##Non asthma as reference
write.table(deltabeta, file="delta beta nonas vs asthma all probes all samples.csv")
head(deltabeta)
##limit to just the 36
sigcgs<-rownames(topt)
sigcgs ###36

head(rownames(deltabeta))
sigdb<-deltabeta[rownames(topt),]
head(sigdb)
sigdbm<-as.matrix(sigdb)

sigdelta <- as.matrix(sigdbm[abs(sigdbm)>0.1,])
dim(sigdelta) ###32 >0.1
sigdelta2 <- as.matrix(sigdbm[abs(sigdbm)>0.2,])
dim(sigdelta2) ###11 >0.2
sigdelta3 <- as.matrix(sigdbm[abs(sigdbm)>0.3,])
dim(sigdelta3) ###3>0.3

head(sigdbm)
plot(hist(sigdbm), xlab = "delta beta", ylim =c(0,20), main=NULL)

write.table(sigdb, file="delta beta for the 36 sig CpGs in all samples.csv")

all(rownames(toptt)==rownames(deltabeta)) ### False
head(toptt)
head(deltabeta) # need to order toptt
toptt2<- toptt[order(rownames(toptt)),]
head(toptt2)
deltabeta2<-as.matrix(deltabeta[order(rownames(deltabeta)),])

all(rownames(toptt2)==rownames(deltabeta2)) ##TRUE :)


head(toptt2)
head(deltabeta2)
forVP<-cbind(toptt2, deltabeta2[,1]) ### new column deltabeta2[, 1]
head(forVP)
colnames(forVP)<-c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "DB")
head(forVP)

forVPdf<- as.data.frame(forVP)


with(forVPdf, plot(DB, -log10(P.Value),pch=19, cex=0.5, xlab="Methylation difference (delta beta)" , ylab= "-log10(pvalue)" ),cex.axis=1.5, cex.lab=1.5)
with(subset(forVPdf, adj.P.Val<.05 ), points(DB, -log10(P.Value), pch=19, cex=0.5, col="darkgray"))
with(subset(forVPdf, adj.P.Val<.05 & DB>=0), points(DB, -log10(P.Value), pch=19, cex=0.5, col="blue"))
with(subset(forVPdf, adj.P.Val<.05 & DB<=0), points(DB, -log10(P.Value), pch=19, cex=0.5, col="red"))
toptt[36:37,] ## to get p value line equivalent
abline(h=5.617, col="gray", lty="longdash")
text(-0.4, 5.617,adj=c(0,0), labels = "2.411791e-06", cex=1.2)
axis(1, at=c(-0.1, -0.3,-0.5, 0.1, 0.3))
grid(nx=8, ny=4, col = "lightgray", lty="dotted")


#Gene expression

####pull out just HASMs 
load("D:/OneDrive - The University of Nottingham/Bioinf MRes/Project/R Objects etc/normData.RData")
load("D:/OneDrive - The University of Nottingham/Bioinf MRes/Project/R Objects etc/all.RData")
load("D:/OneDrive - The University of Nottingham/Bioinf MRes/Project/R Objects etc/my_frame.RData")

colnames(all)
HASM<-c("SMACD04","MMP1-A02","MMP-1A08","MMP-1A12","AZAD03","AZCC14","AZCC15","AZCC23","MMP-1A11","MMP1-A14","MMP-1A16","AZCC01","ALS9007","ALS9008","AZAD01","MMP1-A06","MMP1-A03","AZAC05","AZAC07","AZAC09","AZAC11","AZAC03","AZAC12","AZCC16","AZCC21","AZCC22","AZCC25","AZCC26","MMP-1H10","MMP1-H02","MMP1-H01","HASM0210","HASM0909","2803005","MMP1-H04")

HASMExp<-my_frame[,colnames(my_frame)%in%HASM] ##35 good thats all of them!
save(HASMExp,file=paste("HASMExp.RData",sep=""))

HASMnormdata<-normData[, colnames(normData)%in%HASM]
save(HASMnormdata,file=paste("HASMnormdata.RData",sep=""))
sampleNames(HASMnormdata) ##35 good

# [1] "2803005"  "ALS9007"  "ALS9008"  "AZAC03"   "AZAC05"   "AZAC07"   "AZAC09"   "AZAC11"   "AZAC12"   "AZAD01"   "AZAD03"   "AZCC01"   "AZCC14"   "AZCC15"  
#[15] "AZCC16"   "AZCC21"   "AZCC22"   "AZCC23"   "AZCC25"   "AZCC26"   "HASM0210" "HASM0909" "MMP-1A08" "MMP-1A11" "MMP-1A12" "MMP-1A16" "MMP-1H10" "MMP1-A02"
#[29] "MMP1-A03" "MMP1-A06" "MMP1-A14" "MMP1-H01" "MMP1-H02" "MMP1-H04" "SMACD04" 

design2<-cbind(NonAs=c(0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,1,0,0,0,0,1,1,1,0), Asth=c(1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,1,0,0,1,1,1,1,1,1,0,1,1,1,1,0,0,0,1))
design2
fit<-lmFit(HASMnormdata, design2)

cont.matrix<-makeContrasts(nonAsvsAs=NonAs-Asth, levels=design2)
fit2<-contrasts.fit(fit, cont.matrix)
fit2<-eBayes(fit2)
tt<-topTable(fit2, adjust= "BH", number = Inf, p.value=0.05)
dim(tt) ###0
ttall<- topTable(fit2, adjust= "BH", number = Inf)
write.csv(ttall, file = paste0(HE_dir, 'ttall toptable for all genes.csv')) 
plot(hist(ttall$P.Value))  ## not a bad distribution - might just be power.
n<-hist(ttall[,"P.Value"], main=NULL, xlab="Raw p values", ylab = "# of Genes", freq=TRUE)
box()
abline(h=nrow(ttall)/length(n$breaks),lty=2,col="black")


#### pull out genes assd with the 36 diff methylation CpGs (33 genes)

head(HASMnormdata)
fD<-featureData(HASMnormdata)
featureNames(HASMnormdata)
head(exprs(HASMnormdata))
head(all)

####So should be able to pull probes associated with genes (using SYMBOL) from "all" and then limit the ExpressionSet to these probes for Limma
sigG<-c("ACP5","BFAR","C11orf21","C20orf166","C20orf85","CAMTA1","CLDND1","CPE","EEFSEC","EPS8L1","FMNL2","FOPNL","HHIPL2","KCNS1","LOC613266","MSLNL","PAQR4","RHOQ","RPL13AP17","RTN4RL1","SNX29","SPATA13","THBS1","TNRC6C","TPM1","TSPAN11","WFIKKN1","WNT4","TSPAN32","C20orf200","MACROD2","MIR662","MAGI2")
allsigG<-all[all$SYMBOL%in%sigG,]
allsigG
allsigG$SYMBOL #26 present. RHOQ found 3 times  ##### CAMTA1    WNT4      HHIPL2    TSPAN32   C11orf21  TSPAN11   SPATA13   SPATA13   THBS1     TPM1      WFIKKN1   MIR662    PAQR4     SNX29     BFAR      FOPNL     TNRC6C   
#[18] RTN4RL1   EPS8L1    ACP5      RHOQ      FMNL2     MACROD2   C20orf85  KCNS1     EEFSEC    CLDND1    RPL13AP17 MAGI2     RHOQ      RHOQ     
head(rownames(all))
RNAprobes<-allsigG$Row.names

HASMnDmethsig<- HASMnormdata[rownames(HASMnormdata)%in%RNAprobes] ##31 features, 35 samples 

ExHASMnDmethsig<-exprs(HASMnDmethsig)
head(ExHASMnDmethsig)
save(HASMnDmethsig,file=paste("HASMnDmethsig.RData",sep=""))
write.csv(ExHASMnDmethsig, file = paste0('ExHASMnDmethsig gene expression annotated to diff methyd CpGs.csv')) 
##these were analysed for diferential expression in PRISM



# Restricted sample analysis ----------------------------------------------

load("D:/OneDrive - The University of Nottingham/Bioinf MRes/Project/R Objects etc/HEpcAM(HE post combat, just AZA and MMP1 samples).RData")
##HEpcAM - HASM extended, post combat, AZA and MMP1 samples only
HEpcAM ##dim: 789613 24
sampleNames(HEpcAM)
##  [1] "AZAD01"        "AZAD03"        "AZAD08"        "MMP_1A02_AZ02" "MMP_1A08"      "MMP_1A11"      "MMP_1A12"     
##[8] "MMP_1A14"      "MMP_1A16"      "MMP1_A03"      "MMP1_A06"      "MMP1_A07"      "SMACD04"       "AZAC03"       
##[15] "AZAC05"        "AZAC07"        "AZAC09"        "AZAC11"        "AZAC12"        "MMP_1H10"      "MMP1_H01"     
##[22] "MMP1_H02"      "MMP1_H04"      "MMP1_H05" 

##Demographics
##Table Analyzed	Age

#Column B	Asthmatic
#vs.	vs.
#Column A	Non-Asthmatic

#Unpaired t test	
#P value	0.9943

##Table Analyzed	Sex

#P value and statistical significance	
#Test	Fisher's exact test
#  P value	0.6792

##Data analyzed	Non-Asthmatic	Asthmatic	Total
#Male	8	8	16
#Female	3	5	8
#Total	11	13	24

##Table Analyzed	Smoking		

#P value and statistical significance			
#Test	Fisher's exact test		
# P value	>0.9999		
#P value summary	ns		
#One- or two-sided	Two-sided		
#Statistically significant (P < 0.05)?	No		

#Data analyzed	Non-Asthmatic	Asthmatic	Total
# Never	11	12	23
#Ex	0	1	1
#Total	11	13	24

### No need to adjust

ex<-getM(HEpcAM)
head(ex)

pDatAM<-pData(HEpcAM)
head(pDatAM)


design<-model.matrix(~Asthmatic, pData(HEpcAM))

all(row.names(design)==colnames(ex))  

require("limma")
fit <- lmFit(getM(HEpcAM), design, na.action=na.omit)
fit <-eBayes(fit)
topt <- topTable(fit, coef=c("AsthmaticNon"), adjust="BH", number=Inf, p.value=0.05) # all significant hits
dim(topt) ##BH p<0.05 356, Bonf p<0.05 = 2
head(topt)
toptt <- topTable(fit, coef=c("AsthmaticNon"), adjust="BH", number=Inf)
head(toptt)
toptbonf <- topTable(fit, coef=c("AsthmaticNon"), adjust="bonferroni", number=Inf, p.value=0.05)
toptbonf


plot(hist(toptt$P.Value))  
n<-hist(toptt[,"P.Value"], main=NULL, xlab="Raw p values", ylab = "# of CpGs", freq=TRUE)
box()
abline(h=nrow(toptt)/length(n$breaks),lty=2,col="black")

stripchart(getBeta(HEpcAM)["cg05804504",]~ pData(HEpcAM)[,14], vertical=TRUE, ylab = "cg05804504 Beta Value", xlab= NULL,) 
stripchart(getBeta(HEpcAM)["cg21969101",]~ pData(HEpcAM)[,14], vertical=TRUE, ylab = "cg21969101 Beta Value", xlab= NULL,) 


write.csv(toptt, file = paste0('toptt all probes AM only.csv')) 

### lets have a look at what they are 
sigcgs<-row.names(topt)
anno<-getAnnotation(HEpcAM)
head(anno)
siganno<- anno[rownames(HEpcAM)%in%sigcgs,]
write.csv(siganno, file = paste0('Anno 356 limma p0.05 AM only.csv')) 


## Delta betas and volcano plot of hits
NAbetas<- getBeta(HEpcAM)[,pData(HEpcAM)$"All_Diag"%in%"Healthy"]
head(NAbetas)

Asbetas <- getBeta(HEpcAM)[,pData(HEpcAM)$"All_Diag"%in%"Asthma"]
head(Asbetas)

NAmean<- as.matrix(rowMeans(NAbetas))
Asmean<- as.matrix(rowMeans(Asbetas))
deltabeta<-NAmean-Asmean ##NA as reference
write.table(deltabeta, file="delta beta NA minus asthmatic all probes AM samples.csv")
head(deltabeta)
##limit to just the 356
sigcgs<-rownames(topt)
sigcgs ###356

head(rownames(deltabeta))
sigdb<-deltabeta[rownames(topt),] #356 :)
head(sigdb)
sigdbm<-as.matrix(sigdb)

sigdelta <- as.matrix(sigdbm[abs(sigdbm)>0.1,])
dim(sigdelta) ###287 >0.1
sigdelta2 <- as.matrix(sigdbm[abs(sigdbm)>0.2,])
dim(sigdelta2) ###127 >0.2
sigdelta3 <- as.matrix(sigdbm[abs(sigdbm)>0.3,])
dim(sigdelta3) ###31 >0.3
sigdelta4 <- as.matrix(sigdbm[abs(sigdbm)>0.4,])
dim(sigdelta4) ###1 >0.4 ##cg22854855 0.4176961

write.table(sigdb, file="delta beta for the 356 sig CpGs AM samples.csv")

plot(hist(sigdbm), xlab = "delta beta", ylim =c(0,150), xlim = c(-0.6,0.6), main=NULL)


all(rownames(toptt)==rownames(deltabeta)) ### False
head(toptt)
head(deltabeta) # need to order toptt
toptt2<- toptt[order(rownames(toptt)),]
head(toptt2)
deltabeta2<-deltabeta[order(rownames(deltabeta)),]
all(rownames(toptt2)==rownames(deltabeta2)) ##TRUE :)
head(deltabeta2)

deltabeta2m<-as.matrix(deltabeta2)
head(deltabeta2m)

forVP<-cbind(toptt2, deltabeta2m[,1]) ### new column sigdbm[, 1]
head(forVP)
colnames(forVP)<-c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "DB")
head(forVP)

forVPdf<- as.data.frame(forVP)


with(forVPdf, plot(DB, -log10(P.Value),pch=19, cex=0.5, xlab="Methylation difference (delta beta)" , ylab= "-log10(pvalue)" ),cex.axis=1.5, cex.lab=1.5)
with(subset(forVPdf, adj.P.Val<.05 ), points(DB, -log10(P.Value), pch=19, cex=0.5, col="darkgray"))
with(subset(forVPdf, adj.P.Val<.05 & DB>=0), points(DB, -log10(P.Value), pch=19, cex=0.5, col="blue"))
with(subset(forVPdf, adj.P.Val<.05 & DB<=0), points(DB, -log10(P.Value), pch=19, cex=0.5, col="red"))
toptt[356:357,] ## to get p value line equivalent
abline(h=4.64, col="gray", lty="longdash")
text(-0.4, 4.6,adj=c(0,0), labels = "p=2.271216e-05", cex=1.2)


##### expression of associated genes
##load my_frame
load("D:/OneDrive - The University of Nottingham/Bioinf MRes/Project/R Objects etc/my_frame.RData")
colnames(my_frame) ##all samps
## need just HASMs MMP and AZA
AMEx<-c("SMACD04","AZAC03", "AZAC05", "AZAC07","AZAC09","AZAC11","AZAC12","AZAD01","AZAD03","MMP.1A08","MMP.1A11","MMP.1A12","MMP.1A16","MMP.1H10","MMP1.A02","MMP1.A03","MMP1.A06","MMP1.A14","MMP1.H01","MMP1.H02","MMP1.H04")
## 21 samples - we are missing some compared to DNA methylation as the RNA quality was not suitable
sampleNames(HEpcAM)
#AZAD08, MMP1A07, MMP1-HO5

AMExp<-my_frame[,colnames(my_frame)%in%AMEx] ##21, good
load("D:/OneDrive - The University of Nottingham/Bioinf MRes/Project/R Objects etc/normData.RData")
sampleNames(normData) #D:\OneDrive - The University of Nottingham\UoN Box Migration\Asthma UK_MRF\post funding\new cells added nov 2016\RNA\Affy\Affy Gene expression All cells Healthy plus HASM asthmatic
AMEx2<- c("AZAC03","AZAC05","AZAC07","AZAC09","AZAC11","AZAC12","AZAD01","AZAD03","MMP-1A08","MMP-1A11","MMP-1A12","MMP-1A16","MMP-1H10","MMP1-A02","MMP1-A03","MMP1-A06","MMP1-A14","MMP1-H01","MMP1-H02","MMP1-H04","SMACD04")
AMnormdata<-normData[, colnames(normData)%in%AMEx2]
sampleNames(AMnormdata) ##21, super

save(AMnormdata, file=paste0('AMnormdata.rda'))

designEx<-cbind(Asthma=c(0,0,0,0,0,0,1,1,1,1,1,1,0,1,1,1,1,0,0,0,1), Healthy=c(1,1,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,1,1,1,0))
designEx

fitEx<-lmFit(AMnormdata, designEx)
cont.matrix<-makeContrasts(NAvsA=Healthy-Asthma, levels=designEx)
fitEx<-contrasts.fit(fitEx, cont.matrix)
fitEx2<-eBayes(fitEx)
ttEx<-topTable(fitEx2, adjust= "BH", number = Inf, p.value=0.05)
dim(ttEx) ##0
ttExall<- topTable(fitEx2, adjust= "BH", number = Inf)

plot(hist(ttExall$P.Value))  
n<-hist(ttExall[,"P.Value"], main=NULL, xlab="Raw p values", ylab = "# of probes", freq=TRUE)
box()
abline(h=nrow(ttExall)/length(n$breaks),lty=2,col="black")

head(ttExall)
write.csv(ttExall, file = paste0('ttExall all HASM extended genes AM only.csv')) 
getwd()

###nominal p<0.001, 158 probes

Expressiom_probes_nominal_p_0_001_AM_samples_ <- read_csv("Expressiom probes nominal p 0.001 AM samples .csv")
Parsed with column specification:
  cols(
    probe = col_double()
  )

pnames<-Expressiom_probes_nominal_p_0_001_AM_samples_
head(pnames)

pnamesm<-as.matrix(pnames)
head(pnamesm)
rownames(pnamesm)<-pnamesm[,1]
head(pnamesm)
load("D:/OneDrive - The University of Nottingham/Bioinf MRes/Project/R Objects etc/all.RData")
head(all)
head(row.names(all))
all2<-all
rownames(all2)<-all2[,1]
head(all2)


sigEX<- all2[rownames(all2)%in%pnamesm,] ### 
head(sigEX)
write.csv(sigEX, file = paste0('Expression and fData for the 158 nominally diff expressed genes Am.csv')) 

#### Look at genes associated with the 356 diff methylation CpGs (263 genes)

head(AMnormdata)
head(all)

####So should be able to pull probes associated with genes (using SYMBOL) from "all" and then limit the ExpressionSet to these probes for Limma
Gene_names_Sig_CpGs_AM <- read_csv("Gene names Sig CpGs AM.csv")
Parsed with column specification:
  cols(
    Gene = col_character()
  )

sigG<-Gene_names_Sig_CpGs_AM
sigGm<-as.matrix(sigG)
rownames(sigGm)<-sigGm[,1]
head(sigGm)

sigG<-c("ACP5","BFAR","C11orf21","C20orf166","C20orf85","CAMTA1","CLDND1","CPE","EEFSEC","EPS8L1","FMNL2","FOPNL","HHIPL2","KCNS1","LOC613266","MSLNL","PAQR4","RHOQ","RPL13AP17","RTN4RL1","SNX29","SPATA13","THBS1","TNRC6C","TPM1","TSPAN11","WFIKKN1","WNT4","TSPAN32","C20orf200","MACROD2","MIR662","MAGI2")
allsigG<-all[all$SYMBOL%in%rownames(sigGm),]
allsigG ##231 observations #### missing 31 genes
save(allsigG,file=paste("allsigG gene expression object for genes assd with the 356 DMPs.RData",sep=""))
allsigG$SYMBOL #231 present.      
head(rownames(all))
RNAprobes<-allsigG$Row.names

HASMnDmethsig<- AMnormdata[rownames(AMnormdata)%in%RNAprobes] ##231 features, 21 samples 

### Pull the DNAm assd probes from TT of full expression set and see p vals
##RNAprobes from ttExall
ttExDm<- ttExall[rownames(ttExall)%in%RNAprobes,] ###231 probes
head(ttExDm)
write.csv(ttExDm, file = paste0('Expression TT subset by DNAm assd 231p.csv'))

#22 probes have a nominal p<0.05
pnames3<-c("16699932","16812982","16769598","16854540","16720852","16849238","17047744", "16822296","16947337","16750154","16720165","16975803","16828696","16911804","16849098","16879721","16988423","16975671","17073495","17054852","17000180","17014257")
sigEX3<- all2[rownames(all2)%in%pnames3,] ### 
head(sigEX3)
write.csv(sigEX3, file = paste0('Expn fData 22 nom DNAm assd diff expd full TT .csv')) 

# eTQM with lm ------------------------------------------------------------
## create function
limma.tt.exprs.nocor <- function(eDat, pDat, var.interest, form, fDat, Subject){
  library(limma)
  form <- as.formula(form)
  des <- model.matrix(form, data=pDat)
  
  fit <- lmFit(eDat, des,na.action=na.omit)
  eb.fit <- eBayes(fit)
  
  tt <- topTable(eb.fit, coef = var.interest, number = nrow(eDat), adjust.method = 'BH')
  cat(paste0('\nreturning coef of ', var.interest,' \n'))
  
  tt$Gene <- fDat$gene_symbol[match(rownames(tt), fDat$probe_id)]
  tt$Probe <- rownames(tt)
  tt %>% dplyr::select(Gene, Probe, logFC, AveExpr, t, P.Value, adj.P.Val)
  tt
}


#####Prepare your data and make the following objects 
#eDat = log2-transformed gene expression matrix of your differentially expressed genes

eDat<-GE_final
head(eDat)

#bDat = M values matrix of your differentially methylated sites
bDat<-Meth_sig #356
head(bDat)
##hang on this needs to be significant ones only


# pDat = phenotype table ##needs sample names to match...... e and bDat :s....
pDat<-pDatAM
head(pDat)
rownames(pDat)

### need to exlude samples that are not in the matched datasets first
nonmatch<-c("AZAD08", "MMP1_A07","SMACD04", "MMP1_H05")

pDatmatch<- pDat[!rownames(pDat)%in%nonmatch,] 
rownames(pDatmatch)
colnames(eDat)
# this is the order the rows need to be in with row naes from the pDat object "AZAC03","AZAC05","AZAC07", "AZAC09","AZAC11","AZAC12", "AZAD01", "AZAD03", "MMP_1A02_AZ02""MMP1_A03","MMP1_A06", "MMP_1A08", "MMP_1A11", "MMP_1A12","MMP_1A14","MMP_1A16", "MMP1_H01","MMP1_H02","MMP1_H04","MMP_1H10"         
test<-pDatmatch[c("AZAC03","AZAC05","AZAC07", "AZAC09","AZAC11","AZAC12", "AZAD01", "AZAD03", "MMP_1A02_AZ02","MMP1_A03","MMP1_A06", "MMP_1A08", "MMP_1A11", "MMP_1A12","MMP_1A14","MMP_1A16", "MMP1_H01","MMP1_H02","MMP1_H04","MMP_1H10"),]

rownames(test) ###cool this is in the right order!
head(test) # and it's taken teh other columns with it
rownames(test)==colnames(eDat) ## OK, happy these are in correct order just format of sample name a bit different. So if I do rownames(test)=colnames(eDat) does that work?
test2<-test
rownames(test2)=colnames(eDat)
rownames(test2)
colnames(eDat)
###sorted finally! Save pDat!

pDat_final<-test2
save(pDat_final, file=paste0('pDat_final.rda'))
pDat<-pDat_final


# var.interest = the name of the variable of interest, append the non-reference level to the name if it is a factor, eg. SexF
# form = formula for limma's design matrix

# fDat = gene expression annotation. Need to have two columns: 'probe_id' and 'gene_symbol'
fDat<-all_anno[,c(1,3)]
head(fDat)
names(fDat)<-c('probe_id', 'gene_symbol')
#Subject=Vector of donor ID for your samples

# assign each CpG site in bDat to a block - results of each block will be saved as 1 file. Assign size <50,000 for each block or it will take ages to run. If your hits are <50000 then you can add them all to one block. 
require("data.table")
size=356
nBlock=nrow(bDat)/size+1
block=1
for (block in 1:nBlock){
  # starting and ending index of each block
  start=(block-1)*size+1
  end=min(block*size, nrow(bDat))
  # do 'expression ~ methylation' for all the CpGs in a block
  index=1
  result_hits=list()
  for (i in start:end){
    pDat = pDat %>% mutate(BetaVal=as.numeric(bDat[i,])) %>% as.data.table() 
    cpgsite = rownames(bDat)[i]
    form_hits = '~BetaVal'
    result_hits[[index]] = limma.tt.exprs.nocor(eDat,pDat,'BetaVal',form_hits,fDat,Subject1) %>% mutate(cpg = cpgsite)
    index = index+1
    rm(pDat_1);rm(cpgsite)
  }
  # save the results after completing a block - have applied an adj.pval cutoff of 0.05
  result_hits = rbindlist(result_hits)
  result_hits<-result_hits[which(result_hits$adj.P.Val<0.05),]
  fwrite(result_hits, paste0('toptable_hits0.05',block,'.csv'))
  print(paste0('Block ',block, ' completed!'))
  # memory release
  rm(result_hits);gc();gc()
}

### saves a csv that goes into the WD contains 279 gene:CpG pairs. 145 unique genes, 24CpGs. Flot1 appears the most frequently so going to have a look at it's epxression

write.table(all_anno, file = "all_anno annotation and gene expression data.txt", sep = "\t",
            row.names = TRUE, col.names = NA)



###determine if genes correlated by eQTM are differentially expressed..
eTQMgenes<-as.list(eQTM_genes_201222)
head(eDat)
eQTMprobes<-eQTM_probes_201222$probe
eQTMprobes
ExpneQTMgenes<-eDat[rownames(eDat)%in%eQTMprobes,]
head(ExpneQTMgenes)
write.table(ExpneQTMgenes, file = "Expression values of eQTM genes.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

###loook at expression of these by lm
ExpneQTMgenes<-Expression_values_of_eQTM_probes
ExpneQTMgenes2 <- data.frame(ExpneQTMgenes, row.names = 1)
head(ExpneQTMgenes2)
colnames(ExpneQTMgenes2)==rownames(pdata) ###now all true, super!

pdata<-pDat
head(pdata)
head(ExpneQTMgenes)

lmHE<-sapply(1:nrow(ExpneQTMgenes2), function(x){z<-lm(unlist(ExpneQTMgenes2[x,]) ~ pdata$Asthmatic)
coef(summary(z))["pdata$AsthmaticNon","Pr(>|t|)"]})


head(lmHE)
plot(hist(lmHE)) ### looks good 


lmHE<-as.matrix(lmHE)
head(lmHE)
head(exHEAMv)
rownames(lmHE) <- rownames(ExpneQTMgenes) ## makes column one the row names
head(lmHE)

paHEd<-p.adjust(lmHE[,1], method="BH") 
head(paHEd)

paHEd<- cbind(lmHE, paHEd)
head(paHEd)
colnames(paHEd) <- c("pval","p.adj") 
head(paHEd)
paHEd2 <- paHEd[order(paHEd[,2]),] 

write.table(paHEd2, file = "lm results from eQTM hit probes.txt", sep = "\t",
            row.names = TRUE, col.names = NA)


head(paHEd2)
dim(paHEd2)
n<-hist(paHEd2[,"pval"], main=NULL, xlab="Raw p values", ylab = "# of Probes", freq=TRUE)
box()
abline(h=nrow(paHEd2)/length(n$breaks),lty=2,col="black")


sum(paHEd[,2]<0.05) # 147

sum(paHEd[,1]<0.05) #153 with nominal p<0.05

#What are the 147 probes and associated gene names
paHEd2<-read.table(file = "lm results from eQTM hit probes.txt", sep = "\t")
paHEd2<-as.data.frame(paHEd2)
head(paHEd2)

eQTMsigprobes<-with(paHEd2, V1[V3<0.05]) ##147 nice! What genes are they. What is my gene annotation called....agggh!
load("D:/OneDrive - The University of Nottingham/UoN Box Migration/Asthma UK_MRF/post funding/new cells added nov 2016/R Ananlysis/AZA MMP1 only/all_anno all annotated genes with probe as row name.rda")
head(all_anno)
colnames(all_anno)
eQTMsigGenes<-all_anno[rownames(all_anno)%in%eQTMsigprobes,]
eQTMsigGenesnames<-(eQTMsigGenes$SYMBOL)

write.table(eQTMsigGenesnames, file = "gene names for 147 probes from eQTM that displayed sigdifexpn.txt", sep = "\t",
            row.names = TRUE, col.names = NA)


###get fold changes of these probes
tt_test<-ttExall_all_HASM_extended_genes_AM_only
head(tt_test)
tt2<-as.data.frame(tt_test)
rownames(tt2)<-tt2$...1
head(tt2) ##good rownames are probe IDS


FCeQTMprobes<-tt2[rownames(tt2)%in%eQTMsigprobes,]
write.table(FCeQTMprobes, file = "TopTable 147 probes from eQTM that displayed sigdifexpn.txt", sep = "\t",
            row.names = TRUE, col.names = NA)




# Replication in GEO dataset ----------------------------------------------


##Pull my sig CpGs from GSE146376 data (Cytokine-induced molecular responses in airway smooth muscle cells inform genome-wide association studies of asthma)
#Read in data from GEO - in expression set format. Extract pDat from this and then use minfi to read in as RGSet
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery", force = TRUE)
require(GEOquery)
gset_EPIC <- getGEO("GSE146376", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset_EPIC) > 1) idx <- grep("GPL21145", attr(gset_EPIC, "names")) else idx <- 1
gset_EPIC <- gset_EPIC[[idx]]
pDat_EPIC<-pData(gset_EPIC)
write.table(pDat_EPIC,file="pDat_EPIC.csv",sep=",")
save(gset_EPIC,file="ExpressionSet_GSE146376_n.rda")

#Downloaded idats into a folder but they are .gz files
#list.files("GSE146376_IDAT/", pattern = "idat")
#idatFiles <- list.files("GSE201872_IDAT/", pattern = "idat.gz$", full = TRUE)
#sapply(idatFiles, gunzip, overwrite = TRUE)

#Reading RGSet for EPIC samples -30 samples
pDat_EPIC<-read.csv("pDat_EPIC.csv",sep=",")
write.table(pDat_EPIC,file="pDat_EPIC.csv",sep=",")
require(minfi)
targets_EPIC <- read.metharray.sheet("GSE201872_IDAT/",pattern = "EPIC.csv$", ignore.case = TRUE,
                                     recursive = TRUE, verbose = TRUE)
RGset_GSE146376_EPIC<-read.metharray.exp(base="GSE201872_IDAT/",targets = targets_EPIC,extended = TRUE,verbose = TRUE)




load("D:/OneDrive - The University of Nottingham/UoN Box Migration/Asthma UK_MRF/post funding/new cells added nov 2016/R Ananlysis/AZA MMP1 only/MSet_Swan_GSE146376.rda")
GSEbetas<-getBeta(MSet_swan)
head(GSEbetas)
GSEMs<-getM(MSet_swan)
head(GSEMs)

load("D:/OneDrive - The University of Nottingham/UoN Box Migration/Asthma UK_MRF/post funding/new cells added nov 2016/R Ananlysis/AZA MMP1 only/pData_GSE146376.rda")

head(pData_GSE146376)
pData_GSE146376$Treatment  ## "Vehicle is the control and thats all i want for now
#Filter pData by Vehicle to get sample names while retaining all other pDat

pDatGSE_Vehicle<-pData_GSE146376[pData_GSE146376$Treatment%in%c("Vehicle"),]
pDatGSE_Vehicle
write.csv(pDatGSE_Vehicle, file = paste0('Pdata GSE vehicle samps.csv'))
vehicle_samps<-pDatGSE_Vehicle$Accession
vehicle_samps
## pull from betas
vehbetas<-GSEbetas[,colnames(GSEbetas)%in%vehicle_samps]
head(vehbetas)
vehMs<-GSEMs[,colnames(GSEMs)%in%vehicle_samps]
head(vehMs)

###now only want sig 356 cgs.....
sigcgs
sigcgs_betas<-vehbetas[rownames(vehbetas)%in%sigcgs,]
sigcgs_betas
write.csv(sigcgs_betas, file = paste0('betas for sig cgs in GSE vehicle samps.csv')) ####291 probes available....
sigcgs_Ms<-vehMs[rownames(vehMs)%in%sigcgs,]
write.csv(sigcgs_Ms, file = paste0('M for sig cgs in GSE vehicle samps.csv'))

colnames(sigcgs_Ms)==rownames(pDatGSE_Vehicle) ##all true, good so can use pDat for model?

design<-model.matrix(~Asthma+Age, pDatGSE_Vehicle)

all(row.names(design)==colnames(sigcgs_Ms))  ##TRUE



lmHE<-sapply(1:nrow(sigcgs_Ms), function(x){z<-lm(unlist(sigcgs_Ms[x,]) ~ pDatGSE_Vehicle$Asthma)
coef(summary(z))["pDatGSE_Vehicle$AsthmaYes","Pr(>|t|)"]})


head(lmHE)
plot(hist(lmHE)) ### looks wonky!


lmHE<-as.matrix(lmHE)
head(lmHE)
head(sigcgs_Ms)
rownames(lmHE) <- rownames(sigcgs_Ms) ## makes column one the row names
head(lmHE)

paHEd<-p.adjust(lmHE[,1], method="BH") 
head(paHEd)

paHEd<- cbind(lmHE, paHEd)
head(paHEd)
colnames(paHEd) <- c("pval","p.adj") 
head(paHEd)
paHEd2 <- paHEd[order(paHEd[,2]),] 

head(paHEd2)
dim(paHEd2)
n<-hist(paHEd2[,"pval"], main=NULL, xlab="Raw p values", ylab = "# of CpGs", freq=TRUE)
box()
abline(h=nrow(paHEd2)/length(n$breaks),lty=2,col="black")


sum(paHEd[,1]<0.05) # nominal p<0.05 =18 of 291
write.csv(paHEd2, file = paste0('lm results from 356 sig hit in GSE146376.csv'))

df<-as.data.frame(paHEd)

#rep_hits<-(paHEd[,1]<0.05)
rep_hits<- df[df$pval < 0.05,]
rep_hitscg<-rownames(rep_hits)

rep_betasGSE<-vehbetas[rownames(vehbetas)%in%rep_hitscg,]
betasHASM_AM<- getBeta(HEpcAM)
re_betasHASM_AM<-betasHASM_AM[rownames(betasHASM_AM)%in%rep_hitscg,]
write.csv(rep_betasGSE, file = paste0('betas of replicated cgs in GSE.csv'))
write.csv(re_betasHASM_AM, file = paste0('betas of replicated cgs in HASMe_AM.csv'))

## replication of eQTM genes with diff expression
##Get the gene expression data
library(GEOquery)
gse146374<-getGEO("GSE146374",GSEMatrix=TRUE)
pdat146374<-pData(phenoData(gse146374[[1]]))
write.csv(pdat, file = paste0('Pdat of gse146374.csv')) 
exprs146374<-exprs(gse146374[[1]])
head(exprs146374)
rownames(exprs146374) ## gene names ....except the ones that are dates!
colnames(exprs146374) ##GSM IDs - did i limit to baseline for methylation and can i use that here? yes it is called vehicle_samps - 

vehexpr<-exprs146374[,colnames(exprs146374)%in%vehicle_samps]
dim(vehexpr)
head(vehexpr)

#### OK doesn't look like the methylation and the expn use the same GSM IDs - so go to expression pData
head(pdat146374)
### get just vehicle - characteristics_ch1.2 treatment: Vehicle
pdat146374_Vehicle<-pdat146374[pdat146374$characteristics_ch1.2%in%c("treatment: Vehicle"),] ##67 samples - so less that in methylation
write.csv(pdat146374_Vehicle, file = paste0('pDat of just the vehicle treated gene expression samples.csv')) ### good.
exprs_veh_samps<-rownames(pdat146374_Vehicle) # super
##limit the expression matrix by these samples
expr_vehicle<-exprs146374[,colnames(exprs146374)%in%exprs_veh_samps]
dim(expr_vehicle) #18279    67
### now need to limit by the genes i'm interested in
eQTM_gene_names_for_val <- read_csv("eQTM_gene_names_for_val.csv", 
                                    +     col_types = cols(genes = col_character()))
genes<-eQTM_gene_names_for_val
genes<-genes$genes
genes

expr_vehicle_eQTMgenes<-expr_vehicle[rownames(expr_vehicle)%in%genes,]
dim(expr_vehicle_eQTMgenes) ##78  67 - only 78 of 148 available (NB 148 includes duplicate gene names - 121 uniques, of which 78 are present in the val data)
val_genes<-rownames(expr_vehicle_eQTMgenes)
write.csv(val_genes, file = paste0('eQTM gene names present in validation data.csv')) 
save(expr_vehicle_eQTMgenes, file=paste0('expr_vehicle_eQTMgenes.rda'))

#### run lm on the genes by asthma vs normal

colnames(expr_vehicle_eQTMgenes)==rownames(pdat146374_Vehicle) ##all true, good so can use pDat for model?

head(pdat146374_Vehicle) #characteristics_ch1.5

design<-model.matrix(~characteristics_ch1.5, pdat146374_Vehicle)

all(row.names(design)==colnames(expr_vehicle_eQTMgenes))  ##TRUE



lmHE<-sapply(1:nrow(expr_vehicle_eQTMgenes), function(x){z<-lm(unlist(expr_vehicle_eQTMgenes[x,]) ~ pdat146374_Vehicle$characteristics_ch1.5)
coef(summary(z))["pdat146374_Vehicle$characteristics_ch1.5asthma: Yes","Pr(>|t|)"]})

head(lmHE)
plot(hist(lmHE)) ### looks good 

lmHE<-as.matrix(lmHE)
head(lmHE)
head(expr_vehicle_eQTMgenes)
rownames(lmHE) <- rownames(expr_vehicle_eQTMgenes) ## makes column one the row names
head(lmHE)

paHEd<-p.adjust(lmHE[,1], method="BH") 
head(paHEd)

paHEd<- cbind(lmHE, paHEd)
head(paHEd)
colnames(paHEd) <- c("pval","p.adj") 
head(paHEd)
paHEd2 <- paHEd[order(paHEd[,2]),] 

head(paHEd2)
dim(paHEd2)
n<-hist(paHEd2[,"pval"], main=NULL, xlab="Raw p values", ylab = "# of CpGs", freq=TRUE)
box()
abline(h=nrow(paHEd2)/length(n$breaks),lty=2,col="black")


sum(paHEd[,1]<0.05) # nominal p<0.05 =6 of 78
write.csv(paHEd2, file = paste0('lm results from 78 sig hit genes in 146374.csv'))

df<-as.data.frame(paHEd)

GErep_hits<- df[df$pval < 0.05,]
GErep_hitsgene<-rownames(GErep_hits)

rep_expnGSE<-expr_vehicle_eQTMgenes[rownames(expr_vehicle_eQTMgenes)%in%GErep_hitsgene,]
write.csv(rep_expnGSE, file = paste0('expn of replicated genes in GSE.csv'))

