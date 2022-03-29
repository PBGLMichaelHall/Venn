# Venn
```r
#Install all necessary library packages

install.packages("devtools")
library(devtools)

install.packages("BiocManager")
library(BiocManager)

install_github("js229/Vennerable",force = TRUE)
library(Vennerable)

install.packages("remotes")
remotes::install_github("jjellis/GenomicVis",force = TRUE)
library(GenomicVis)

BiocManager::install("VariantAnnotation",force = TRUE)
library(VariantAnnotation)

install.packages("ggplot2")
library(ggplot2)

install.packages("stringi")
library(stringi)

install.packages("chromoMap")
library(chromoMap)

install.packages("QTLseqr")
library(QTLseqr)

install.packages("vcfR")
library(vcfR)

#Set Working Directory
setwd("/home/michael/Desktop/GenomicVis")
getwd()

#Confirm R Session has all libraries and packages available
sessionInfo()
```



# Edit readVcf function on line 56 to include indexing from the GenomicVis otherwise you will generate an error in evaluating the argument 'i' in selecting a method for function could not find function elementlengths

# I have modified the function to correct for error
```r

z<-function (f, genome, exclude.filtered = FALSE, read.info = FALSE, 
             read.geno = FALSE, filter.func = NULL) 
{
  params <- if (read.info == FALSE & read.geno == FALSE) 
    ScanVcfParam(info = NA, geno = NA)
  else if (read.info == FALSE & read.geno == TRUE) 
    ScanVcfParam(info = NA)
  else if (read.info == TRUE & read.geno == FALSE) 
    ScanVcfParam(geno = NA)
  else ScanVcfParam()
  vcf <- readVcf(f, genome, params)
  if (is.function(filter.func)) 
    vcf <- filter.func(vcf)
  if (exclude.filtered) 
    vcf <- vcf[fixed(vcf)$FILTER %in% c("PASS", ".")]
  vcf <- vcf[width(ref(vcf)) == 1]
  vcf <- vcf[elementNROWS(alt(vcf)) == 1]
  vcf <- vcf[width(unlist(alt(vcf))) == 1]
  vcf
}
```



# Call in your VCF File
# The original file was freebayes~bwa~IRGSP-1.0~all-mutants-minus-S14~QUAL1000-S15-HOMREF.vcf.gz"

```r
f1 <- 'freebayes~bwa~IRGSP-1.0~all-mutants-minus-S14~QUAL1000-S15-HOMREF.vcf'
```
# We want to plot the quality scores per SNP for each Chromosome

```r
vcf <- read.vcfR(file = f1)
```

# VCF attributes

![Screenshot from 2022-03-18 15-25-37](https://user-images.githubusercontent.com/93121277/159021255-6ee11a89-f906-4620-93dd-2bf26d653fa1.png)


```r
#Convert to tidy data frame

vcf <- vcfR2tidy(vcf)

```
![Screenshot from 2022-03-18 14-15-42](https://user-images.githubusercontent.com/93121277/159009597-2002b451-a108-43ba-a35f-cc109e49d473.png)


# Make a Data frame calling the function below
```r
# Call ChromQual function from QTLseqr package
Data <- QTLseqr::ChromQual(vcf = vcf, SampleName = "S14")

#Refer to DataFrame and find fixed length of Quality vector
l <- length(Data$Qual)

#Make an object equal to the number of SNPs called
SNP<-seq(from = 0,to = l - 1, by = 1)

#Extract quality value per SNP called
Quality<- Data$Qual

#Extract Chromosome Key
Chrom <- Data$CHROM

#Extract Position
POS <- Data$POS

#Extract Depth
DP <- Data$DP

#Extract Reference Allele
REF <- Data$REF

#Make a dataframe with both objects
df <- data.frame(x = POS, y = Quality, z = Chrom, z2 = DP, z3 = REF)

#Call ggplot package to plot the data
d<-ggplot(data = df, aes(x = SNP, y = Quality))
d + geom_point() + labs(title = "Quality per SNP Sample 14") + facet_grid(Chrom) + geom_smooth()

```

# The plot
![Screenshot from 2022-03-18 15-23-10](https://user-images.githubusercontent.com/93121277/159020586-2f46749a-8f94-4234-a7e8-5da77ac4a9f8.png)


# With the function you can specify any sample in the VCF
```r
Data <- QTLseqr::ChromQual(vcf = vcf, SampleName = "S13")

#Call ggplot package to plot the data
d<-ggplot(data = df, aes(x = SNP, y = Quality))
d + geom_point() + labs(title = "Quality per SNP Sample 15") + facet_grid(Chrom) + geom_smooth()
```

# The Plot etc. ...

![Screenshot from 2022-03-18 15-32-06](https://user-images.githubusercontent.com/93121277/159022328-c4032f87-078d-4aca-a355-aa5d776fa15c.png)



# Call the command from the command line to view the VCF header

![Screenshot from 2022-03-17 13-17-09](https://user-images.githubusercontent.com/93121277/158806871-b714b995-82b3-49c5-8fc1-9d2733e597b7.png)

# See how I am comparing samples 15 and 14

![Screenshot from 2022-03-17 13-26-16](https://user-images.githubusercontent.com/93121277/158808252-a43390f9-11d7-402b-8703-8c0266015fe4.png)

```r
#Call the subsettable vcf
f1 <- 'freebayes~bwa~IRGSP-1.0~S14~HOM-VAR.vcf.gz'
```

# The reference genome is on line 5 of the header followed by chromosomes/contigs and their corresponding lengths
![Screenshot from 2022-03-17 13-15-26](https://user-images.githubusercontent.com/93121277/158806497-0fb347a8-e7da-4736-a130-f0b575a4486f.png)



# Your reference Genome

```r
genome <- 'GCF_001433935.1_IRGSP-1.0_genomic.fna'

#Read Your VCF File

vcf_z2<-z(f1, genome ,read.info = TRUE, read.geno = TRUE, exclude.filtered = FALSE, filter.func = FALSE)


#Refer to DFrame object vcf_z2 in global environment and extract Quality data and number of SNPs called
l <- length(vcf_z2@fixed$QUAL)

#Make an object equal to the number of SNPs called
SNP<-seq(from = 0,to = l - 1, by = 1)

#Extract quality value per SNP called
Quality<-vcf_z2@fixed$QUAL

#Unlist Alternate Allele
AlternateAllele<- unlist(vcf_z2@fixed$ALT)

#Convert to character data vector
AlternateAllele<- as.character(AlternateAllele)

#Make a dataframe with both objects
df <- data.frame(x = SNP, y = Quality)

#Call ggplot package to plot the data
d<-ggplot(data = df, aes(x = SNP, y = Quality))
d + geom_point() + labs(title = "Quality per SNP Sample 14 vs. Sample 15") + geom_line(mapping = aes(x = SNP, y = 40)) + geom_smooth()

```

# The plot of quality scores against SNPs called
![Screenshot from 2022-03-17 13-27-39](https://user-images.githubusercontent.com/93121277/158808390-67a943f2-92aa-4c99-a23d-56035826df39.png)

```r
#Color code each observation by Alternate Allele
e<-ggplot(data = df, aes(x = SNP, y = Quality, color = AlternateAllele))
e + geom_point() + labs(title = "Quality per SNP by Alternate Allele")  + geom_smooth()

```
# The plot of quality scores again
![Screenshot from 2022-03-17 13-28-18](https://user-images.githubusercontent.com/93121277/158808514-b084a0d5-057a-461e-9630-d3a2644cfc38.png)


```r

#Set the working directory
setwd("/home/michael/Desktop/SNPREL")
getwd()

#Decompress gzip file if necessary
#vcf.fn <- gunzip("freebayes~bwa~IRGSP-1.0~all-mutants-minus-S14~QUAL1000-S15-HOMREF.vcf.gz")
#Use your own VCF File
vcf.fn <- "Rename.vcf"
snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")
``` 

# R Console output
![Screenshot from 2022-03-18 09-33-12](https://user-images.githubusercontent.com/93121277/158965095-aacf87e1-8a4a-4a1e-be09-7415f4149133.png)

```r
snpgdsSummary("test.gds")

```
# Summary

![Screenshot from 2022-03-18 09-34-30](https://user-images.githubusercontent.com/93121277/158965202-6ff2de9c-0197-4eeb-ba9a-eaaadee31af6.png)


# Open the GDS file

```r
genofile1 <- snpgdsOpen("test.gds")

#Set the seed for consistent reproducible results
set.seed(1000)

# Pruning the set at different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(genofile1, autosome.only=TRUE,ld.threshold=0.7)

```
# SNP Pruning output
![Screenshot from 2022-03-18 09-35-24](https://user-images.githubusercontent.com/93121277/158965405-54b0d403-27a8-4854-9dc2-595830d42479.png)

```r

names(snpset)
snpset.id <- unlist(snpset)

#Open a new development graphics frame
dev.new()
diss <- snpgdsDiss(genofile1)
```
# Dissimilarity Output
![Screenshot from 2022-03-18 09-37-08](https://user-images.githubusercontent.com/93121277/158965652-bf822e77-c76d-4e49-aa72-f66c49eb15cd.png)

```r
hc <- snpgdsHCluster(diss)

rv <- snpgdsCutTree(hc, label.H = TRUE, label.Z = TRUE)
dev.new()
png(file="Rice Dendrogram.png")
snpgdsDrawTree(rv,type = "dendrogram",outlier.col = "red", main = "SNP Relate Rice Dendrogram",y.label.kinship = TRUE,leaflab = "textlike")
dev.off()
```
# Rice Dendrogram Tree

![Screenshot from 2022-03-18 09-38-47](https://user-images.githubusercontent.com/93121277/158965871-b397f9ec-1709-4ba0-b0a4-168057ba30df.png)


```r

#Edit vcf.venn function from GenomicVis, change names to rownames in line 102  
z1 <- function (vcf.files, genome, sample.names = NULL, ...) 
{
  if (length(vcf.files) > 9) 
    stop("do.venn: too many VCF files (maximum is 9)")
  if (is.null(sample.names)) {
    sample.names <- basename(vcf.files)
  }
  else {
    if (length(vcf.files) != length(sample.names)) 
      stop("do.venn: length(sample.name) must equal length(vcf.files)")
  }
  x <- lapply(vcf.files, function(f) rownames((rowData(readVcf(TabixFile(f), genome, ...)))))
  names(x) <- sample.names
  list(venn = Vennerable::Venn(x), data = x)
}


#Choose two samples preferably from the same species in this case it is Rice
sample.names <- c('S2','S4')

#Again call both VCF files
f1 <- 'freebayes~bwa~IRGSP-1.0~S2~HOM-VAR.vcf.gz'
f2 <- "freebayes~bwa~IRGSP-1.0~S4~HOM-VAR.vcf.gz"

#Save the VCF files into a data vector
vcf.files <- c(f1,f2)

#Call the vcf.venn function and provide appropriate arguments 
v1<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)

#Remove anything that is not numerical
v1 <- gsub("[^0-9.-]","",v1$data$S2)

#Repeat both steps for the second VCF File
v2<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v2 <- gsub("[^0-9.-]","",v2$data$S4)

#Create an list object labeling them appropriately by Sample name variables
z99 <- list("Sample2" = v1, "Sample4" = v2)

#Confirm the data structure is in the right form
str(z99)

#Create a formal class Venn Object
w1<-Venn(Sets=z99)

#Now plot it
plot(w1)
```

# A Venn Diagram comparing samples 2 and 4 and showing which SNPs they share in common
![Screenshot from 2022-03-17 14-42-55](https://user-images.githubusercontent.com/93121277/158821765-7f3cb9f7-bebd-444d-80fe-c8990cc733e9.png)

# Choose three samples preferably from the same species in this case it is Rice

```r
sample.names <- c('S1','S3','S13')

#Again call both VCF files
f1 <- 'freebayes~bwa~IRGSP-1.0~S1~HOM-VAR.vcf.gz'
f2 <- "freebayes~bwa~IRGSP-1.0~S3~HOM-VAR.vcf.gz"
f3 <- "freebayes~bwa~IRGSP-1.0~S13~HOM-VAR.vcf.gz"


#Save the VCF files into a data vector
vcf.files <- c(f1,f2,f3)

#Call the vcf.venn function and provide appropriate arguments 
v1<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)

#Remove anything that is not numerical
v1 <- gsub("[^0-9.-]","",v1$data$S1)

#Repeat both steps for the second VCF File
v2<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v2 <- gsub("[^0-9.-]","",v2$data$S3)

#Repeat both steps for the third VCF File
v3<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v3 <- gsub("[^0-9.-]","",v3$data$S13)

#Create an list object labeling them appropriately by Sample name variables
z99 <- list("Sample1" = v1, "Sample3" = v2, "Sample13" = v3)

#Confirm the data structure is in the right form
str(z99)

#Create a formal class Venn Object
w1<-Venn(Sets=z99)

#Now plot it
plot(w1)
```

# The plot comparing sample1, sample3, and sample13 showing which SNPs they share in common
![Screenshot from 2022-03-17 14-44-48](https://user-images.githubusercontent.com/93121277/158822191-b336f464-9700-43c9-b221-55dc8a754ae6.png)



# Choose four samples preferably from the same species in this case it is Rice

```r
sample.names <- c('S1','S3','S4','S13')

#Again call both VCF files
f1 <- 'freebayes~bwa~IRGSP-1.0~S1~HOM-VAR.vcf.gz'
f2 <- "freebayes~bwa~IRGSP-1.0~S3~HOM-VAR.vcf.gz"
f3 <- "freebayes~bwa~IRGSP-1.0~S4~HOM-VAR.vcf.gz"
f4 <- "freebayes~bwa~IRGSP-1.0~S13~HOM-VAR.vcf.gz"

#Save the VCF files into a data vector
vcf.files <- c(f1,f2,f3,f4)

#Call the vcf.venn function and provide appropriate arguments 
v1<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)

#Remove anything that is not numerical
v1 <- gsub("[^0-9.-]","",v1$data$S1)

#Repeat both steps for the second VCF File
v2<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v2 <- gsub("[^0-9.-]","",v2$data$S3)

#Repeat both steps for the third VCF File
v3<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v3 <- gsub("[^0-9.-]","",v3$data$S4)

#Repeat both steps for fourth VCF file
v4<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v4 <- gsub("[^0-9.-]","",v4$data$S13)

#Create an list object labeling them appropriately by Sample name variables
z99 <- list("Sample1" = v1, "Sample3" = v2, "Sample4" = v3, "Sample13" = v4)

#Confirm the data structure is in the right form
str(z99)

#Create a formal class Venn Object
w1<-Venn(Sets=z99,Weight = FALSE)

#Now plot it
plot(w1)

```

# A Plot showing which SNPs sample1, sample3, sample4, and sample 13 share in common
![Screenshot from 2022-03-17 14-46-38](https://user-images.githubusercontent.com/93121277/158822648-419d366b-2fbb-4e61-83eb-5a5a99cf8f0c.png)



# Group the smaller cluster with samples

```r

f10 <-"freebayes~bwa~IRGSP-1.0~S1~HOM-VAR.vcf.gz"
f11 <-"freebayes~bwa~IRGSP-1.0~S2~HOM-VAR.vcf.gz"
f12 <-"freebayes~bwa~IRGSP-1.0~S3~HOM-VAR.vcf.gz"
f13 <-"freebayes~bwa~IRGSP-1.0~S4~HOM-VAR.vcf.gz"
f14 <-"freebayes~bwa~IRGSP-1.0~S13~HOM-VAR.vcf.gz"

#Create a vcf file vector with each sample
vcf.files <- c(f10,f11,f12,f13,f14)
#Provide appropriate sample names
sample.names <- c('S1','S2','S3','S4','S13')


#Get rid of noise
v6<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v6 <- gsub("[^0-9.-]","",v6$data$S1)
v7<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v7 <- gsub("[^0-9.-]","",v7$data$S2)
v8<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v8 <- gsub("[^0-9.-]","",v8$data$S3)
v9<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v9 <- gsub("[^0-9.-]","",v9$data$S4)
v10<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v10<- gsub("[^0-9.-]","",v10$data$S13)

z99 <- list("Sample1" = v6, "Sample2" = v7, "Sample3" = v8, "Sample4" = v9, "Sample13" = v10)
str(z99)
w1<-Venn(Sets=z99,Weight = FALSE)
plot(Venn(Sets=z99,Weight = F))

```
# A Plot showing shared SNPs for sample1, sample2, sample3, sample4, and sample 13
![Screenshot from 2022-03-17 14-50-29](https://user-images.githubusercontent.com/93121277/158823357-0c8db59d-ce70-4126-b407-1c31b8654082.png)


```r

All_Int <- Reduce(intersect, list(v6,v7,v8,v9,v10))


Pos<-lapply(All_Int, function(f) substr(f, 9:16, 21))


Chr<-lapply(All_Int, function(f) substr(f,6:8,8))


Chr<-lapply(Chr, function(f) gsub("6.1", "Chr1", f))
Chr<-lapply(Chr, function(f) gsub("7.1", "Chr2", f))
Chr<-lapply(Chr, function(f) gsub("8.1", "Chr3", f))
Chr<-lapply(Chr, function(f) gsub("9.1", "Chr4", f))
Chr<-lapply(Chr, function(f) gsub("0.1", "Chr5", f))
Chr<-lapply(Chr, function(f) gsub("1.1", "Chr6", f))
Chr<-lapply(Chr, function(f) gsub("2.1", "Chr7", f))
Chr<-lapply(Chr, function(f) gsub("3.1", "Chr8", f))
Chr<-lapply(Chr, function(f) gsub("4.1", "Chr9", f))
Chr<-lapply(Chr, function(f) gsub("5.1", "Chr10", f))
Chr<-lapply(Chr, function(f) gsub("6.1", "Chr11", f))
Chr<-lapply(Chr, function(f) gsub("7.1", "Chr12", f))

#Generate random alphanumeric code for each SNP
unique_id<-stri_rand_strings(n = 134, length = 4, pattern = "[A-Za-z0-9]")

DF7 <- as.data.frame(cbind(unique_id,Chr,Pos))



df_unlist<-function(df){
  
  df<-as.data.frame(df)
  
  nr<-nrow(df)
  
  c.names<-colnames(df)
  
  lscols<-as.vector(which(apply(df,2,is.list)==TRUE))
  
  if(length(lscols)!=0){
    
    for(i in lscols){
      
      temp<-as.vector(unlist(df[,i]))
      
      if(length(temp)!=nr){
        
        adj<-nr-length(temp)
        
        temp<-c(rep(0,adj),temp)
        
      }
      
      df[,i]<-temp
      
    } #end for
    
    df<-as.data.frame(df)
    
    colnames(df)<-c.names
  }
  return(df)
}


#Prepare an annotations file
DF7<-df_unlist(DF7)
write.table(DF7, file = "SNP.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
```


# A text file with each unique SNP shared in common with the group of samples is written to the current directory.

![Screenshot from 2022-03-18 12-19-44](https://user-images.githubusercontent.com/93121277/158994173-e104ef2e-5dd9-407c-8a7c-3f430ba39f3e.png)



# A Chromosome file with starting positions and lengths is needed for annotation purposes

![Screenshot from 2022-03-18 08-45-29](https://user-images.githubusercontent.com/93121277/158958572-b2fed242-f251-4106-bd91-51860d86b168.png)



```r

Chr1Len <- 43270923
Chr2Len <- 35937250
Chr3Len <- 36413819
Chr4Len <- 35502694
Chr5Len <- 29958434
Chr6Len <- 31248787
Chr7Len <- 29697621
Chr8Len <- 28443022
Chr9Len <- 23012720
Chr10Len <- 23207287
Chr11Len <- 29021106
Chr12Len <- 27531856

ChrStartPos <- c(1,1,1,1,1,1,1,1,1,1,1,1)
ChrEndPos <- c(Chr1Len, Chr2Len, Chr3Len, Chr4Len, Chr5Len, Chr6Len, Chr7Len, Chr8Len, Chr9Len, Chr10Len, Chr11Len, Chr12Len)
length(ChrStartPos)
length(ChrEndPos)
Chr <- c("Chr1", "Chr2", "Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12")
DF8 <- as.data.frame(cbind(Chr,ChrStartPos,ChrEndPos))
DF8 <- df_unlist(DF8)
write.table(DF8, file = "chr.txt", sep = "\t", row.names = FALSE, col.names = FALSE)


chromoMap("chr.txt", "SNP.txt", ploidy = 1)



```

![Screenshot from 2022-03-18 12-31-59](https://user-images.githubusercontent.com/93121277/158995826-7e2f9398-af1f-4584-8a93-9f17b4ba857e.png)




# Group the smaller cluster with samples
```r
f10 <-"freebayes~bwa~IRGSP-1.0~S5~HOM-VAR.vcf.gz"
f11 <-"freebayes~bwa~IRGSP-1.0~S6~HOM-VAR.vcf.gz"
f12 <-"freebayes~bwa~IRGSP-1.0~S7~HOM-VAR.vcf.gz"
f13 <-"freebayes~bwa~IRGSP-1.0~S8~HOM-VAR.vcf.gz"
f14 <-"freebayes~bwa~IRGSP-1.0~S9~HOM-VAR.vcf.gz"
f15 <-"freebayes~bwa~IRGSP-1.0~S10~HOM-VAR.vcf.gz"
f16 <-"freebayes~bwa~IRGSP-1.0~S11~HOM-VAR.vcf.gz"
f17 <-"freebayes~bwa~IRGSP-1.0~S12~HOM-VAR.vcf.gz"


#Create a vcf file vector with each sample
vcf.files <- c(f10,f11,f12,f13,f14,f15,f16,f17)
#Provide appropriate sample names
sample.names <- c('S5','S6','S7','S8','S9','S10','S11','S12')


#Get rid of noise
v11<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v11 <- gsub("[^0-9.-]","",v11$data$S5)
v12<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v12 <- gsub("[^0-9.-]","",v12$data$S6)
v13<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v13 <- gsub("[^0-9.-]","",v13$data$S7)
v14<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v14 <- gsub("[^0-9.-]","",v14$data$S8)
v15<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v15<- gsub("[^0-9.-]","",v15$data$S9)
v16<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v16 <- gsub("[^0-9.-]","",v16$data$S10)
v17<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v17 <- gsub("[^0-9.-]","",v17$data$S11)
v18<-z1(vcf.files, 'GCF_001433935.1_IRGSP-1.0_genomic.fna',sample.names)
v18<- gsub("[^0-9.-]","",v18$data$S12)





All_Int <- Reduce(intersect, list(v11,v12,v13,v14,v15,v16,v17,v18))


Pos<-lapply(All_Int, function(f) substr(f, 9:16, 21))


Chr<-lapply(All_Int, function(f) substr(f,6:8,8))


Chr<-lapply(Chr, function(f) gsub("6.1", "Chr1", f))
Chr<-lapply(Chr, function(f) gsub("7.1", "Chr2", f))
Chr<-lapply(Chr, function(f) gsub("8.1", "Chr3", f))
Chr<-lapply(Chr, function(f) gsub("9.1", "Chr4", f))
Chr<-lapply(Chr, function(f) gsub("0.1", "Chr5", f))
Chr<-lapply(Chr, function(f) gsub("1.1", "Chr6", f))
Chr<-lapply(Chr, function(f) gsub("2.1", "Chr7", f))
Chr<-lapply(Chr, function(f) gsub("3.1", "Chr8", f))
Chr<-lapply(Chr, function(f) gsub("4.1", "Chr9", f))
Chr<-lapply(Chr, function(f) gsub("5.1", "Chr10", f))
Chr<-lapply(Chr, function(f) gsub("6.1", "Chr11", f))
Chr<-lapply(Chr, function(f) gsub("7.1", "Chr12", f))
Chr

#Generate random alphanumeric code for each SNP
unique_id<-stri_rand_strings(n = 45, length = 4, pattern = "[A-Za-z0-9]")

DF7 <- as.data.frame(cbind(unique_id,Chr,Pos,Pos))



df_unlist<-function(df){
  
  df<-as.data.frame(df)
  
  nr<-nrow(df)
  
  c.names<-colnames(df)
  
  lscols<-as.vector(which(apply(df,2,is.list)==TRUE))
  
  if(length(lscols)!=0){
    
    for(i in lscols){
      
      temp<-as.vector(unlist(df[,i]))
      
      if(length(temp)!=nr){
        
        adj<-nr-length(temp)
        
        temp<-c(rep(0,adj),temp)
        
      }
      
      df[,i]<-temp
      
    } #end for
    
    df<-as.data.frame(df)
    
    colnames(df)<-c.names
  }
  return(df)
}


#Prepare an annotations file
DF7<-df_unlist(DF7)
write.table(DF7, file = "SNP2.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
```

# The SNP2.txt file 
![Screenshot from 2022-03-18 12-37-30](https://user-images.githubusercontent.com/93121277/158996508-4d24123b-489a-4c1a-aaf3-4f3d35fd749b.png)



# A chromosome file is needed to annotate the SNPs
![Screenshot from 2022-03-18 08-45-29](https://user-images.githubusercontent.com/93121277/158961964-d26eaec9-2cc2-48d4-8094-26d387797d20.png)


```r

Chr1Len <- 43270923
Chr2Len <- 35937250
Chr3Len <- 36413819
Chr4Len <- 35502694
Chr5Len <- 29958434
Chr6Len <- 31248787
Chr7Len <- 29697621
Chr8Len <- 28443022
Chr9Len <- 23012720
Chr10Len <- 23207287
Chr11Len <- 29021106
Chr12Len <- 27531856

ChrStartPos <- c(1,1,1,1,1,1,1,1,1,1,1,1)
ChrEndPos <- c(Chr1Len, Chr2Len, Chr3Len, Chr4Len, Chr5Len, Chr6Len, Chr7Len, Chr8Len, Chr9Len, Chr10Len, Chr11Len, Chr12Len)
length(ChrStartPos)
length(ChrEndPos)
Chr <- c("Chr1", "Chr2", "Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12")
DF8 <- as.data.frame(cbind(Chr,ChrStartPos,ChrEndPos))
DF8 <- df_unlist(DF8)
write.table(DF8, file = "chr2.txt", sep = "\t", row.names = FALSE, col.names = FALSE)


chromoMap("chr2.txt", "SNP2.txt", ploidy = 1)
getwd()


chromoMap("chr2.txt", "SNP2.txt",segment_annotation = T)

chromoMap("chr.txt", "SNP.txt", data_based_color_map = T, data_type = "categorical")

chromoMap("chr2.txt", "SNP2.txt",labels = T, label_font = 12, label_angle = -65)


```

# The annotated chromosomes for second related group  
![Screenshot from 2022-03-18 12-36-37](https://user-images.githubusercontent.com/93121277/158996526-077a0031-9567-4439-9652-66126193fe35.png)

```r
All_Int <- Reduce(intersect, list(v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18))
Pos<-lapply(All_Int, function(f) substr(f, 9:16, 21))
Chr<-lapply(All_Int, function(f) substr(f,6:8,8))


Chr<-lapply(Chr, function(f) gsub("6.1", "Chr1", f))
Chr<-lapply(Chr, function(f) gsub("7.1", "Chr2", f))
Chr<-lapply(Chr, function(f) gsub("8.1", "Chr3", f))
Chr<-lapply(Chr, function(f) gsub("9.1", "Chr4", f))
Chr<-lapply(Chr, function(f) gsub("0.1", "Chr5", f))
Chr<-lapply(Chr, function(f) gsub("1.1", "Chr6", f))
Chr<-lapply(Chr, function(f) gsub("2.1", "Chr7", f))
Chr<-lapply(Chr, function(f) gsub("3.1", "Chr8", f))
Chr<-lapply(Chr, function(f) gsub("4.1", "Chr9", f))
Chr<-lapply(Chr, function(f) gsub("5.1", "Chr10", f))
Chr<-lapply(Chr, function(f) gsub("6.1", "Chr11", f))
Chr<-lapply(Chr, function(f) gsub("7.1", "Chr12", f))


set.seed(99)
unique_id<-stri_rand_strings(n = 45, length = 4, pattern = "[A-Za-z0-9]")

DF8 <- as.data.frame(cbind(unique_id,Chr,Pos,Pos))
head(DF8)
DF8

df_unlist<-function(df){
  
  df<-as.data.frame(df)
  
  nr<-nrow(df)
  
  c.names<-colnames(df)
  
  lscols<-as.vector(which(apply(df,2,is.list)==TRUE))
  
  if(length(lscols)!=0){
    
    for(i in lscols){
      
      temp<-as.vector(unlist(df[,i]))
      
      if(length(temp)!=nr){
        
        adj<-nr-length(temp)
        
        temp<-c(rep(0,adj),temp)
        
      }
      
      df[,i]<-temp
      
    } #end for
    
    df<-as.data.frame(df)
    
    colnames(df)<-c.names
  }
  return(df)
}


#Prepare an annotations file
DF8<-df_unlist(DF8)
write.table(DF8, file = "SNP3.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(DF8, file = "SNP3.csv", sep = ",")
```

# A list of SNPs shared by both groups
![Screenshot from 2022-03-18 12-37-30](https://user-images.githubusercontent.com/93121277/158997314-0f03ca48-b721-4535-abbb-1bbdda8a5d5e.png)


# Annotated chromosome map of shared SNPs for both groups and all samples
```r
chromoMap("chr2.txt", "SNP3.txt",labels = T, label_font = 12, label_angle = -65)
```
![Screenshot from 2022-03-18 12-40-46](https://user-images.githubusercontent.com/93121277/158997409-59ea965e-2342-4917-b134-513b8279a322.png)

# Install the ape package and use the Gene Annotation file to determine if any of the shared SNPs on chromomsome 6 is in a gene

```r
install.packages("ape")
library(ape)
library(dplyr)
setwd("/home/michael/Desktop/GenomicVis/")
GFF <- read.gff(file = "GCF_001433935.1_IRGSP-1.0_genomic.gff")
GFF <- GFF %>% mutate_if(is.factor, as.character)
# Fiter GFF file by Gene and Chromosome 6
GFF_Gene <- filter(GFF,type == "gene" & seqid == "NC_029261.1")
#There are 2725 genes on Chromosome 6




#Call the SNP Table of Common SNPs
SNP_Table<- read.table(file = "SNP3.csv", header = TRUE, sep = ",")
#Filter by Chromosome 6
SNP_Table <- filter(SNP_Table, Chr == "Chr6")
#Order the Position by Descending
SNP_Table <- SNP_Table %>% arrange(desc(Pos), unique_id, Chr, Pos.1)
SNP_Table

GF <- filter(GFF_Gene, start > 18250740 & start < 19202286)
GF
GF[29,]
```
![Screenshot from 2022-03-28 08-34-01](https://user-images.githubusercontent.com/93121277/160339779-b7972c3f-ed37-4d23-bfc9-d649950fd89a.png)



![Screenshot from 2022-03-25 15-06-50](https://user-images.githubusercontent.com/93121277/160136349-9469533d-d8c6-4256-a2f1-58a5f7a82250.png)

# There are four particular SNPs shared by all Samples on Chromosome 6 intersecting Gene LOC107276770
# G ----> T

![Screenshot from 2022-03-28 08-36-39](https://user-images.githubusercontent.com/93121277/160340097-287d9e19-e1a1-47be-bdf2-f806717f76fb.png)


![Screenshot from 2022-03-29 10-03-37](https://user-images.githubusercontent.com/93121277/160563580-e61177ed-cbe4-4e74-8b35-c7681d4d6805.png)


