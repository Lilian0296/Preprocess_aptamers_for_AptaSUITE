library(Biostrings)
library(fuzzywuzzyR)
library(QuasR)
library(tidyr)
library(dplyr)
library("optparse")


option_list = list(make_option("--output", type = "character", default = NULL, help  = "The path of output files"),
                   make_option("--data", type = "character", default = NULL, help  = "The path of data folder (fastq files)"))  


args <- parse_args(OptionParser(option_list=option_list))

pathout<-args$outpath

dir.create(paste(pathout,"/raw",sep=""))
dir.create(paste(pathout,"/Clean_data",sep=""))

# Merge fastq -------------------------------------
message("Merge fastq files")

datafolder <-args$data
files=list.files(datafolder)
message(files," will be preocessed")

merge_fastq=DNAStringSet()

for (i in 1:length(files)){
  if (endsWith(files[i], "fastq.gz")){
    data_path=paste(datafolder,"/",files[i],sep="")
    print(paste("read ",files[i],sep=""))
    data=readDNAStringSet(data_path,format="fastq",with.qualities=FALSE)
    merge_fastq=append(merge_fastq,data)
    #append(path)
  }else{
    message("Skip non-fastq files")
  }
}

df <- as.data.frame(merge_fastq)
# Read datatable----------------------------

message("Read Primers data")
Primer <- read.csv(paste(pathout,"/Primer.csv",sep = ''))
message("Read Barcodes data")
Barcodes <- read.csv(paste(pathout,"/Barcodes.csv",sep = ''))

Primer=Primer %>% rowwise() %>% mutate(Sequences_reverse=as.character(reverseComplement(DNAString(Sequences))))
Barcodes=Barcodes%>% rowwise() %>% mutate(Search_seed=paste(Sequences,Primer$Sequences[1],sep=("")))

message("Process forward and reverse reads")

fd <- data.frame(agrep(Primer["Sequences"][1,],df$x, max.distance = 2, value = TRUE))
rv <- data.frame(agrep(Primer["Sequences_reverse"][1,], df$x, max.distance = 2, value = TRUE))
rv <- data.frame(reverseComplement(DNAStringSet(rv$agrep.Primer..Sequences_reverse...1.....df.x..max.distance...2..)))
colnames(fd)<-"sequences"
colnames(rv)<-"sequences"
df=rbind(fd,rv)

message("Split reads based on barcodes")

for (i in 1:length(Barcodes$Search_seed)){
  message(paste("Process ",Barcodes$Name[i],sep=""))
  barcodes_data <-agrep(Barcodes$Search_seed[i],df$sequences, max.distance = 1, value = TRUE)
  data=DNAStringSet(barcodes_data)
  filename_fasta <- paste(pathout,"/raw/(+)_",Barcodes$Name[i],".fasta",sep = '')
  writeXStringSet(data,filename_fasta,format='fasta')
  filename_fasta_removeprimers <- paste(pathout,"/raw/(+)_removeprimers_",Barcodes$Name[i],".fasta",sep = '')
  removeprimers <-QuasR::preprocessReads(filename_fasta,
                                         outputFilename=filename_fasta_removeprimers,
                                         minLength=30,Lpattern = Primer$Sequences[1],Rpattern = Primer$Sequences[2],
                                         truncateEndBases = 3)
  print(removeprimers)
  clean_data=data.frame(readDNAStringSet(filename_fasta_removeprimers,format="fasta",with.qualities=FALSE))
  colnames(clean_data)<-"sequences"
  random_clean=data.frame(clean_data[-grep(Primer$Sequences[1],clean_data$sequences),])
  random_clean=random_clean[-grep(Primer$Sequences[2],random_clean$clean_data..grep.Primer.Sequences.1...clean_data.sequences...),]
  random_clean=DNAStringSet(random_clean)
  Fakeid=sprintf("ID%07d", seq_len(length(random_clean)))
  names(random_clean) <-Fakeid
  writeXStringSet(random_clean,paste(pathout,"/Clean_data/(+)_",Barcodes$Name[i],".fastq",sep = ''),format='fastq')
}






# Supplementary-----------------------------------------------------------------
## create dataframe------
#Primer=data.frame(Primer=c("Primer_forward","Primer_reverse"),
#                  Sequences=c("ACGCTCGGATGCCACTACAG","CTCATGGACGTGCTGGTGAC"))


#Barcodes=data.frame(Name=c("B1","B2","B3","B4"),
#                    Sequences=c("CATGTCA",
#                                "GTAAAGTCAA",
#                                "CGGT",
#                                "AACAGTTG"))

#write.csv(Primer,"/Primer.csv",row.names=FALSE)
#write.csv(Barcodes,"/Barcodes.csv",row.names=FALSE)
