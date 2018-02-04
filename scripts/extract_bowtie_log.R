#Extract read mapping information from bowtie log files

samples = read.table("barcodes.tsv", stringsAsFactors = F, header = F, sep="\t")

results_number=NULL
results_perc=NULL

#This loop inputs each file using read.delim and extracts the number of reads from each row.
#results_number stores the number of reads
#results_perc stores the same as a percentage of the initial reads
for(j in 1:nrow(samples))
{
  filename=paste("logs/bowtie/bowtie-align-",samples[j,1],".log",sep="")
  logfile = read.delim(filename, stringsAsFactors = F)
  number = NULL
  perc = NULL
  for (i in 1:4)
  {
    end = regexpr("\\(",logfile[i,])-1
    number = c(number,as.numeric(substring(logfile[i,],1,end)))
    perc = c(perc,round(number[i]/number[1]*100,2))
  }
  number=c(samples[j,1], number)
  perc=c(samples[j,1], perc)
  
  results_number=rbind(results_number,number)
  results_perc=rbind(results_perc,perc)
}
colnames(results_number)= c("Sample","Total", "Unmapped", "Mapped exactly 1 time", "Mapped >1 time")
colnames(results_perc) = colnames(results_number)

write.table(results_number, snakemake@output[["number"]], row.names=F, quote=F, sep="\t")
write.table(results_perc, snakemake@output[["percentage"]], row.names=F, quote=F, sep="\t")