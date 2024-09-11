suppressMessages(library(dplyr))
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(tibble))
#install.packages("svMisc")

#setwd("./Desktop/CancerResponseData/src/")
#getwd()
args = commandArgs(trailingOnly=TRUE)

drug = args[1]
tissue = args[2]
n_iters = as.numeric(args[3])
countFile = args[4]
colDataFile = args[5]


counts = as.matrix(read.csv(countFile,row.names = 'Gene_ID',check.names=FALSE))
counts = round(counts)

mode(counts) <- "integer"


col_data = read.csv(colDataFile)

col_data$Response <- factor(col_data$Binned_Response)


N = length(colnames(counts))




arg_echo = paste("\nWorking on:\n\tdrug:\t",drug,"\n\ttissue:\t", tissue,"\n\titerations:\t",n_iters,"\n")
cat(arg_echo)



fdr_thresholds <- c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1,0.15)

names(fdr_thresholds) <- c('narrow', "very.tight", "tight", 
                           "medium", "loose", "very.loose","easy")


count_table <- matrix(0,nrow(counts),length(fdr_thresholds))
colnames(count_table) <- names(fdr_thresholds)
rownames(count_table) <- rownames(counts)



set.seed(1234)
pb = txtProgressBar(0,n_iters,style = 3)

for(i in seq(n_iters)){
  setTxtProgressBar(pb,i)
  
  samples = as.vector(sample(colnames(counts),
                           size = round(0.75*N), 
                           replace = FALSE))
  
  sample_data = col_data[match(samples,col_data$Run_ID),]
 
  sample_counts = counts[,samples]
  
  rownames(sample_data) <- sample_data[,1]
  sample_data[,1] <- NULL
  
  
  dds <- DESeqDataSetFromMatrix(countData = sample_counts,
                                colData = sample_data,
                                design = ~ Response)
   
  
  smallestGroupSize <- 3
  keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
  dds <- dds[keep,]
  
  dds <- estimateSizeFactors(dds)
  
  dds <- DESeq(dds,quiet = TRUE)
  
  res <- results(dds)
  res <- res[!is.na(res$padj),]
  
  for (thresh.name in names(fdr_thresholds)){
    thresh.name
    temp <- res[res$padj <= fdr_thresholds[thresh.name],]

    for (gene in rownames(temp)){
      count_table[gene,thresh.name] <-  count_table[gene,thresh.name]+1
    }
  }

}


res = data.frame()

for(thresh.name in names(fdr_thresholds)){
  count_vec <- count_table[,thresh.name]
  count_vec <- count_vec[count_vec>0]
  df = data.frame(Gene = names(count_vec),
                  Count = count_vec, 
                  Thresh.Name = rep(thresh.name,length(count_vec)),
                  Thresh.Value = rep(fdr_thresholds[thresh.name],length(count_vec)))
  rownames(df) <- NULL
  res = rbind(res,df)
 
  
}


write.csv(x= res,file = paste0("../data/processed/genesets/",drug,"_",tissue,"_DE.csv"))




