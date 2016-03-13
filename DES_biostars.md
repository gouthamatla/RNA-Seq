```
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggrepel")

# matrix of 22348 rows and 60 columns. Treated and untreated samples are adjacent to each other.
# treated untreated treated untreated treated untreated . . . . 

x <- read.table("filt_counts.txt", header=T, row.names=1)

#scrap the column names and give numbers in x.1 and x.2 format to indicate treated and untreated respectively.
a <- seq(1.1, 30.1, by=1)
b <- seq(1.2, 30.2, by=1)
colnames(x) <- sort(as.numeric(c(a,b)))

#create the paired design.
subjects=factor(c(rep(1:30, each=2)))
treat <- as.factor(rep(c("treat","untreat"),30))
design <- model.matrix(~subjects+treat)

colData <- data.frame(colnames(x),subjects=subjects, treat=treat, row.names=1)

dds <- DESeqDataSetFromMatrix(countData = x, colData = colData, design = ~ subjects + treat)
design(dds) <- formula(~ subjects + treat)
dds <- DESeq(dds)
resdds <- results(dds)
resdds=resdds[order(resdds$padj),]

summary(resdds)

out of 22348 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 5724, 26% 
LFC < 0 (down)   : 5643, 25% 
outliers [1]     : 0, 0% 
low counts [2]   : 0, 0% 
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)         


data <- plotPCA(vsd, intgroup=c("treat", "subjects"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=treat)) +
   geom_point(size=3) +
   geom_text_repel(aes(label=row.names(data))) +
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance"))
```
![PCA plot of VST of all genes]
(http://i.imgur.com/ZSTOhNX.png?1)

```
#PCA plot of VSD of Differentially expressed genes at 0.05. 

de <- rownames(resdds[na.omit(resdds)$padj<0.05,])
data <- plotPCA(vsd[de,], intgroup=c("treat", "subjects"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=treat)) +
   geom_point(size=3) +
   geom_text_repel(aes(label=row.names(data))) +
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance"))
```
![PCA plot of VST of DE genes]
(http://i.imgur.com/TXZSN41.png)

```
# create the Fold change matrix
norm <- assay(normTransform(dds))
i <- seq.int(1L,60,by = 2L)

#substract the every treat from corresponding untreated samples. This creates a matrix of 22348 genes and  30 columns.
norm.fc <- norm[,i]-norm[,i+1]

#Create a DESeq2 object for plotting.
subjects=factor(c(rep(1:30)))
colData <- data.frame(colnames(norm.fc),subjects=subjects, row.names=1)
se <- SummarizedExperiment(norm.fc,colData=DataFrame(colData))
data <- plotPCA(DESeqTransform(se), intgroup=c("subjects"), returnData=TRUE)

#Plot a PCA for FC matrix.
ggplot(data, aes(PC1, PC2)) +
   geom_point(size=3) +
   geom_text_repel(aes(label=data$group)) +
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar[2],"% variance"))
```

![PCA plot of Fold Change Matrix]
(http://i.imgur.com/QAYH6EI.png)
   
   
   
