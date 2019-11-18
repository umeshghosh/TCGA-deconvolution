
### Compare DemixT purity estimates with Tumeric and other methods

# Download TCGA LUAD rna seq data, and unzip
# https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/HiSeqV2.gz

# Download TCGA HNSC rna seq data, and unzip
# https://tcga.xenahubs.net/download/TCGA.HNSC.sampleMap/HiSeqV2.gz

### Read purity estimates of other methods
methods.purity = read.csv('data/purity_separate.csv',row.names = 1)
row.names(methods.purity) = make.names(row.names(methods.purity))

### Read gene expression data
ctype = 'LUAD'
x = read.table('~/Downloads/HiSeqV2.LUAD', stringsAsFactors = F)
nsamples = 30
#ctype = 'HNSC'
#x = read.table('~/Downloads/HiSeqV2.HNSC', stringsAsFactors = F)

### Format data
rownames(x) = x[,1]; colnames(x) = make.names(unlist(x[1,])); x = data.matrix(x[-1,-1])
# move to linear scale, log2(x+1)
x=2**x-1
# filter non and lowly-expressed genes
x = x[apply(x,1,median)>2,] # median > 2 FPKM

tumor = (x[,grep('.01$',colnames(x))])[,1:nsamples] # 515 samples
normal = (x[,grep('.11$',colnames(x))])[,1:nsamples] # 59 samples

### Demix
#devtools::install_github("wwylab/DeMixTallmaterials/DeMixT_0.2.1")
library(DeMixT)
res <- DeMixT(data.Y = tumor, data.comp1 = normal)
demix.purity = res$pi[1,]

### Compare purity estimates

all.purity = t(t(1-demix.purity))
all.purity = cbind(all.purity,methods.purity[names(demix.purity),'Purity'])
all.purity = cbind(all.purity,methods.purity[names(demix.purity),'AbsCN.seq'])
all.purity = cbind(all.purity,methods.purity[names(demix.purity),'ASCAT'])
all.purity = cbind(all.purity,methods.purity[names(demix.purity),'ESTIMATE'])
all.purity = cbind(all.purity,methods.purity[names(demix.purity),'PurBayes'])
all.purity = cbind(all.purity,methods.purity[names(demix.purity),'CPE'])
all.purity[all.purity>0.98] = NA
colnames(all.purity) = c('DeMix','TUMERIC','AbsCN-seq','ASCAT','ESTIMATE','PurBayes','CPE')

boxplot(all.purity,ylab='estimated purity',main=ctype)
cor.pur = cor(all.purity,use='complete.obs')
cor.pur[cor.pur>0.98] = NA
boxplot(cor.pur,ylab='Correlation (r) with other methods',main=ctype)
pairs(all.purity,xlim=c(0,1),ylim=c(0,1))

