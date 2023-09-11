##### 1. Setup #####

# 1.1 Load necessary packages
library(limma)
library(dplyr)
library(ggplot2)
library(tidyr)
library(EnhancedVolcano)

# 1.2 Set directory
setwd('E:\\2023\\August 2023\\R data') 

# 1.3 Load counts matrix and extract relevant data
all_counts = read.csv('GeneCountsAllSamplesNegNormalised.csv', 
                      row.names=1, check.names=FALSE)

all_counts <- all_counts[-1,]   # remove empty row
all_counts <- log2(all_counts)  # transform to log2count
all_counts <- all_counts[       # remove POS and NEG control probes
  rownames(all_counts)
  [!startsWith(rownames(all_counts),'POS') 
    & !startsWith(rownames(all_counts),'NEG')],] 

# 1.4 Load the meta data and extract relevant data
pdata = read.csv('GeneCountLabelsAllSamples.csv', row.names=1, check.names=FALSE) 

pdata <- pdata[pdata$Bubble %in% c('OVA','MSA'), ] # only interested in the OVA and MSA data
pdata = pdata %>%
  filter(!(Label %in% c('6BT','7BM')))

##### 2. Relate the data columns to the experiment design #####
matching_columns <- intersect(colnames(all_counts), rownames(pdata)) 
counts <- all_counts[,matching_columns]                               # get relevant data
pdata <- pdata[matching_columns,]                                     # make sure counts and pdata are in the same order
colnames(counts) <- pdata$Label                                       # change column names to mouse label

##### 3. Create different "transformations" of cavitation dose #####
pdata = pdata %>%
  rename('cd' = 'Cavitation dose',
         'vol' = 'Tumour volume',
         'weight' = 'Tumour weight') %>%
  mutate(log_cd = log(cd*1e-17),
         log_cd_vol = log(cd*1e-17/vol),
         vol_weight_mean = (vol+weight)/2,
         cd_vol = cd/vol,
         cd_weight = cd/weight,
         cd_weight_sens = ifelse(weight<15,cd_weight/2,cd_weight),
         cd_vol_weight = cd/vol_weight_mean,
         US=`Ultrasound pressure`>0)


##### 4. Fit model #####

# 4.1 Create design matrix
design <- model.matrix(~ pdata$cd + pdata$Bubble
                       + pdata$Bubble:pdata$cd, data = counts)

# 4.2 Fit model
fit <- lmFit(counts, design, weights=arrayWeights(counts,design))

colnames(coef(fit))

# 4.3 Create contrast matrix
fit = limma::contrasts.fit(
  fit,
  cbind(
    c(0,1,0,0), # US on MSA
    c(0,1,0,1), # US on OVA
    c(0,0,0,1) # interaction (effect of type of bubble on effect of ultrasound)
  ))

# 4.4 compute statistics from fit
fit_eBayes <- eBayes(fit, trend = TRUE, robust = T)


# 4.5 Generate tables of top genes
MSA_topTable = topTable(fit_eBayes, coef = 1, number = 'inf', 
                        p.value=0.05, sort.by = "logFC", adjust.method="BH")
OVA_topTable = topTable(fit_eBayes, coef = 2, number = 'inf', 
                        p.value=0.05, sort.by = "logFC", adjust.method="BH")
MSAandOVA_topTable = topTable(fit_eBayes, coef = 3, number = 'inf', 
                              p.value=0.05, sort.by = "logFC", adjust.method="BH")

##### 5. Create volcano plots #####

# 5.1 Define a threshold for the LogFC boundary
threshold = quantile(pdata$cd,0.9)
sort(pdata$cd)

# 5.2 Generate table for data to be plotted 
# -- coef needs to be changed depending on variable of interest
all_genes <- topTable(fit_eBayes, coef = 2, number = Inf, sort.by = "P", adjust.method="BH") # Change coef based on MSA, OVA, OVA vs MSA

# 5.3 Table for genes to be labeled in the volcano plot
# -- coef needs to be changed depending on variable of interest
top_genes <- topTable(fit_eBayes,coef = 2, sort.by='logFC', lfc=log2(2)/threshold, # 2-fold change for 5 nJ difference
                      p.value=0.05,n=Inf)

# 5.4 Volcano plot
EnhancedVolcano(all_genes,
                lab = rownames(all_genes),
                selectLab = rownames(top_genes),
                labSize=3,
                x='logFC',
                y='P.Value',
                FCcutoff=1/threshold,
                pCutoff=max(all_genes[all_genes$adj.P.Val<0.05,]$P.Value),
                xlab=bquote(~Log[2]~'FC/unit cavitation dose'),
                xlim=c(-0.01,0.01),
                ylim=c(0,8),
                legendPosition = 'none',
                max.overlaps=Inf,
                title = "OVA MBs",
                subtitle = bquote(italic('')),
)

##### 6. Plot venndiagram #####
dt <- decideTests(fit_eBayes)
vennDiagram(dt[,1:2],names = c('MSA MBs','OVA MBs'),circle.col = c('turquoise', 'salmon'))



##### 7. Gene enrichment analysis #####

# 7.1 Create csv file of all analysed genes
ref_genes = rownames(counts)
write.table(matrix(ref_genes,ncol=1),
            file.path('ReferenceGenes.csv'),row.names=F,
            col.names=F,quote=F)

# 7.2 Create csv file of all genes with adj. P-value less than 0.05
genes_overrep = all_genes %>% 
  filter(adj.P.Val<0.05) %>%
  dplyr::select(logFC) 

write.table(genes_overrep,
            file.path('GenesOVAOverrep.csv'),row.names=T,
            col.names=F,quote=F,sep='\t')

# 7.3 Upload the csv files to pantherdb.org to do the GE analysis

# 7.4 Extract the relevant data from the analysis file
overrep_raw=read_json(file.path(basedir,'analysis.json'))
overrep_df = data.frame()
for(group in overrep_raw$overrepresentation[[3]]) {
  overrep_df = rbind(
    overrep_df,
    data.frame(
      number_in_list = group$input_list$number_in_list,        
      fdr = group$input_list$fdr,
      pValue = group$input_list$pValue,    
      plus_minus = group$input_list$plus_minus,  
      fold_enrichment = group$input_list$fold_enrichment,  
      expected = group$input_list$expected,  
      #term_id = group$term$id,
      term_label = group$term$label
    )
  )
}

# 7.5 Plot the biological processes
ggplot(overrep_df %>% filter(pValue<0.005, fold_enrichment > 1),
       aes(x=fold_enrichment,
           y=forcats::fct_reorder(as.factor(stringr::str_to_sentence(term_label)), fold_enrichment),
           size=number_in_list,color=pValue)) +
  geom_point() +
  #ggrepel::geom_text_repel(aes(label=number_in_list)) +
  #scale_size_continuous(limits=c(8,25),breaks=c(10,15,20,25)) +
  scale_color_gradient(low='blue',high='red') +
  xlab('Fold enrichment') +
  ylab('GO Biological Process') +
  labs(size='Number',colour='p value') +
  #facet_wrap(~group,ncol=1) +
  theme_minimal()

sum(p.adjust(overrep_df$pValue,method='fdr') < 1)
