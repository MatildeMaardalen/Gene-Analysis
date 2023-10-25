##### 1. Setup #####
# 1.1 Load necessary packages
library(limma)
library(dplyr)
library(ggplot2)
library(tidyr)
library(progress)

# 1.2 Set directory
setwd('D:\\2023\\August 2023\\R data') 

# 1.3 Load counts matrix and extract relevant data
all_counts = read.csv('GeneCountsAllSamplesNegNormalised.csv', 
                      row.names=1, check.names=FALSE)

ctrl_probes <- all_counts[grepl("^POS|^NEG", rownames(all_counts)),]  # positive and negative control probes
housekeeping = c('Abcf1','Alas1','Edc3','Eef1g','Eif2b4','G6pdx',     # housekeeping genes
                 'Hdac3','Hprt','Nubp1','Oaz1','Polr1b','Polr2a',
                 'Ppia','Rpl19','Sap130','Sdha','Sf3a3','Tbd',
                 'Tubb5')

all_counts <- all_counts %>%
  filter(
    rownames(all_counts) != "",         # remove empty rows
    !(rownames(all_counts) %in%         # remove control probes
        rownames(ctrl_probes)),         
    !(rownames(all_counts) %in%         # remove housekeeping genes
        housekeeping)
  )%>%
  mutate(across(everything(), 
                ~ log2(.)))             # transform all columns to log2

# 1.4 Load the meta data and extract relevant data
pdata = read.csv('GeneCountLabelsAllSamples.csv', row.names=1, check.names=FALSE) 

pdata = pdata %>%
  filter(
    Bubble == 'OVA',                    # only interested in the OVA data
    !(Label %in% c('6BT','7BM')),       # outliers/missing data points
    `Ultrasound pressure` %in% c(0,2)   # only look at the 0 MPa and 2 MPa US group
  ) %>%        
  rename('cd' = 'Cavitation dose')      # to simplify

##### 2. Relate the data columns to the experiment design #####
matching_columns <- intersect(colnames(all_counts), rownames(pdata)) 
counts <- all_counts[,matching_columns]                               # get relevant data
pdata <- pdata[matching_columns,]                                     # make sure counts and pdata are in the same order
colnames(counts) <- pdata$Label                                       # change column names to mouse label

##### 3. Fit model #####
# 3.1 Create design matrix
design <- model.matrix(~pdata$cd, data = counts)

# 3.2 Fit model
fit <- lmFit(counts, design, weights=arrayWeights(counts,design))

colnames(coef(fit))

# 3.3 compute statistics from fit
fit_eBayes <- eBayes(fit, trend = TRUE, robust = T)

##### 4. Get values for power simulation #####

# 4.1 Plot Normal Q-Q plot to check the assumption that the residuals are normally distributed
baseline_counts= as.numeric(
  pracma::Reshape(
    apply(counts[,which(pdata$`Ultrasound pressure`==0)],1,       
          function(x) scale(x,center=TRUE,scale=FALSE)),1,752*6))
qqnorm(baseline_counts)

# 4.2 Compute mean and standard deviation of cavitation dose
mean_cd = mean(pdata[pdata$`Ultrasound pressure` == 2,]$cd)
sd_cd = sd(pdata[pdata$`Ultrasound pressure` == 2,]$cd)

# 4.3 Get standard deviations of the residuals
sd_residuals = sd(residuals(fit_eBayes,counts))

threshold = quantile(pdata$cd,0.90)

# 4.4 Simulate log2 counts and cavitation dose values
pressure_to_cd = function(p) {
  ifelse(p==0,2.4,rnorm(length(p),mean_cd,sd_cd))
}

cd_to_log2counts = function(cd,slope) {
  slope*cd + rnorm(length(cd),0,sd_residuals)
}

do_experiment_limma = function(p,n_mice_ctrl,n_mice_trt,n_genes,n_deg,pb) {
  pressure = c(rep(0,n_mice_ctrl),rep(p,n_mice_trt))
  pdata_sim = data.frame(cd = pressure_to_cd(pressure))
  counts_sim = data.frame()
  for (i in 1:n_genes) {
    if (i < n_deg + 1) {
      slope = log2(2)/threshold
    } else {
      slope = 0
    }
    counts_sim = rbind(
      counts_sim,
      cd_to_log2counts(pdata_sim$cd,slope)
    )
  }
  design_sim <- model.matrix(~pdata_sim$cd, data = counts_sim)
  fit_sim <- lmFit(counts_sim, design_sim, weights=arrayWeights(counts_sim,design_sim))
  fit_eBayes_sim <- eBayes(fit_sim, trend = TRUE, robust = T)
  
  # focuses on power to detect differential expression of any of the DEG
  pb$tick()
  return(topTable(fit_eBayes_sim, 2, adjust.method = 'BH', sort.by='none')$adj.P.Val[1]<0.05)
}

n_sim = 10000 
exp_design = expand.grid(pressure = 2, n_mice_ctrl = c(6,8), tot_mice = c(12,16,20,24)) %>%
  mutate(n_mice_trt = tot_mice - n_mice_ctrl)
pb=progress_bar$new(total=nrow(exp_design)*n_sim, format='[:bar] :percent eta: :eta')
results = exp_design %>%
  rowwise() %>%
  mutate(power = 100*sum(replicate(n_sim, do_experiment_limma(pressure, n_mice_ctrl, n_mice_trt, 752, 200, pb)))/n_sim) %>%
  ungroup()

ggplot(results,aes(x=tot_mice,y=power,color=as.factor(n_mice_ctrl),group=n_mice_ctrl)) +
  geom_point() +
  geom_line() +
  xlim(0,25) +
  ylim(0,100) +
  geom_hline(yintercept=80) +
  #facet_wrap(~n_mice_ctrl) +
  theme_minimal()

