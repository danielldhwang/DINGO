#Input format:
#Tab-separated .gz file with: CHR	SNP BP	A1	A2	FREQ	BETA	SE	P   N

library(data.table)
library(lubridate)

args=commandArgs(trailingOnly=T)
fetal_file <- args[1] #fetal gwas file name
maternal_file <- args[2] #maternal gwas file name
int <- as.numeric(args[3]) #LD score rg intercept
out_file <- args[4] #output file name

message(date(),": Reading GWAS files")

fetal <- fread(fetal_file, header=T) 
maternal <- fread(maternal_file, header=T) 

fetal$A1 <- toupper(fetal$A1)
fetal$A2 <- toupper(fetal$A2)
maternal$A1 <- toupper(maternal$A1)
maternal$A2 <- toupper(maternal$A2)

message(date(),": Read GWAS files")

colnames(fetal) <- c("CHR", "SNP","BP", "ea_fetal", "nea_fetal", "eaf_fetal", "Beta_fetal", "SE_fetal","p_fetal","n_fetal")
colnames(maternal) <- c("CHR", "SNP","BP","ea_maternal","nea_maternal","eaf_maternal","Beta_maternal","SE_maternal","p_maternal","n_maternal")

#Quick merge
setkey(fetal, "SNP")
setkey(maternal, "SNP")
data<-fetal[maternal]
data<-data[complete.cases(data), ]

message(date(),": Merged GWAS files")
message(date(),": Number of SNPs: ", nrow(data))

#Check the effect alleles are the same
#If alleles do not match, switch maternal alleles and change the sign of the beta coefficient
data$Beta_maternal_new <- ifelse(data$ea_fetal == data$nea_maternal, data$Beta_maternal*(-1), data$Beta_maternal)
data$ea_maternal_new <- ifelse(data$ea_fetal == data$nea_maternal, data$nea_maternal , data$ea_maternal)
data$nea_maternal_new <- ifelse(data$ea_fetal == data$nea_maternal, data$ea_maternal , data$nea_maternal)
data$eaf_maternal_new <- ifelse(data$ea_fetal == data$nea_maternal, 1-data$eaf_maternal , data$eaf_maternal)


data <- subset(data, data$ea_fetal == data$ea_maternal_new & data$nea_fetal == data$nea_maternal_new)
message(date(),": Number of SNPs after removing non matching alleles: ", nrow(data))

#Check strandness / palindromic SNPs 
#Identify palindromic SNPs; 
data$palindromic <- ifelse(data$ea_maternal_new=="A" & data$nea_maternal_new =="T",1,
                           ifelse(data$ea_maternal_new=="T" & data$nea_maternal_new =="A",1,
                                  ifelse(data$ea_maternal_new=="C" & data$nea_maternal_new =="G",1,
                                         ifelse(data$ea_maternal_new=="G" & data$nea_maternal_new =="G",1,0)))) 

#Identify the palindromic SNPs + maternal and fetal SNPs have different minor allele; 
data$palindromic2 <- ifelse(data$palindromic==1 & data$eaf_maternal_new > 0.5 & data$eaf_fetal < 0.5, 1,
                            ifelse(data$palindromic==1 & data$eaf_maternal_new < 0.5 & data$eaf_fetal > 0.5, 1,0))

#Identify the palindromic SNPs with different MAF from maternal and fetal + the difference between eaf_maternal and eaf_fetal > 0.2; n=26
data$palindromic3 <- ifelse(data$palindromic2==1 & abs(data$eaf_maternal_new - data$eaf_fetal)>0.2, 1, 0)

data <- subset(data, palindromic3==0) 
message(date(),": Number of SNPs after removing palindromic SNPs with different MAF from maternal and fetal: ", nrow(data))

data$Beta_maternal <- data$Beta_maternal_new
data$ea_maternal <- data$ea_maternal_new
data$nea_maternal <- data$nea_maternal_new
data$eaf_maternal <- data$eaf_maternal_new

#Identify SNPs with standard errors = 0. SE is used as denominator in DINGO anlysis ana cannot be 0.
data <- subset(data, SE_fetal !=0 & SE_maternal !=0)
message(date(),": Number of SNPs after removing SNPs with SE=0 from maternal and fetal: ", nrow(data))

message(date(),": Cleaned GWAS files")

#NB Insert intercept from LD score regression
data$fetal_beta_adjusted <- ((4/3)*data$Beta_fetal) - ((2/3)*data$Beta_maternal)
data$maternal_beta_adjusted <- ((4/3)*data$Beta_maternal) - ((2/3)*data$Beta_fetal)

data$fetal_var_adjusted <- ((16/9)*(data$SE_fetal)^2) + ((4/9)*(data$SE_maternal)^2) - ((16/9)*int*data$SE_fetal*data$SE_maternal)
data$maternal_var_adjusted <- ((16/9)*(data$SE_maternal)^2) + ((4/9)*(data$SE_fetal)^2) - ((16/9)*int*data$SE_fetal*data$SE_maternal)

data$covar <- ((20/9)*int*data$SE_fetal*data$SE_maternal) - ((8/9)*(data$SE_fetal)^2)-((8/9)*(data$SE_maternal)^2)

effects <- matrix(nrow = 1, ncol = 2)
sigma <- matrix(nrow = 2, ncol = 2)

chisq_2df <- vector(length= nrow(data))

message(date(), ": Computing 2df pvalue for all SNPs")

pb <- txtProgressBar(min = 0,
                     max = nrow(data),
                     style = 3,
                     width = 100, # Needed to avoid multiple printings
                     char = "=")

for (i in 1:nrow(data)) {
  effects <- matrix(c(data$fetal_beta_adjusted[i], data$maternal_beta_adjusted[i]), nrow=1, ncol=2, byrow=TRUE)
  sigma <- matrix(c(data$fetal_var_adjusted[i], data$covar[i], data$covar[i], data$maternal_var_adjusted[i]), nrow=2, ncol=2, byrow=TRUE)
  chisq_2df[i] <- effects%*%solve(sigma) %*%t(effects)
  setTxtProgressBar(pb, i)
}

pval_2df <- pchisq(chisq_2df, df = 2, ncp = 0, lower.tail = FALSE, log.p = FALSE)
chisq_fetal_adj <- (data$fetal_beta_adjusted^2)/(data$fetal_var_adjusted)
chisq_maternal_adj <- (data$maternal_beta_adjusted^2)/(data$maternal_var_adjusted)
pval_fetal_adj <- pchisq(chisq_fetal_adj, df = 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
pval_maternal_adj <- pchisq(chisq_maternal_adj, df = 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)

data$pval_fetal_adj <- pval_fetal_adj
data$fetal_se_adjusted <- sqrt(data$fetal_var_adjusted)
data$pval_maternal_adj <- pval_maternal_adj
data$maternal_se_adjusted <- sqrt(data$maternal_var_adjusted)
data$chi_2df <- chisq_2df
data$pval_2df <- pval_2df
data$pval_best_1df <- ifelse(data$p_fetal < data$p_maternal, data$p_fetal, data$p_maternal)



fwrite(data, file=args[4], quote=F, col.names=T, sep="\t", row.names=F)

message(date(), ": Finished!")