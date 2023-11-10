##拆分下载好的肠道菌群数据
set.seed(123)
load("D:/研究课题/肠道菌群/Gut Microbiota-abdominal aneurysm/拆分数据/exp_dat.Rdata")  ##提前运行了


setwd("D:/研究课题/肠道菌群/Gut Microbiota-abdominal aneurysm")#设置工作空间
Q1<-read.csv("mibiogen_tophits_2023-06-05_11_35_04.csv")#读入数据
name_city<-unique(Q1[,2])#根据第四列（县级行政单位），数据去重，得到各个县的字符串
n<-length(name_city)#获得数据中县级行政单位的个数
out1<-as.character()#定义out1为字符串类型
out_filePath<-as.character()#定义out_filePath<为字符串类型

dir.create("拆分数据")
setwd("D:/研究课题/肠道菌群/Gut Microbiota-abdominal aneurysm/拆分数据")
for (j in 1:n) 
{
  print(name_city[j]);
  outPath = "D:/研究课题/肠道菌群/Gut Microbiota-abdominal aneurysm/拆分数据"   ##输出路径
  out1[j]=paste(outPath,name_city[j],sep='/') ##输出路径名
  out_filePath[j]=paste(name_city[j],".csv",sep='') #最终的输出文件整个路径
  Q2<-subset(Q1, Q1[, 2] == name_city[j])#根据县分类，提取县名为name_city[j]的数据单独成表
  write.csv(Q2,file=out_filePath[j])#输出表格
}
##输出文件

#* ---install or library-------------------------------------
if (!require("devtools")) { install.packages("devtools") } else {}
if (!require("data.table")) { install.packages("data.table") } else {}
if (!require("TwoSampleMR")) { devtools::install_github("MRCIEU/TwoSampleMR") } else {}
#install.packages("remotes")
#remotes::install_github("MRCIEU/MRInstruments")
library(MRInstruments)
library(plyr)
library(dplyr)  


FileNames<-list.files(paste0(getwd()),pattern=".csv")

#* ---Set working directory-------------------------------------
`%+%` <- function(x,y) paste0(x,y)
dir.create(getwd()%+%"/data")
path_data <- getwd()%+%"/data"
dir.create(getwd()%+%"/result")
path_res <- getwd()%+%"/result"

for(i in c(1:length(FileNames))){
  #d1<- try(fread(paste0(getwd(),"/",FileNames[1]),sep = "\t"),silent = T)
  d1<- try(read.csv(paste0(getwd(),"/",FileNames[i])))
  d2<-subset(d1,d1$P<1e-5)
  d3 <- d2[,c(3,6,7,8,9,10,12,13)]
  colnames(d3) <- c("bac","rsID","ref.allele","eff.allele","beta","SE","P.weightedSumZ","N")
  write.table(d3,paste0(path_data,"/",FileNames[i]),sep=",", col.names = T, row.names = F, quote = F)}
######################上面的循环是需要下载全部的肠道菌群数据，
#d2 <-  read_csv("mibiogen_tophits_2023-06-05_11_35_04.csv") ##这是在线下载所有菌群小于1e-5.
#d3 <- d2[,c(2,5,6,7,8,9,11,12)]
#colnames(d3) <- c("bac","rsID","ref.allele","eff.allele","beta","SE","P.weightedSumZ","N")

#* ---export IVs-------------------------------------
if(F){
  exp_dat <- list()
  ex_pore <- c()
  for(i in c(1:length(FileNames))){
    IV <- fread(paste0(path_data,"/",FileNames[i]))
    IV$PHENO <- FileNames[i]
    IV1<-format_data(IV,
                     type="exposure",
                     phenotype_col = "PHENO",
                     snp_col = "rsID",
                     beta_col = "beta",
                     se_col = "SE",
                     pval_col = "P.weightedSumZ",
                     samplesize_col = "N",
                     effect_allele_col = "eff.allele",
                     other_allele_col = "ref.allele")
    IV1 <- clump_data(IV1,clump_kb = 10000,clump_r2 = 0.001)
    exp_dat[[i]] <- IV1
    ex_pore<-c(ex_pore,FileNames[i])
  }
  save.image("exp_dat.Rdata")
  load("exp_dat.Rdata")
}
#######################################################

#IEU open gwas在线数据
#allSNP  <- do.call(rbind, exp_dat)
#out_data <- extract_outcome_data(allSNP$SNP,"finn-b-I9_CEREBATHER")  ##"Cerebral atherosclerosis"


#######################################################
##Outcome data(本地数据)
 # dataname1="D:/Codeprogram/overall.tsv"
library("data.table")
  
GWAS_1 <- fread('finngen_R9_I9_ABAORTANEUR.gz', sep = '\t', header = TRUE)
  allSNP  <- do.call(rbind, exp_dat)
  GWAS_2 <- subset(GWAS_1,GWAS_1$rsids %in% allSNP$SNP)
  rm(GWAS_1)
  GWAS_2$PHENO<-"abdominal aneurysm"
  out_data     <- format_data(GWAS_2,
                              type="outcome",
                              phenotype_col = "PHENO",
                              snp_col = "rsids",
                              beta_col = "beta",
                              se_col = "sebeta",
                              eaf_col = "af_alt",    ##"effect_allele_frequency",
                              pval_col = "pval",
                              effect_allele_col = "alt",    ##"effect_allele",
                              other_allele_col = "ref",         ###"other_allele",
                              )

#################################################################
out_dat <- list()
out_dat[[1]] <- out_data

#define the circle  ################
out_come<-c("abdominal aneurysm")
#####################
##储存结果
save.image("IV_SETUP.Rdata")
load("IV_SETUP.Rdata")


######################
results <- list()
het <- list()
pleio<- list()
single <- list()
res_single <- list()
library("ggplot2")

for (i in c(1:length(ex_pore))){
  for (j in c(1:length(out_come))){
    dat <- harmonise_data(
      exposure_dat = exp_dat[[i]],
      outcome_dat = out_dat[[j]],
      action = 2
    )
    
    res <- mr(dat)
    res$exposure=ex_pore[i]
    res$outcome=out_come[j]
    print(paste0("------", ex_pore[i], " & ",out_come[j],"------"))
    
    print(generate_odds_ratios(res))
    
    #primary results
    results[[length(out_come)*(i-1)+j]]    <- generate_odds_ratios(res)
    
    ###异质性检验
    het[[length(out_come)*(i-1)+j]] <- mr_heterogeneity(dat)
    #多效性检验
    pleio[[length(out_come)*(i-1)+j]] <- mr_pleiotropy_test(dat)
    #逐个剔除检验
    single[[length(out_come)*(i-1)+j]] <- mr_leaveoneout(dat)
    # write.csv(single,file = "mr_leaveoneout_single.csv")
    # pdf(paste0("leaveoneout",ex_pore[i],".pdf"),width = 12,height = 8)
    mr_leaveoneout_plot(single[[length(out_come)*(i-1)+j]])
    ggsave(paste0("leaveoneout",ex_pore[i],".pdf"),width = 12,height = 8)
    #dev.off()
    
    #散点图
    #pdf(paste0("scatter",ex_pore[i],".pdf"),width = 12,height = 8)
    mr_scatter_plot(res,dat)
    ggsave(paste0("scatter",ex_pore[i],".pdf"),width = 12,height = 8)
    #dev.off()
    #森林图
    res_single[[length(out_come)*(i-1)+j]] <- mr_singlesnp(dat)
    # write.csv(res_single,file = "res_single_forest_plot.csv")
    #pdf(paste0("forest",ex_pore[i],".pdf"),width = 12,height = 8)
    mr_forest_plot(res_single[[length(out_come)*(i-1)+j]])
    ggsave(paste0("forest",ex_pore[i],".pdf"),width = 12,height = 8)
    #dev.off()
  }
}
#################################################
results_allIV  <- do.call(rbind, results)
het_allIV <- do.call(rbind, het)
pleio_allIV<- do.call(rbind, pleio)
single_allIV <- do.call(rbind, single)
res_single_allIV <- do.call(rbind, res_single)




#################################################
#format(round(x, 2), nsmall = 2)
results_allIV$estimate <- paste0(format(round(results_allIV$or, 2), nsmall = 2), " (", 
                                 format(round(results_allIV$or_lci95, 2), nsmall = 2), "-",
                                 format(round(results_allIV$or_uci95, 2), nsmall = 2), ")")

row_x <- rownames (results_allIV[which(results_allIV$pval > 0.05), ])

results_allIV$pvalue          <- format(results_allIV$pval, scientific = TRUE, digits = 2)
results_allIV[row_x, ]$pvalue <- format(round(results_allIV[row_x, ]$pval, 2), nsmall = 2)


###################### output ###############
outname1="/results.csv"
write.table(results_allIV[,c(3:ncol(results_allIV))],path_res%+%outname1,sep = ",",row.names =FALSE,col.names =TRUE,quote =TRUE)
outname2="/het_allIV.csv"
write.table(het_allIV,path_res%+%outname2,sep = ",",row.names =FALSE,col.names =TRUE,quote =TRUE)
outname3="/pleio_allIV.csv"
write.table(pleio_allIV,path_res%+%outname3,sep = ",",row.names =FALSE,col.names =TRUE,quote =TRUE)
outname4="/single_allIV.csv"
write.table(single_allIV,path_res%+%outname4,sep = ",",row.names =FALSE,col.names =TRUE,quote =TRUE)
outname5="/res_single_allIV.csv"
write.table(res_single_allIV,path_res%+%outname5,sep = ",",row.names =FALSE,col.names =TRUE,quote =TRUE)
########################################################################


save.image("all_result.Rdata")
load("all_result.Rdata")


