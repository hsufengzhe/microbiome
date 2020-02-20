#R 3.5.3
library(TwoSampleMR)
library(MRPRESSO)

result<-data.frame()

for (i in 1:5 ){
 ab<-paste0(i,".txt")
  a<-paste0(i,sep="_micro_",ab)
  c<-read.csv(a)
  nn<-62
  if (nrow(c)>0){
    ii<-i
    for (m in 1:nn){
      ab<-paste0(m,".txt")
      a<-paste0(i,sep="_micro_",ab)
      exposure<-read_exposure_data(a,sep=" ")
      b<-paste0(m,".txt")
      bb<-paste0(i,sep="_outcome_",b)
      row<-(i-1)*nn+m
      result[row,1]<- i
      cc<-read.csv(bb)
      if (nrow(cc)>0 && nrow(exposure)<3 ){
        outcome<-read_outcome_data(bb,sep=" ")
        dat <- harmonise_data(exposure, outcome, action=2)
        res <-  mr(dat)
        egger<- mr_egger_regression(dat$beta.exposure,dat$beta.outcome,dat$se.exposure ,dat$se.outcome,T)
        for (j in 1:1) {
          aa<-paste0(round(res[j,7],3),sep="(",round(res[j,8],3))
          aa<-paste0(aa,")",sep="")
          result[row,6]<-aa
          result[row,7]<-round(res[j,9],3)
        }
        result[row,15]<-egger$pval_i
        result[row,16]<-egger$b_i
        result[row,17]<-egger$se_i
        result[row,18]<-IVW$Q_pval
      }
      else if(nrow(cc)>0 && nrow(exposure)>3 ){
        outcome<-read_outcome_data(bb,sep=" ")
        dat <- harmonise_data(exposure, outcome, action=2)
        res <-  mr(dat)
        presso<-mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure ="se.exposure", data= dat,OUTLIERtest = TRUE, DISTORTIONtest = TRUE,NbDistribution = 1000,  SignifThreshold = 0.05)
        egger<- mr_egger_regression(dat$beta.exposure,dat$beta.outcome,dat$se.exposure ,dat$se.outcome,T)
        res_single <- mr_singlesnp(dat, all_method=c("mr_ivw", "mr_two_sample_ml"))
        write.table(res_single,paste0(bb,sep="_","res_single"),row.names = F,sep=" ")
        mr_forest_plot(res_single)
        png(file=paste(bb,".png"), bg="transparent")
        dev.off()
        IVW <- mr_ivw(dat$beta.exposure,dat$beta.outcome,dat$se.exposure ,dat$se.outcome)
        for (j in 1:5) {
          aa<-paste0(round(res[j,7],3),sep="(",round(res[j,8],3))
          aa<-paste0(aa,")",sep="")
          result[row,j*2]<-aa
          result[row,j*2+1]<-round(res[j,9],3)
        }
        ccc<-paste0(round(presso$`Main MR results`$`Causal Estimate`[1],3),sep="(",round(presso$`Main MR results`$Sd[1],3))
        result[row,12]<-paste0(ccc,")",sep="")
        result[row,13]<-presso$`Main MR results`$`P-value`[1]
        result[row,14]<-presso$`MR-PRESSO results`$`Global Test`$Pvalue
        result[row,15]<-egger$pval_i
        result[row,16]<-egger$b_i
        result[row,17]<-egger$se_i
        result[row,18]<-IVW$Q_pval 
      }
      else {
        result[((ii-1)*nn+1):(ii*nn),1]<- ii
      }
    }
  }
  else {
    result[((i-1)*nn+1):(i*nn),1]<- ii
  }
}
