#R3.5.3
  
#analyses of microbiome risk scores
len<-nrow(score)
aaa[]<-0
mydata[,aa]<-0
for (i in 1:len) {
  for (j in 1:1919) {
    qq<-quantile(mydata[,score[i,1]-420],0.05)
    if(mydata[j,score[i,1]-420]>qq){
      mydata[j,4]<-mydata[j,4]+1
      aaa[j,i]<-1
    }
  }
}
for (i in 1:nrow(score2)) {
  for (j in 1:1919) {
    qq<-quantile(mydata[,score[i,1]-420],0.05)
    if(mydata[j,score[i,1]-420]>qq){
      mydata[j,4]<-mydata[j,4]-1
      aaa[j,16]<--1
    }

  }
}
q1<-quantile(mydata[,4],0.5)
mydata[,4][mydata[,4] <q1]<-0
mydata[,4][mydata[,4] >=q1]<-1
res<-glm(dm_c~run+sex0+age_c+bmi_c+alc_c+energy+bristol_scale,data = mydata,family = binomial(link = "logit"))
summary(res)
 
# analyses of clustering
library(cluster)
library(clusterSim)
library(factoextra)
data2<-scale(data)
df=t(data2)
t<-pamk(df,krange=2:20,criterion="asw",usepam=T,ns=20,critout=TRUE,diss=F)
pam.res<-pam(df,the_optimum_number)

 #scripts for figures
fviz_cluster(pam.res,df,repel =T, pointsize =4,labelsize = 20, main = " ",show.clust.cent = F)+labs(title = "Guangzhou Nutrition and Health Study Cohort")+theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),plot.title=element_text(size=24))+theme(legend.text = element_text(size = 18))+theme(legend.title = element_text(size = 18))+ theme(axis.text.x = element_text(size = 18))+ theme(axis.text.y = element_text(size = 18))
fviz_cluster(pam.res,df,repel =T, pointsize =4,labelsize = 20, main = " ",show.clust.cent = F)+labs(title = "Replication Cohort")+theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),plot.title=element_text(size=24))+theme(legend.text = element_text(size = 18))+theme(legend.title = element_text(size = 18))+ theme(axis.text.x = element_text(size = 18))+ theme(axis.text.y = element_text(size = 18))
fviz_cluster(pam.res,df,repel =T, pointsize =4,labelsize = 20, main = " ",show.clust.cent = F)+labs(title = "Antibiotic-taking Group")+theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),plot.title=element_text(size=24))+theme(legend.text = element_text(size = 18))+theme(legend.title = element_text(size = 18))+ theme(axis.text.x = element_text(size = 18))+ theme(axis.text.y = element_text(size = 18))

#scripts for the input of heat maps

result[]<-NA
t<-names(mydata)
for (j in 2:72){
  for (i in 73:94){
    su<-cor.test(mydata[,j],mydata[,i],method = "spearman")
    result[i-72,j]<-as.numeric(su$estimate)
    result[i-72,j+72]<-su$p.value
    result[i-72,1]<-t[i]
    rm(su)
  }
}

data=mydata[,73:94]
data2<-scale(data)
df=t(data2)
res.hc = eclust(df, "hclust") 
path<-t(result[,1:71])
cluster<-data.frame(1:71)
cluster[]<-NA
for (i in 1:22) {
  a<-res.hc$order[i]
  cluster<-cbind(cluster,path[,a])
}
write.csv(cluster,"heatmapr.csv", row.names = F)
