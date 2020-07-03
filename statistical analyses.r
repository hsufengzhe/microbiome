#R3.5.3
  
#analyses of microbiome risk scores
len<-nrow(score)
aaa[]<-0
mydata$taxon<-0
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
res<-glm(dm_c~taxon+sex0+age_c+bmi_c+alc_c+energy+bristol_scale,data = mydata,family = binomial(link = "logit"))
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

#scripts for Jaccard similarity coefficient
data=read.csv("discovery_genus_pathway.csv", header=T, row.names=1)
df=t(data)
pam.res<-pam(df,5)
c1<-pam.res$clustering

data=read.csv("replication_genus_pathway.csv", header=T, row.names=1)
df=t(data)
pam.res<-pam(df,5)
c2<-pam.res$clustering

data=read.csv("antibiotic.csv", header=T, row.names=1)
df=t(data)
pam.res<-pam(df,5)
c3<-pam.res$clustering

c11<-as.numeric(c1)
c22<-as.numeric(c2)
c33<-as.numeric(c3)

m<-1
result1<-data.frame()
for (i in 1:21) {
  for (j in (i+1):22) {
    result1[m,1]<-i
    result1[m,2]<-j
    if(c11[i]==c11[j])
      result1[m,3]<-1
    else result1[m,3]<-0
    m<-m+1
  }
  
}

m<-1
result2<-data.frame()
for (i in 1:21) {
  for (j in (i+1):22) {
    result2[m,1]<-i
    result2[m,2]<-j
    if(c22[i]==c22[j])
      result2[m,3]<-1
    else result2[m,3]<-0
    m<-m+1
  }
  
}

m<-1
result3<-data.frame()
for (i in 1:21) {
  for (j in (i+1):22) {
    result3[m,1]<-i
    result3[m,2]<-j
    if(c33[i]==c33[j])
      result3[m,3]<-1
    else result3[m,3]<-0
    m<-m+1
  }
  
}
ss<-0
sd<-0
ds<-0
dd<-0
for (i in 1:231) {
 if(result1[i,3]==1 && result2[i,3]==1) ss<-ss+1
 if(result1[i,3]==1 && result2[i,3]==0) sd<-sd+1
 if(result1[i,3]==0 && result2[i,3]==1) ds<-ds+1
 if(result1[i,3]==0 && result2[i,3]==0) dd<-dd+1
}
ss+sd+ds+dd
jaccard_fh=ss/(ss+sd+ds)

#scripts for spieceasi
metabolites<-data.matrix(round(data[,c(3022,3026,3032,3033,3061,3079,3084,3087,3135,3138,3139,3141,3188,3189,3216,3219,3230,3261,3276,3278,3279,3287)]*100000))
se.hmp2 <- spiec.easi(metabolites, method='mb', nlambda=1000,
                      lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))

dtype <- c(rep(1,22), rep(2,end-start+1))
plot(adj2igraph(getRefit(se.hmp2)), vertex.color=dtype+1, vertex.size=10,vertex.attr = list( cc[1:22,1]))
write.csv(data.matrix( se.hmp2[["refit"]][["stars"]]),"speci.csv",quote = F)

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
