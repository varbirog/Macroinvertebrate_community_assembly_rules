load("phyto.traits.RData")

load("phyto.plots.RData")
load("phyto.env.RData")
trait= read.table(file= "clipboard", sep="\t", header=T,encoding=c("windows-1250"),row.names = 1)
plots= read.table(file= "clipboard", sep="\t", header=T,encoding=c("windows-1250"),row.names = 1)
x= read.table(file= "clipboard", sep="\t", header=T,encoding=c("windows-1250"))

# function calculating Rao's quadratic entropy
f_rao<-function(x,d=dist.mat.Flagellated) {
  x=as.matrix(x)
  1-(x/sum(x))*d*t(x/sum(x))
}
# function calculating Rao's quadratic entropy simplyfied
f<-function(x,d) 1/(1-(x/sum(x))%*%d%*%t(x/sum(x)))

n<-nrow(plots)
nr<-9999 #999 #999 
#dist.mat
t_names=names(trait)
nt<-length(t_names)

for (i in t_names){
  
  assign(paste('dist.mat.', i, sep=''),as.matrix(dist(trait[,i])))
  assign(paste('RaoQ.', i, sep=''),rep(NA,n))
}
# ez a jo
for (k in 1:nt){
  nam=paste('dist.mat.', t_names[k], sep='')
  adat=eval(parse(text=nam))
  nam2 <- paste("RaoQ.", t_names[k], sep = "")
  adat2=eval(parse(text=nam2))

for (i in 1:n){
	a=as.matrix(plots[i,])
	adat2[i] = 1/(1-(a/sum(a))%*%adat%*%t(a/sum(a)))
}
  
  assign(paste('RaoQ.', t_names[k], sep=''),adat2)
  
}



for (i in t_names){
  
  assign(paste('random.RaoQ.', i, sep=''),matrix(NA,nrow=nr,ncol=n))
  
}

#Flagellated
	start.time <- Sys.time()
	for (k in t_names){
	  print(k)
	  
	  for (j in 1:nr){
	    print(j)
	    
	    nam <- paste("dist.mat.", k, sep = "")  
	    nam2 <- paste("random.RaoQ.", k, sep = "")
	    adat=eval(parse(text=nam)) 
	    adat2=eval(parse(text=nam2))
	    W<-sample(1:nrow(adat))
	    random.dist.mat<-adat[W,W]
	    #random.dist.mat<-dist.mat.aff[W,W]
	    for (i in 1:n){
	      adat2[j,i]<-f(as.matrix(plots[i,]),random.dist.mat)
	      
	    }
	    assign(nam2,adat2)
	  }
	}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


for (j in 1:nr){
  print(j)
  
  nam <- paste("dist.mat.", k, sep = "")  
  nam2 <- paste("random.RaoQ.", k, sep = "")
  adat=eval(parse(text=nam)) 
  adat2=eval(parse(text=nam2))
  W<-sample(1:nrow(adat))
  random.dist.mat<-adat[W,W]
  #random.dist.mat<-dist.mat.aff[W,W]
  for (i in 1:n){
    adat2[j,i]<-f(as.matrix(plots[i,]),random.dist.mat)
  
  }
  assign(nam2,adat2)
}



for (i in t_names){
  
  assign(paste('pvalue.', i, sep=''),-rep(0,n))
  
}




for (k in t_names){
  
  nam <- paste("RaoQ.", k, sep = "")  
  nam2 <- paste("random.RaoQ.", k, sep = "")  
  nam3<- paste("pvalue.", k, sep = "")
  
  adat1=eval(parse(text=nam))
  adat2=eval(parse(text=nam2))
  
  adat3=eval(parse(text=nam3))
  
  
  
for (i in 1:n){
  
    adat3[i]<-(sum(adat1[i]>=adat2[,i])+1)/(nr+1)
    
  }  
  
  adat3= qnorm(adat3)
  adat3[adat3==Inf]<-qnorm(0.999999)
  adat3[adat3==-Inf]<-qnorm(0.00001)
  assign(nam3,adat3)
  
  
  
  }





kimenet=rbind(
pvalue.gra,
pvalue.min,
pvalue.xyl,
pvalue.shr,
pvalue.gat,
pvalue.aff,
pvalue.pff,
pvalue.pre,
pvalue.par,
pvalue.oth,
pvalue.s1,
pvalue.s2,
pvalue.s3,
pvalue.s4,
pvalue.s5,
pvalue.s6,
pvalue.s7,
pvalue.cd1,
pvalue.cd2,
pvalue.cy1,
pvalue.cy2,
pvalue.cy3,
pvalue.dis1,
pvalue.dis2,
pvalue.dis3,
pvalue.dis4,
pvalue.life1,
pvalue.life2,
pvalue.life3,
pvalue.life4,
pvalue.fwl1,
pvalue.fwl2,
pvalue.fwl3,
pvalue.fwl4,
pvalue.fwl5,
pvalue.fwl6,
pvalue.fwl7,
pvalue.fwl8,
pvalue.wnb1,
pvalue.wnb2,
pvalue.wnb3,
pvalue.wnb4,
pvalue.wnb5,
pvalue.egg1,
pvalue.egg2,
pvalue.egg3,
pvalue.egg4,
pvalue.drift1,
pvalue.drift2,
pvalue.drift3)

# END STEPS 
save.image("data_9999.RData")

 write.table(kimenet, file="clipboard-16384",sep="\t",row.names=T,col.names=T)


dat= read.table(file= "kimenet_with_states.csv", sep="\t", header=T,encoding=c("windows-1250"),row.names = 1)


# Generating Plots

for(i in 4:ncol(dat)){
  
  y= dat %>%  select(sate,satus,colnames(dat)[i])

  ggplot(data=y,aes(x=satus,y=y[,3]))+
    geom_jitter(width = 0.10,alpha=0.5)+
    geom_boxplot(width=0.4,alpha=0.6)+
    geom_hline(yintercept =0)+
    ylim(-3, 3)+
    theme(plot.title = element_text(hjust = 0.5),axis.ticks.x=element_blank())+
    #annotate("text",x=1.5,y=0.8,label = paste("trait divergence\n  biotic interactions"))+
    #annotate("text",x=1.5,y=-0.8,label = paste("trait convergence\n environmental filtering"))+
    labs(y="Effect size",title=colnames(dat[i]),x="Group")
    
  
  ggsave(paste(colnames(dat)[i],".jpg"), device = "jpg")
  
}

