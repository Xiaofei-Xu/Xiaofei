library(mvtnorm)
library(foreach)
library(doParallel)
library(Matrix)
library(numDeriv)
library(stringr)
library(numDeriv)
# run first 3 section data 43-45 session
X<-readRDS("Xs13.rds")  
Y<-readRDS("Cs13.rds")
Ytemp=Y[,c(1,3,4,10)]# Only use 4 bigrams:  ``american people'', ``member congress'',  ``men women'' and ``year ago''
bigramName=str_replace_all(colnames(Ytemp), fixed(" "), "")
Mu <- rowSums(Ytemp) # compute m_t

Xt=X[,c(1,3,4,5)]# Only use 4 covariates:republican, gender, chamber, majority
bigramid=1 # run code for one bigram, change the index to run code for other bigram counts
Yt=Ytemp[,bigram] 
nT=length(Yt)
z_r<-sum(Yt==0)/length(Yt) # initial value of zero-inflation paramater \rho

no_cores<-detectCores()
cl<-makeCluster(no_cores)
registerDoParallel(cl)

weidth=350 # rolling window size
SampleSize = nT
TestWindowStart      = weidth+1
TestWindow           = TestWindowStart:SampleSize
testsample           = length(TestWindow)

start_time <- Sys.time()

file=paste(bigramName[bigramid],".csv", sep="")
save<-foreach(k = 1:testsample, .combine='rbind',.multicombine = TRUE,.packages='mvtnorm') %dopar%
  {
    nob=TestWindowStart+k-1
    trainperiod=(nob-weidth+1):nob
    yt=Yt[trainperiod]
    xt=Xt[trainperiod,]
    mu=Mu[trainperiod]
    loglik = function(star){
      lik2=0
      lam_n	=	NULL
      lamr_n	=	NULL
      pcov=dim(xt)[2]
      #star=c(ome,beta,rho,phi,g1,g2,g3,g4)
      lam_n[1]=yt[1]
      lamr_n[1]=mu[1]*(1-star[4])/(1-star[3])*lam_n[1]
      for (t in 2:weidth)
      { 
        lam_n[t] 	= exp(star[1] + star[2]*lam_n[t-1]+star[5:8]%*%t(xt[t,]))
        lamr_n[t]   =mu[t]*(1-star[4])/(1-star[3])*lam_n[t]
        if (yt[t]==0) 
        {lik2= lik2+log(star[3]+(1-star[3])*exp(-lamr_n[t]))}
        else         {lik2= lik2+log(1-star[3])+log(lamr_n[t])+(yt[t]-1)*log(lamr_n[t]+star[4]*yt[t])-(lamr_n[t]+star[4]*yt[t])}
        
      }
      neglogliki=-lik2
      return(neglogliki)
    }
    # 0<=rho, phi<1, -1<beta<1
    star0= c(-0.2,0.1,z_r,0.4,rep(0.01,4))
    
    ui=matrix(c(0,0,1,rep(0,5),0,0,0,1,rep(0,4), 0,0,-1,rep(0,5),
                0,0,0,-1,rep(0,4)),nrow=4,byrow=T)
    ci=matrix(c(0,0,-0.9999,-0.9999),nrow=4)
 
    fit=constrOptim(star0, loglik, NULL,  ui=ui, ci=ci,control = list(maxit = 5000))
   
    Flam_t=NULL
    Flamr_n=0
    Flam_t[1]=yt[1]
    
    for (t in 2:weidth)
    { 
      Flam_t[t] 	= exp(theta[1] + theta[2]*Flam_t[t-1]+theta[5:8]%*%t(xt[t,]))
    }
    Flamr_n  =mu[t]*(1-theta[4])/(1-theta[3])*Flam_t[weidth]
    Flamr_mean=mu[t]*Flam_t[weidth]
    resultvec=c(theta, Flam_t[weidth],Flamr_n,Flamr_mean)
    
    df = data.frame(t(resultvec))
    rownames(df) = k
    write.table(df,file,sep=',',col.names=FALSE,append = TRUE)
    df
  }

end_time <- Sys.time()
end_time - start_time

#finalreult=save
#saveRDS(finalreult, "congress34.rds")

stopCluster(cl)



