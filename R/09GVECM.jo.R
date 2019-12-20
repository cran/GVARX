GVECM.jo <- function (data,p=2,FLag=3,ecdet = "const", type = "eigen",spec = "longrun", season = NULL,weight.matrix)
    {
  ID<-NULL
DIFF=TRUE
NAME=unique(as.character(data[,"ID"]))
endo.no=ncol(data)-2
N=length(NAME)
#weight.matrix=weight.matrix

##==Create country-specific new dataset

Ft=GVAR_Ft(data,weight.matrix)

JO.test=VECMoutputs=RESID=list()

for(jj in 1:N) {
  datai=Ft[[jj]]

  endcol=ncol(datai)
  ytemp=datai
  exo=embed(datai,FLag)

  if (DIFF==TRUE) { exo0=as.matrix(diff(datai)) }   else {exo0=as.matrix(datai)}

  newEXOname.tmp=colnames(ytemp)
  newEXOname=paste(newEXOname.tmp,"Exo",sep=".")

  Dexo=embed(as.matrix(exo0),FLag)

  embed_names=NULL
  for (j in 1:FLag){
    embed_names_tmp=paste(newEXOname,".L",j-1,sep="")
    embed_names=c(embed_names,embed_names_tmp)
  }
 colnames(Dexo)=embed_names

  y1=subset(data[,-2],ID==NAME[jj])[,-1]
  y1=y1[-(1:(nrow(y1)-nrow(Dexo))),]

#  p=p+1   # lags for endogenous regressors
#  Plags=vars::VARselect(y1, lag.max = 3, type = c("const", "trend", "both", "none")[1])$selection
#  p=Plags[1]
  if(p==1) {p=p+1}

  out_vecm1 <- urca::ca.jo(y1, K=p, ecdet = ecdet, type = type,spec = spec, season = NULL,dumvar =Dexo)
   JO.test[[jj]]<-urca::summary(out_vecm1)
   VECMoutputs[[jj]]<-out_vecm1

     obj1=vars::vec2var(out_vecm1, r = (endo.no-1))
    RSD0=resid(obj1)
   colnames(RSD0)=paste0("RSD.",colnames(y1))
     RESID[[jj]]=RSD0

  } # jj loop ends
results <-list(JO.test=JO.test,RESID=RESID,VECMoutputs=VECMoutputs)


}



