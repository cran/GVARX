#dat=Data
GVECMest <- function (data,p=2,FLag, lag.max=NULL, type="const", ic,weight.matrix=NULL){
  ID<-NULL
  lag.max =lag.max
  DIFF=c(FALSE,TRUE)[2]
  p=p
  type=type

  idCol=which(colnames(data) == "ID")
  timeCol=which(colnames(data)=="Time")

  variates=colnames(data[,-c(idCol,timeCol)]) #names of column variables
  dat=as.data.frame(data[,-timeCol]) #Data with ID column only

  NAME=as.character(unique(dat[,"ID"]))
  N=length(NAME) #number of countries

  timeID=as.character(subset(data, ID==NAME[1])[,timeCol])
  Year=unique(as.character(lubridate::year(timeID)))

  output1=output2=output3=outRSD=list()
  outputP=outputP1=NULL

  Ft=GVAR_Ft(data,weight.matrix=weight.matrix)

for(i in 1:N) {

  if (DIFF==TRUE) { exo0=as.matrix(diff(Ft[[i]])) }   else {exo0=as.matrix(Ft[[i]])}
  exo=embed(exo0,FLag)   # Foreign variables
  ytmp=as.matrix(subset(data,ID==NAME[i])[,-c(1:2)])
  varnames=paste(NAME[i],paste(".",variates,sep=""),sep="")
  colnames(ytmp)=varnames
  y0=as.matrix(diff(ytmp)) #1-st order difference
  DIFF.no1=nrow(y0)-nrow(exo)
  yDIFF=y0[-c(1:DIFF.no1),]

  exoNAMES=NULL
  for (j in 1:FLag) {
    exoNAMES=c(exoNAMES,paste("F.",variates,".Lag",j-1,sep=""))
  }

  colnames(exo)=c(exoNAMES)

  #library()
  VECM.EG<-tsDyn::lineVar(ytmp, lag=p, model="VECM",estim=c("2OLS","ML")[2])
  beta0=summary(VECM.EG)$model.specific$beta
  ECT.r0=ytmp %*% beta0
  ECT.r1=as.matrix(ECT.r0)

  ECT.r1=as.matrix(ECT.r1[-length(ECT.r1),])

  DIFF.no2=length(ECT.r1)-nrow(yDIFF) #print(DIFF.no2)
  if (DIFF.no2==0) {ECT.r1=as.matrix(ECT.r1)} else {
  ECT.r1=as.matrix(ECT.r1[-c(1:DIFF.no2),])}

#  print(nrow(ECT.r1)); print(nrow(exo))
  EXO=as.matrix(cbind(exo,ECT.r1))
  colnames(EXO)=c(colnames(exo),"ECT")

# print(nrow(yDIFF));print(length(ECT.r1))

    ##==Compute Country-specific VAR results
  temp.out=.GVAR0(y=yDIFF, p = p, type = type, exogen = EXO, lag.max = lag.max, season = NULL,ic=ic)
  tmp <- temp.out
  for (q in 1:length(variates)) {
    output1[[i]]=tmp
    output2[[i]]=lmtest::coeftest(temp.out,vcov = vcovHC)
    output3[[i]]=temp.out$NWHAC
  }
  outputP=rbind(outputP,temp.out$p)

  ##===Estimate ordinary VAR and get residuals
 ECT.r1=as.matrix(ECT.r1)
 colnames(ECT.r1)="ECT"
 #print(ECT.r1)
rsdOUTtmp=.GVAR0(yDIFF, p=p, type=type, exogen = ECT.r1, season=NULL, lag.max=lag.max,ic=ic)
  outRSD[[i]]=rsdOUTtmp
  outputP1=rbind(outputP1,rsdOUTtmp$p)

} # end of 1 st i country-specific n loop



lagmatrix=data.frame(NAME,outputP)
colnames(lagmatrix)=c("Country","lags")

lagmatrix1=data.frame(NAME,outputP1)
colnames(lagmatrix1)=c("Country","lags")
# lagmatrix1

maxlag=max(outputP)
maxlag1=max(outputP1)

## 2-nd i loop begins
gvecmRSD=NULL;vecmRSD=NULL
for (i in 1:N)  {

  ##== Compute GVAR residuals
  thislag=outputP[i]
  rmobs=maxlag-thislag
  gvarrsd=resid(output1[[i]])
  if(rmobs==0) {gvarrsd=gvarrsd} else {gvarrsd=gvarrsd[-c(1:rmobs),]}
  gvecmRSD=cbind(gvecmRSD,gvarrsd)

  ##== Compute VAR residuals*
  thislag1=outputP1[i]
  rmobs1=maxlag1-thislag1
  varrsd=resid(outRSD[[i]])
  if(rmobs1==0) {varrsd=varrsd} else {varrsd=varrsd[-c(1:rmobs1),]}
  vecmRSD=cbind(vecmRSD,varrsd)

}  ## end of 2nd i country loop

results <- list(gvecm=output1,White=output2,NWHAC=output3,p=temp.out$p,K=temp.out$K,type=temp.out$type,datamat=data,lagmatrix =lagmatrix,lagmatrix1 =lagmatrix1, exoLag=FLag,Ft= Ft,NAMES=NAME,gvecmRSD=gvecmRSD,vecmRSD=vecmRSD)

return(results)

} ## end of function GVECMest()






getCOEF <- function(out,sheet) {
  save=FALSE
  filename=NULL

  k=sheet
  endo=out$K

  tmp0=coef(out$gvar[[k]])
  my0=data.frame(tmp0);
  coef.no=length(rownames(my0))
  newNAMES=c("Estimate","Std.Error","tstat","Pvalue")

  temp0=NULL
  for (j in 1:endo) {
    myexport0=my0[,((j-1)*4+1):(j*4)]
    colnames(myexport0)=newNAMES
    myexport=rbind(myexport0,newNAMES)
    temp0=rbind(temp0,myexport)
  }
  d=nrow(temp0)
  temp0=temp0[-d,]

  if (save==TRUE) {write.csv(temp0,file=paste(filename,".csv",sep="")) }
  return(coef=temp0)
}


getCOEFexo <- function(out) {
  save=FALSE
  filename=NULL

  endo=out$K
  N=length(out$NAMES)
  newNAMES=c("Estimate","Std.Error","tstat","Pvalue")
  coefEXO=NULL
  for (i in 1:N) {

    tmp0=coef(out$gvar[[i]])
    my0=data.frame(tmp0)
    coef.no=length(rownames(my0))
    exo.no=out$exoLag*endo
    a=coef.no-exo.no+1
    tempEXO=NULL
    for (j in 1:endo) {
      myexport0=my0[,((j-1)*4+1):(j*4)]
      colnames(myexport0)=newNAMES

      myexportF=rbind(myexport0[a:(a+endo-1),],newNAMES)
      tempEXO=rbind(tempEXO,myexportF)

    }

    coefEXO=rbind(coefEXO,i,tempEXO)

    if (save==TRUE) {write.csv(coefEXO,file=paste(filename,".csv",sep="")) }

  }
  return(coefEXO)
}



getWhiteCOEF <- function(out,sheet) {
  save=FALSE
  filename=NULL

  k=sheet
  endo=out$K
  temp1=out$White[[k]]
  my1=data.frame(temp1[,1:4])
  tmpRowNAMES=rownames(my1)
  newRowNAMES=NULL
  for (j in 1:(length(tmpRowNAMES)/out$K)) {
    newRowNAMES=c(newRowNAMES,strsplit(tmpRowNAMES,":")[[j]][2])
  }


  rhs.no=nrow(my1)/out$K
  newcolNAMES=c("Estimates","White Std.","tstat","Pvalue")
  colnames(my1)=newcolNAMES
  endo=out$K
  temp0=NULL
  for (j in 1:endo) {
    myexport0=my1[((j-1)*rhs.no+1):(j*rhs.no),]
    rownames(myexport0)=newRowNAMES
    myexport1=rbind(newcolNAMES,myexport0)

    temp0=rbind(temp0,myexport1)
  }
  colnames(temp0)=NULL
  if (save==TRUE) { write.csv(temp0,file=paste(filename,".csv",sep="")) }
  return(temp0)
}



getWhiteCOEFexo <- function(out) {
  save=FALSE
  filename=NULL

  N=length(out$NAMES)
  endo=out$K
  all=NULL
  for (i in 1:N) {
    temp1=out$White[[i]]
    my1=data.frame(temp1[,1:4])
    tmpRowNAMES=rownames(my1)
    newRowNAMES=NULL
    for (j in 1:(length(tmpRowNAMES)/out$K))
    {
      newRowNAMES=c(newRowNAMES,strsplit(tmpRowNAMES,":")[[j]][2])
    }

    coef.no=nrow(my1)/out$K
    newcolNAMES=c("Estimates","White Std.","tstat","Pvalue")
    colnames(my1)=newcolNAMES
    endo=out$K
    exo.no=out$exoLag*endo;
    a=coef.no-exo.no+1;
    temp0=NULL
    for (j in 1:endo)
    {
      myexport.0=my1[((j-1)*coef.no+1):(j*coef.no),]
      myexport0=myexport.0[a:(a+endo-1),]
      newRowNAMES1=newRowNAMES[a:(a+endo-1)]
      rownames(myexport0)=newRowNAMES1
      myexport1=rbind(myexport0,newcolNAMES)
      temp0=rbind(temp0,myexport1)
    }

    all=rbind(all,i,temp0)

  }

  #colnames(temp0)=NULL

  if (save==TRUE) { write.csv(all,file=paste(filename,".csv",sep="")) }

  return(all)
}

getNWCOEF <- function(out,sheet) {
  save=FALSE
  filename=NULL

  k=sheet
  endo=out$K
  newcolNAMES=c("Estimates","NW.std","tstat","pvalue")

  temp0=NULL
  for (j in 1:endo)
  {
    ACOV.tmp=out$NWHAC[[k]][j]
    ACOV=diag(as.matrix(as.data.frame(ACOV.tmp)))
    NW.std=sqrt(ACOV)
    ACOEF=as.data.frame(coef(out$gvar[[k]])[j])
    colnames(ACOEF)=NULL
    newrowNAMES=rownames(ACOEF)
    Pvalue=1-pnorm(abs(ACOEF[,1]/NW.std))
    NEW=data.frame(ACOEF[,1],NW.std,ACOEF[,1]/NW.std, Pvalue)
    colnames(NEW)=newcolNAMES
    rownames(NEW)=newrowNAMES
    New=rbind(newcolNAMES,NEW)
    temp0=rbind(temp0,New)
  }
  colnames(temp0)=NULL
  if (save==TRUE) { write.csv(temp0,file=paste(filename,".csv",sep="")) }
  return(temp0)
}



getNWCOEFexo <- function(out) {
  save=FALSE
  filename=NULL

  N=length(out$NAMES)
  endo=out$K
  newcolNAMES=c("Estimates","NW.std","tstat","pvalue")

  CoefEXO=NULL
  for (i in 1:N)
  {
    NW=out$NWHAC[[i]]
    tmp0=coef(out$gvar[[i]])
    my0=data.frame(tmp0)
    coef.no=length(rownames(my0))
    exo.no=out$exoLag*endo
    a=coef.no-exo.no+1

    tempEXO=NULL
    for (j in 1:endo)
    {
      ACOV.tmp=NW[j]
      ACOV=diag(as.matrix(as.data.frame(ACOV.tmp)))
      NW.std=sqrt(ACOV)
      ACOEF=as.data.frame(tmp0[j])
      #    colnames(ACOEF)=NULL
      newrowNAMES=rownames(ACOEF)
      Pvalue=1-pnorm(abs(ACOEF[,1]/NW.std))
      NEW=data.frame(ACOEF[,1],NW.std,ACOEF[,1]/NW.std, Pvalue)
      colnames(NEW)=newcolNAMES
      rownames(NEW)=newrowNAMES
      New=rbind(newcolNAMES,NEW)
      New=New[(a+1):(a+endo),]
      tempEXO=rbind(tempEXO,New,newcolNAMES)
    }
    CoefEXO=rbind(CoefEXO,i,tempEXO)

  }

  #colnames(CoefEXO)=NULL
  if (save==TRUE) { write.csv(CoefEXO,file=paste(filename,".csv",sep="")) }

  return(CoefEXO)
}


averageCORgvecm <-function (out) {
  var.no=out$K
  varRSD=out$vecmRSD
  V1=t(matrix(colnames(varRSD),var.no,))

  gvarRSD=out$gvecmRSD
  V=t(matrix(colnames(gvarRSD),var.no,))

  varRSD.cor=list();gvarRSD.cor=list()
  for (j in 1:var.no) {
    temp1=(apply(cor(varRSD[,V1[,j]]),2,sum)-1)/(length(V1[,j])-1)
    varRSD.cor[[j]]=temp1

    temp=(apply(cor(gvarRSD[,V[,j]]),2,sum)-1)/(length(V[,j])-1)
    gvarRSD.cor[[j]]=temp
  }
  results=list(vecmRSDcor=varRSD.cor,gvecmRSDcor=gvarRSD.cor)
  return(results)

}
