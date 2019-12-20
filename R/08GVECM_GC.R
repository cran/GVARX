.grangerGVECM <- function (data,p=2,FLag=3, type="const",lag.max=NULL, ic, weight.matrix=NULL, boot.no=0) {
ID <- NULL
  stopifnot(ncol(data)==4)
ic=ic
 p=p
 type=type
DIFF=c(FALSE,TRUE)[2]
 idCol=which(colnames(data) == "ID")
 timeCol=which(colnames(data)=="Time")
 variates=colnames(data[,-c(idCol,timeCol)]) #names of column variables
 dat=data[,-timeCol] #Data with ID column only
 NAME=as.character(unique(dat[,"ID"]) )
 N=length(NAME) #number of countries

 timeID=as.character(subset(data, ID==NAME[1])[,timeCol])
 Year=unique(as.character(lubridate::year(timeID)))

output1=output2=output3=outRSD=EXO=list()
Y1gcY2=Y2gcY1=NULL
Y1gcY2.GVAR=Y2gcY1.GVAR=NULL

Ft=GVAR_Ft(data,weight.matrix)

for(i in 1:N) {

  if (DIFF==TRUE) { exo0=as.matrix(diff(Ft[[i]])) }   else {exo0=as.matrix(Ft[[i]])}
  exo=embed(exo0,FLag)   # Foreign variables
  ytmp=as.matrix(subset(data,ID==NAME[i])[,-c(1:2)])
  varnames=paste(NAME[i],paste(".",variates,sep=""),sep="")
  colnames(ytmp)=varnames
  y0=as.matrix(diff(ytmp)) #1-st order difference
  DIFF.no1=nrow(y0)-nrow(exo)
  yDIFF=y0[-c(1:DIFF.no1),] #

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
  DIFF.no2=length(ECT.r1)-nrow(yDIFF)#;  print(DIFF.no2)
  ECT.r1=ECT.r1[-c(1:DIFF.no2),]
  EXO=as.matrix(cbind(exo,ECT.r1))
  colnames(EXO)=c(colnames(exo),"ECT")

##==Compute Granger Causality based upon ordinary VAR
  ECT.r1=as.matrix(ECT.r1)
  colnames(ECT.r1)="ECT"
  lag.max=lag.max
VAR.GC=vars::VAR(y=yDIFF, p = p, type = type, lag.max = lag.max, ic=ic, exogen = ECT.r1)
if (boot.no==0) {
GCY1=vars::causality(VAR.GC,cause=varnames[1],vcov.=vcovHC(VAR.GC))$Granger
GCY2=vars::causality(VAR.GC,cause=varnames[2],vcov.=vcovHC(VAR.GC))$Granger } else {
GCY1=vars::causality(VAR.GC,cause=varnames[1],boot = TRUE, boot.runs = boot.no)$Granger
GCY2=vars::causality(VAR.GC,cause=varnames[2],boot = TRUE, boot.runs = boot.no)$Granger}

Y1gcY2_tmp=cbind(as.numeric(GCY1$statistic[1]),as.numeric(GCY1$p.value))
Y2gcY1_tmp=cbind(as.numeric(GCY2$statistic[1]),as.numeric(GCY2$p.value))
Y1gcY2=rbind(Y1gcY2,Y1gcY2_tmp)
Y2gcY1=rbind(Y2gcY1,Y2gcY1_tmp)

##==Compute Granger Causality based upon GVARX
GVAR.GC=vars::VAR(y=yDIFF, p = p, type = type, exogen = EXO, lag.max = lag.max, ic=ic)

if (boot.no==0) {GCY1.GVAR=vars::causality(GVAR.GC,cause=varnames[1],vcov.=vcovHC(GVAR.GC))$Granger
GCY2.GVAR=vars::causality(GVAR.GC,cause=varnames[2],vcov.=vcovHC(GVAR.GC))$Granger
} else {
GCY1.GVAR=vars::causality(GVAR.GC,cause=varnames[1], boot = TRUE, boot.runs = boot.no)$Granger
GCY2.GVAR=vars::causality(GVAR.GC,cause=varnames[2], boot = TRUE, boot.runs = boot.no)$Granger
}
Y1gcY2.GVAR_tmp=cbind(as.numeric(GCY1.GVAR$statistic[1]),as.numeric(GCY1.GVAR$p.value))
Y2gcY1.GVAR_tmp=cbind(as.numeric(GCY2.GVAR$statistic[1]),as.numeric(GCY2.GVAR$p.value))
Y1gcY2.GVAR=rbind(Y1gcY2.GVAR,Y1gcY2.GVAR_tmp)
Y2gcY1.GVAR=rbind(Y2gcY1.GVAR,Y2gcY1.GVAR_tmp)

} # end of 1 st jj countryloop

options(digits=4)
## Reporting Granger Causality output based on VECM
rownames(Y1gcY2)=NAME
null1=paste(variates[1],"does not Granger cause",variates[2])
colnames(Y1gcY2)=c(null1,"VECM, pvalue")

rownames(Y2gcY1)=NAME
null2=paste(variates[2],"does not Granger cause",variates[1])
colnames(Y2gcY1)=c(null2,"VECM, pvalue")


## Reporting Granger Causality output based on GVECM
rownames(Y1gcY2.GVAR)=NAME
null1=paste(variates[1],"does not Granger cause",variates[2])
colnames(Y1gcY2.GVAR)=c(null1,"GVECM, pvalue")

rownames(Y2gcY1.GVAR)=NAME
null2=paste(variates[2],"does not Granger cause",variates[1])
colnames(Y2gcY1.GVAR)=c(null2,"GVECM, pvalue")

#print(Y1gcY2);print(Y2gcY1)
#print(Y1gcY2.GVAR);print(Y2gcY1.GVAR)

results <- list(y1GCy2.vecm=Y1gcY2,y2GCy1.vecm=Y2gcY1,y1GCy2.gvecm=Y1gcY2.GVAR,y2GCy1.gvecm=Y2gcY1.GVAR)

return(results)

} ## end of function GC()
