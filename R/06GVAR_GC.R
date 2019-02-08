grangerGVAR <- function (data,p=2,FLag=2, type="const",lag.max, ic, weight.matrix=NULL) {

  stopifnot(ncol(data)==4)

 p=p
 type=type
 weight=weight.matrix

 ID="ID"
 Time="Time"
 idCol=which(colnames(data) == ID)
 timeCol=which(colnames(data)==Time)
 variates=colnames(data[,-c(idCol,timeCol)]) #names of column variables
 dat=data[,-timeCol] #Data with ID column only
 NAME=as.character(unique(dat[,"ID"]) )
 N=length(NAME) #number of countries

 timeID=as.character(subset(data, ID==NAME[1])[,timeCol])
 Year=unique(as.character(lubridate::year(timeID)))

output1=list()
output2=list()
output3=list()
outRSD=list()
EXO=list()
Y1gcY2=NULL
Y2gcY1=NULL
Y1gcY2.GVAR=NULL
Y2gcY1.GVAR=NULL

for (i in 1:N)  {
  h=NAME[i]  # index number of home country

  exo.tmp=NULL  # Compute Exogenous Foreign Variables


  for (j in 2:(length(variates)+1)) {

    if (is.null(weight)){

      F.tmp=apply(matrix(dat[,j], ,N),1,mean)

    }
    if (isTRUE(is.matrix(weight.matrix))) {

      varMatrix=matrix(dat[,j], ,N)
      F.tmp=varMatrix %*% as.matrix(weight[,i])


    } else {

      dat_matrix=matrix(dat[,j], ,N)
      varMatrix=as.xts(dat_matrix,as.Date(timeID))

      F.tmp=NULL
      for (k in 1:length(Year)) {
        F.tmp0=varMatrix[Year[k]] %*% as.matrix(weight[[k]][,i])
        F.tmp=rbind(F.tmp,F.tmp0)
      }

    }

    exo.tmp = cbind(exo.tmp,F.tmp)

 }


 exo=embed(exo.tmp,FLag)   # Foreign variables
 ytmp=subset(dat,ID==h)[,-1];
 varnames=paste(h,paste(".",variates,sep=""),sep="")
colnames(ytmp)=varnames
y=ytmp[-c(1:(FLag-1)),]

exoNAMES=NULL
for (j in 1:FLag) {
exoNAMES=c(exoNAMES,paste("F.",variates,".Lag",j-1,sep=""))
}

colnames(exo)=c(exoNAMES)

##==Compute Granger Causality based upon ordinary VAR

VAR.GC=vars::VAR(y=y, p = p, type = type, lag.max = lag.max, ic=ic)
GCY1=vars::causality(VAR.GC,cause=varnames[1],vcov.=vcovHC(VAR.GC))$Granger
GCY2=vars::causality(VAR.GC,cause=varnames[2],vcov.=vcovHC(VAR.GC))$Granger

Y1gcY2_tmp=cbind(as.numeric(GCY1$statistic[1]),as.numeric(GCY1$p.value))
Y2gcY1_tmp=cbind(as.numeric(GCY2$statistic[1]),as.numeric(GCY2$p.value))
Y1gcY2=rbind(Y1gcY2,Y1gcY2_tmp)
Y2gcY1=rbind(Y2gcY1,Y2gcY1_tmp)

##==Compute Granger Causality based upon GVARX
GVAR.GC=vars::VAR(y=y, p = p, type = type, exogen = exo, lag.max = lag.max, ic=ic)

GCY1.GVAR=vars::causality(GVAR.GC,cause=varnames[1],vcov.=vcovHC(GVAR.GC))$Granger
GCY2.GVAR=vars::causality(GVAR.GC,cause=varnames[2],vcov.=vcovHC(GVAR.GC))$Granger

Y1gcY2.GVAR_tmp=cbind(as.numeric(GCY1.GVAR$statistic[1]),as.numeric(GCY1.GVAR$p.value))
Y2gcY1.GVAR_tmp=cbind(as.numeric(GCY2.GVAR$statistic[1]),as.numeric(GCY2.GVAR$p.value))
Y1gcY2.GVAR=rbind(Y1gcY2.GVAR,Y1gcY2.GVAR_tmp)
Y2gcY1.GVAR=rbind(Y2gcY1.GVAR,Y2gcY1.GVAR_tmp)

} # end of 1 st jj countryloop

options(digits=4)
## Reporting Granger Causality output based on VAR
rownames(Y1gcY2)=NAME
null1=paste(variates[1],"does not Granger cause",variates[2])
colnames(Y1gcY2)=c(null1,"pvalue, VAR")

rownames(Y2gcY1)=NAME
null2=paste(variates[2],"does not Granger cause",variates[1])
colnames(Y2gcY1)=c(null2,"pvalue, VAR")


## Reporting Granger Causality output based on GVAR
rownames(Y1gcY2.GVAR)=NAME
null1=paste(variates[1],"does not Granger cause",variates[2])
colnames(Y1gcY2.GVAR)=c(null1,"pvalue, GVAR")

rownames(Y2gcY1.GVAR)=NAME
null2=paste(variates[2],"does not Granger cause",variates[1])
colnames(Y2gcY1.GVAR)=c(null2,"pvalue, GVAR")

#print(Y1gcY2);print(Y2gcY1)
#print(Y1gcY2.GVAR);print(Y2gcY1.GVAR)

results <- list(y1GCy2.var=Y1gcY2,y2GCy1.var=Y2gcY1,y1GCy2.gvar=Y1gcY2.GVAR,y2GCy1.gvar=Y2gcY1.GVAR)

return(results)

} ## end of function GC()
