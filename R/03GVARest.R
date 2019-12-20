GVARest <- function (data,p=2,FLag, lag.max=NULL, type="const", ic,weight.matrix=NULL){
  ID<-NULL
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

for (i in 1:N) {
  exo=embed(Ft[[i]],FLag)   # Foreign variables
  ytmp=subset(data,ID==NAME[i])[,-c(1:2)];
  varnames=paste(NAME[i],paste(".",variates,sep=""),sep="")
  colnames(ytmp)=varnames
  y=ytmp[-c(1:(FLag-1)),]

  exoNAMES=NULL
  for (j in 1:FLag) {
  exoNAMES=c(exoNAMES,paste("F.",variates,".Lag",j-1,sep=""))
  }

colnames(exo)=c(exoNAMES)


##==Compute Country-specific VAR results
lag.max =lag.max
temp.out=.GVAR0(y=y, p = p, type = type, exogen = exo, lag.max = lag.max, season = NULL,ic=ic)
   tmp <- temp.out
   for (q in 1:length(variates)) {
   output1[[i]]=tmp
   output2[[i]]=lmtest::coeftest(temp.out,vcov = vcovHC)
   output3[[i]]=temp.out$NWHAC
   }
   outputP=rbind(outputP,temp.out$p)

##===Estimate ordinary VAR and get residuals
rsdOUTtmp=.GVAR0(y, p=1, type=type, season=NULL, lag.max=lag.max,ic=ic)
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
gvarRSD=NULL;varRSD=NULL
for (i in 1:N)  {

##== Compute GVAR residuals
thislag=outputP[i]
rmobs=maxlag-thislag
gvarrsd=resid(output1[[i]])
if(rmobs==0) {gvarrsd=gvarrsd} else {gvarrsd=gvarrsd[-c(1:rmobs),]}
gvarRSD=cbind(gvarRSD,gvarrsd)

##== Compute VAR residuals*
thislag1=outputP1[i]
rmobs1=maxlag1-thislag1
varrsd=resid(outRSD[[i]])
if(rmobs1==0) {varrsd=varrsd} else {varrsd=varrsd[-c(1:rmobs1),]}
varRSD=cbind(varRSD,varrsd)

}  ## end of 2nd i country loop

results <- list(gvar=output1,White=output2,NWHAC=output3,p=temp.out$p,K=temp.out$K,type=temp.out$type,datamat=data,lagmatrix =lagmatrix,lagmatrix1 =lagmatrix1, exoLag=FLag,Ft= Ft,NAMES=NAME,gvarRSD=gvarRSD,varRSD=varRSD)

return(results)

} ## end of function GVARest()






#.GVARirf <- function (Obj, impulse, response,ortho,cumulative,boot,ci = 0.95,n.ahead=15, Weights,shockNeg=FALSE,STD=FALSE,G=FALSE) {

#DATA=GVAR_Xt(data=Obj$datamat,p=Obj$p, weight.matrix=Obj$weight,ic="AIC")

#type=Obj$type
#out2=vars::VAR(DATA$Xt,p=2,type=type,lag.max =5,ic="SC")


#myIRF1=irf(out2, impulse =impulse , response=response,ortho = ortho,cumulative = cumulative,boot=boot,ci = 0.95,n.ahead=15,shockNeg=shockNeg,STD=STD)

#if (G==TRUE)
#{

#  if (is.null(Weights))
#  {
#  y=apply(as.data.frame(myIRF1[[1]]),1,mean)
#  Low=apply(as.data.frame(myIRF1[[2]]),1,mean)
#  Up=apply(as.data.frame(myIRF1[[3]]),1,mean)
#  dat=cbind(y,Low,Up)
#  plot.ts(y,ylim=c(min(dat),max(dat)),col=4);abline(h=0)
#  lines(Low,lwd=0.1,lty=4,col=2)
#  lines(Up,lwd=0.1,lty=4,col=2)
#  }
#  else {
#  yw=as.matrix(as.data.frame(myIRF1[[1]]))%*%Weights
#  Loww=as.matrix(as.data.frame(myIRF1[[2]]))%*%Weights
#  Upw=as.matrix(as.data.frame(myIRF1[[3]]))%*%Weights
#  datw=cbind(yw,Loww,Upw)
#  plot.ts(yw,ylim=c(min(datw),max(datw)),col=4);abline(h=0)
#  lines(Loww,lwd=0.1,lty=4,col=2)
#  lines(Upw,lwd=0.1,lty=4,col=2)
#  }
#print(myIRF1)
#} ## boot endif
#else {print(myIRF1);plot(myIRF1)}
#print(out2$p)
#return(result=out2)

#}  ## end of irf function





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


averageCORgvar <-function (out) {
var.no=out$K
varRSD=out$varRSD
V1=t(matrix(colnames(varRSD),var.no,))

gvarRSD=out$gvarRSD
V=t(matrix(colnames(gvarRSD),var.no,))

varRSD.cor=list();gvarRSD.cor=list()
for (j in 1:var.no) {
temp1=(apply(cor(varRSD[,V1[,j]]),2,sum)-1)/(length(V1[,j])-1)
varRSD.cor[[j]]=temp1

temp=(apply(cor(gvarRSD[,V[,j]]),2,sum)-1)/(length(V[,j])-1)
gvarRSD.cor[[j]]=temp
}
results=list(varRSDcor=varRSD.cor,gvarRSDcor=gvarRSD.cor)
return(results)

}
