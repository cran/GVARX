GVAR_Xt <- function(data,p=2,type="const",ic="AIC",weight.matrix=NULL){
#data=Data;type="const";ic="AIC";p=1;FLag=2;lag.max=NULL
ID <- NULL
  type=type
ic=ic

idCol=which(colnames(data) == "ID")
timeCol=which(colnames(data)=="Time")

dat1=data[,-timeCol] #Data with ID column only

endo.no=ncol(dat1)-1 #names of column variables
N=length(unique(dat1[,idCol])) #Number of countries
weight.matrix=weight.matrix
FLag=p+1
lag.max=NULL
myout= GVARest(data=data,p,FLag,lag.max=NULL, type=type,ic=ic, weight.matrix)


pmatrix=myout$lagmatrix[,2]
p=min(pmatrix)

cat("\n","Number of lag for eventual Xt is", p, "\n")


GO.tmp=NULL
G1=NULL
G2=NULL
RESID=NULL

if (p==1) {
if (is.list(weight.matrix)) {
  AVG=matrix(rep(0,N^2),N,N)
    for (i in 1:length(weight.matrix)){ 
      AVG=AVG+as.matrix(weight.matrix[[i]])
      }
      weight.matrix=AVG/(length(weight.matrix)-1)
      } else {weight.matrix=as.matrix(weight.matrix)}
  
  for (jj in 1:N) {
  rsd_tmp=resid(myout$gvar[[jj]])
  diff_jj=max(pmatrix)-(pmatrix[jj])
    if (diff_jj==0) { rsd_tmp1=rsd_tmp}  else {rsd_tmp1=rsd_tmp[-(1:diff_jj),]}

  RESID=cbind(RESID,rsd_tmp1)

##===== Compute G0

  endo_lagk0=NULL
  endo_lagk1=NULL
  endo_lagk2=NULL
  for (k in 1:endo.no)  {
  coeff=coef(myout$gvar[[jj]])[[k]]
  exo_no=myout$exoLag*endo.no
if (myout$type=="const") {
    coeff_EXO=coeff[(which(rownames(coeff)=="const")+1):nrow(coeff),]
} else if (myout$type %in%  c("trend","both")) {
  coeff_EXO=coeff[(which(rownames(coeff)=="trend")+1):nrow(coeff),]
} else {coeff_EXO=coeff[(endo.no*p+1):nrow(coeff),]}
  
    coeff_EXO_LAG0=coeff_EXO[1:endo.no,1]

     exo_lagi0=NULL;exo_lagk1=NULL
     for (i in 1:endo.no)
     {
     exo_lagi0=rbind(exo_lagi0,coeff_EXO_LAG0[i]*as.matrix(weight.matrix)[,jj])
     } #end of i loop
endo_lagk0=cbind(endo_lagk0,c(exo_lagi0))

#===== Compute G1: Lags components
coeff_EXO_LAG=coeff_EXO[-(1:endo.no),1]
coeff.ENDO=coeff[1:(pmatrix[jj]*endo.no),1]
coeff.ENDO.LAG1=coeff.ENDO[1:endo.no]
coeff_EXO_LAG1=coeff_EXO_LAG[1:endo.no]

     endo_lagi1=NULL
     for (i in 1:endo.no)
     {
     endo_lagi1=rbind(endo_lagi1,coeff_EXO_LAG1[i]*as.matrix(weight.matrix)[,jj])
     } #end of i loop
    tmp.lagi1=matrix(0,endo.no,N)
    tmp.lagi1[,jj]=coeff.ENDO.LAG1
    tmp1=tmp.lagi1+endo_lagi1
endo_lagk1=cbind(endo_lagk1,c(tmp1))

} #end of k loop
 GO.tmp=cbind(GO.tmp,endo_lagk0)
 G1=cbind(G1,endo_lagk1)
} #end of jj loop

G0=diag(1,endo.no*N)-GO.tmp
invGO=solve(G0)
F1=invGO%*%G1

newRESID=t(invGO%*%t(RESID))

varnames=colnames(dat1)[-1]
NAME=myout$NAMES

dataNT=vnames=NULL
for (j in 1:N) {
dat=subset(dat1,ID==NAME[j])
vnames=c(vnames,paste(NAME[j],varnames,sep="."))
colnames(dat)=NULL
datz=as.matrix(dat[,-1])
dataNT=cbind(dataNT,datz)
}
colnames(dataNT)=vnames

##== Recursive Procedure
#myXt_tmp1=.GVARfilter(dataNT,p=1,Bcoef=F1)
#myXt1=myXt_tmp1$xfitted;
#t1=nrow(myXt1);t_rsd=nrow(newRESID)
#Xt1=myXt1[-(1:(t1-t_rsd)),]+newRESID
#colnames(Xt1)=vnames
colnames(newRESID)=vnames

results <-list(lagmatrix=myout$lagmatrix,G0=G0,G1=G1,F1=F1,RESID=RESID,newRESID=newRESID)
# end of if (p=1)


} else if (p>=2) {
  if (is.list(weight.matrix)) {
    AVG=matrix(rep(0,N^2),N,N)
    for (i in 1:length(weight.matrix)){ 
      AVG=AVG+as.matrix(weight.matrix[[i]])
    }
    weight.matrix=AVG/(length(weight.matrix)-1)
    
    } else {weight.matrix=as.matrix(weight.matrix)}
  
  for (jj in 1:N) {

    #Collect residuals
    rsd_tmp=resid(myout$gvar[[jj]])
    diff_jj=max(pmatrix)-(pmatrix[jj])
    if (diff_jj==0) {rsd_tmp1=rsd_tmp } else {rsd_tmp1=rsd_tmp[-(1:diff_jj),]}

    RESID=cbind(RESID,rsd_tmp1)

    ##===== Compute G0

    endo_lagk0=endo_lagk1=endo_lagk2=NULL
    for (k in 1:endo.no) {
      coeff=coef(myout$gvar[[jj]])[[k]]
      exo_no=myout$exoLag*endo.no
      if (myout$type=="const") {
        coeff_EXO=coeff[(which(rownames(coeff)=="const")+1):nrow(coeff),]
      } else if (myout$type %in%  c("trend","both")) {
        coeff_EXO=coeff[(which(rownames(coeff)=="trend")+1):nrow(coeff),]
      } else {coeff_EXO=coeff[(endo.no*p+1):nrow(coeff),]}
      coeff_EXO_LAG0=coeff_EXO[1:endo.no,1]

      exo_lagi0=NULL
      exo_lagk1=NULL
      for (i in 1:endo.no)
      {
        exo_lagi0=rbind(exo_lagi0,coeff_EXO_LAG0[i]*as.matrix(weight.matrix)[,jj])
      } #end of i loop
      endo_lagk0=cbind(endo_lagk0,c(exo_lagi0))

      #===== Compute G1: Lags components
      coeff_EXO_LAG=coeff_EXO[-(1:endo.no),1]
      coeff.ENDO=coeff[1:(pmatrix[jj]*endo.no),1]
      coeff.ENDO.LAG1=coeff.ENDO[1:endo.no]
      coeff_EXO_LAG1=coeff_EXO_LAG[1:endo.no]

      endo_lagi1=NULL
      for (i in 1:endo.no)  {
        endo_lagi1=rbind(endo_lagi1,coeff_EXO_LAG1[i]*as.matrix(weight.matrix)[,jj])
      } #end of i loop

      tmp.lagi1=matrix(0,endo.no,N)
      tmp.lagi1[,jj]=coeff.ENDO.LAG1
      tmp1=tmp.lagi1+endo_lagi1
      endo_lagk1=cbind(endo_lagk1,c(tmp1))

      #===== Compute G2: Lags components
      coeff.ENDO.LAG2=coeff.ENDO[(endo.no+1):(endo.no*p)]
      coeff_EXO_LAG2=coeff_EXO_LAG[(endo.no+1):(endo.no*p)]

            endo_lagi2=NULL
      for (i in 1:endo.no)  {
        endo_lagi2=rbind(endo_lagi2,coeff_EXO_LAG2[i]*as.matrix(weight.matrix)[,jj])
      } #end of i loop

      tmp.lagi2=matrix(0,endo.no,N)
      tmp.lagi2[,jj]=coeff.ENDO.LAG2
      tmp2=tmp.lagi2+endo_lagi2
      endo_lagk2=cbind(endo_lagk2,c(tmp2))

    } #end of k loop

    GO.tmp=cbind(GO.tmp,endo_lagk0)
    G1=cbind(G1,endo_lagk1)
    G2=cbind(G2,endo_lagk2)
  } #end of jj loop

  G0=diag(1,endo.no*N)-GO.tmp
  invGO=solve(G0)

  F1=invGO%*%G1
  F2=invGO%*%G2

  newRESID=t(invGO %*% t(RESID))
  varnames=colnames(dat1)[-1]
  NAME=myout$NAMES

  dataNT=NULL;vnames=NULL
  for (j in 1:N) {
    dat=subset(dat1,ID==NAME[j])
    vnames=c(vnames,paste(NAME[j],varnames,sep="."))
    colnames(dat)=NULL
    datz=as.matrix(dat[,-1])
    dataNT=cbind(dataNT,datz)
  }
  colnames(dataNT)=vnames

  ##== Recursive Procedure
#  myXt_tmp1=.GVARfilter(dataNT,p=1,Bcoef=F1)
#  myXt_tmp2=.GVARfilter(dataNT,p=2,Bcoef=cbind(F1,F2))
#  myXt1=myXt_tmp1$xfitted
#  myXt2=myXt_tmp2$xfitted
#  t1=nrow(myXt1)
#  t2=nrow(myXt2)
#  t_rsd=nrow(newRESID)
#  Xt1=myXt1[-(1:(t1-t_rsd)),]+newRESID
#  Xt2=myXt2[-(1:(t2-t_rsd)),]+newRESID
#  colnames(Xt1)=vnames
#  colnames(Xt2)=vnames
  colnames(newRESID)=vnames
  results <-list(G0=G0,G1=G1,G2=G2, F1=F1,F2=F2,lagmatrix=myout$lagmatrix,RESID=RESID,newRESID=newRESID)

} # if p=2



return(results)
}






.GVARfilter <- function(X, p, Bcoef, exogen = NULL, postpad = c("none", "constant", "zero", "NA"))
{
	X = as.matrix(X)
	if(any(is.na(X))) stop("\nvarxfilter:-->error: NAs in X.\n")
	if(ncol(X) < 2) stop("\nvarxfilter:-->error: The matrix 'X' should contain at least two variables.\n")
	if(is.null(colnames(X))) colnames(X) = paste("X", 1:ncol(X), sep = "")
	colnames(X) = make.names(colnames(X))
	postpad = tolower(postpad[1])
	if(any(colnames(Bcoef)=="const")){
		constant = TRUE
		ic = 1
	} else{
		constant = FALSE
		ic = 0
	}
	obs = dim(X)[1]
	K = dim(X)[2]
	xsample = obs - p
	Xlags = embed(X, dimension = p + 1)[, -(1:K)]
	temp1 = NULL
	for (i in 1:p) {
		temp = paste(colnames(X), ".l", i, sep = "")
		temp1 = c(temp1, temp)
	}
	colnames(Xlags) = temp1
	Xend = X[-c(1:p), ]
	if(constant){
		rhs = cbind( Xlags, rep(1, xsample))
		colnames(rhs) <- c(colnames(Xlags), "const")
	} else{
		rhs = Xlags
		colnames(rhs) <- colnames(Xlags)
	}
	if( !(is.null(exogen)) ) {
		exogen = as.matrix(exogen)
		if (!identical(nrow(exogen), nrow(X))) {
			stop("\nvarxfit:-->error: Different row size of X and exogen.\n")
		}
		XK = dim(exogen)[2]
		if (is.null(colnames(exogen))) colnames(exogen) = paste("exo", 1:ncol(exogen), sep = "")
		colnames(exogen) = make.names(colnames(exogen))
		tmp  = colnames(rhs)
		rhs =  cbind(rhs, exogen[-c(1:p), ])
		colnames(rhs) = c(tmp, colnames(exogen))
	} else{
		XK = 0
	}
	datamat = as.matrix(rhs)
	colnames(datamat) = colnames(rhs)
	xfitted = t( Bcoef %*% t( datamat ) )
	xresiduals = tail(X, obs - p) - xfitted
	if(postpad!="none"){
		if(postpad == "constant"){
			# pre-pad values with the constant
			xfitted = t( Bcoef %*% t( rbind(matrix(c(rep(0, p*K), if(constant) 1 else NULL, if(XK>0) rep(0, XK) else NULL), nrow = p, ncol=dim(Bcoef)[2], byrow = TRUE), datamat ) ) )
			xresiduals = X - xfitted
		} else if(postpad == "zero"){
			xfitted = t( Bcoef %*% t( rbind(matrix(rep(0, dim(Bcoef)[2]), nrow = p, ncol=dim(Bcoef)[2], byrow = TRUE), datamat ) ) )
			xresiduals = X - xfitted
		} else if(postpad == "NA"){
			xfitted = t( Bcoef %*% t( rbind(matrix(rep(NA, dim(Bcoef)[2]), nrow = p, ncol=dim(Bcoef)[2], byrow = TRUE), datamat ) ) )
			xresiduals = X - xfitted
		} else{
			# do nothing
			xfitted = t( Bcoef %*% t( datamat ) )
			xresiduals = tail(X, obs - p) - xfitted
		}
	}
	ans = list( Bcoef = Bcoef, xfitted = xfitted, xresiduals = xresiduals, lag = p, constant = constant)
	return( ans )
}



