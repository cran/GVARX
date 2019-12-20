GVAR_Ft <- function (data, weight.matrix=NULL)  {

ID<-NULL
NAME=unique(as.character(data[,"ID"]))
timeCol=which(names(data)=="Time")
idCol=which(names(data)=="ID")

timeID=as.character(subset(data, ID==NAME[1])[,"Time"])
Year=unique(as.character(lubridate::year(timeID)))
dat=data[,-c(1:2)]
variates=colnames(dat) #names of column variables
N=length(NAME) #Number of countries

Ft=list()

for (i in 1:N) {

 exo=NULL  # Compute Exogenous Foreign Variables
 for (j in 1:length(variates)) {

    if (is.null(weight.matrix)){

    F.tmp=apply(matrix(dat[,j], ,N),1,mean)

   }
   if (isTRUE(is.matrix(weight.matrix))) {

        varMatrix=matrix(dat[,j], ,N)
        F.tmp=varMatrix %*% as.matrix(weight.matrix[,i])


      } else {

        dat_matrix=matrix(dat[,j], ,N)
        varMatrix=as.xts(dat_matrix,as.Date(timeID))

        F.tmp=NULL
        for (k in 1:length(Year)) {
        F.tmp0=varMatrix[Year[k]] %*% as.matrix(weight.matrix[[k]][,i])
        F.tmp=rbind(F.tmp,F.tmp0)
        }

        }
 exo=cbind(exo, F.tmp) # Foreign variables
 }
 h=NAME[i]  # index number of home country

exoNAMES=paste0(h,variates,"_FL0")

 colnames(exo)=exoNAMES

Ft[[i]]=exo
}

 return(Ft)

}
