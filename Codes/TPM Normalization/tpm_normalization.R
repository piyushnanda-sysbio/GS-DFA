library(dplyr)
library(plyr)
library(tidyselect)

enslen=read.csv('EnsembleIDsLungIsoformsLength.txt')
colnames(enslen)=c('GENEID','Length')
enslen=aggregate(Length ~ GENEID,enslen,mean)

expr=read.csv('GSE147507.csv')

colnames(expr)[1]='SYMBOL'

expr_matrix=join(expr,enslen,type='right')

r_tpm <- function(dfr,len)
{
  dfr1 <- sweep(dfr,MARGIN=1,(len/10^4),`/`)
  scf <- colSums(dfr1,na.rm=TRUE)/(10^6)
  return(sweep(dfr1,2,scf,`/`))
}

tpmexpr=r_tpm(expr_matrix[3:80],expr_matrix[,81])

tpm=cbind(tpmexpr[,1:78],'GENEID'=expr_matrix[,1])

#Subset metabolic genes

mapping=read.csv('ensemble_id_mapping.csv')

mettpm=join(tpm,mapping,by='GENEID',type='right')

write.csv(mettpm,'MetTPMAll.csv')