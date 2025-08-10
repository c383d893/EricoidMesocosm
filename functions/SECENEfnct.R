#Creating Function to calculate CE, SE, and NE; also SE and N functions----------------
#****************************************************************
#SECENE calculation is described in Box 1:

#Loreau, M., Hector, A., 2001. Partitioning selection and complementarity in biodiversity experiments. Nature 412, 72–76.

SECE<-function(df,PLOTNO="",Mi="",Yoi="",Yo="",RYej=""){
  PLOTNOcolno<-which(colnames(df)==PLOTNO)
  SECE<-data.frame(df[,PLOTNOcolno])
  names(SECE)[1]<-paste("PLOTNO")
  Micolno<-which(colnames(df)==Mi)   
  SECE$M.i<-df[,Micolno]
  
  #Y[oj] is the observed yield of speces i in mixture
  Yoicolno<-which(colnames(df)==Yoi)   
  SECE$Y.oj<-df[,Yoicolno]
  
  #Y[o] total observed yield of the mixture
  Yocolno<-which(colnames(df)==Yo)   
  SECE$Y.o<-df[,Yocolno]
  
  #RY[Ej] is the expected relative yield of species i in the mixture, which is simply its proportion seeded/planted
  RYejcolno<-which(colnames(df)==RYej)   
  SECE$RY.ej<-df[,RYejcolno]
  
  #RY[oj] is the observed relative yield of species i in the mixture
  SECE$RY.oj<-SECE$Y.oj/SECE$M.i
  
  #Y[ej] is the expected yield of species i in the mixture
  SECE$Y.ej<-SECE$RY.ej*SECE$M.i
  
  #Y[e] is the total expected yield of the mixture
  Y.e<-aggregate(SECE$Y.ej,by=list(SECE$PLOTNO),sum)
  names(Y.e)[1:2]<-paste(c("PLOTNO","Y.e"))
  SECE<-left_join(SECE,Y.e,by="PLOTNO")
  remove(Y.e)
  
  #delta Y is the deviation from total expeced yield in the mixture
  SECE$deltaY<-SECE$Y.o-SECE$Y.e
  
  #delta RY[i] is the deviation from expected relative yield of species i in the mixture
  SECE$deltaRY.i<-SECE$RY.oj-SECE$RY.ej
  
  #N is the number of species in the mixture
  SECE$N<-df$SPPNO
  
  #deltaRYbar is the average deltaRY for all species in mixture
  deltaRYbar<-aggregate(SECE$deltaRY.i,by=list(SECE$PLOTNO),mean)
  names(deltaRYbar)[1:2]<-paste(c("PLOTNO","deltaRYbar"))
  SECE<-left_join(SECE,deltaRYbar,by="PLOTNO")
  remove(deltaRYbar)
  
  #Mbar is the average monoculture yield for all species in mixture
  Mbar<-aggregate(SECE$M.i,by=list(SECE$PLOTNO),mean)
  names(Mbar)[1:2]<-paste(c("PLOTNO","Mbar"))
  SECE<-left_join(SECE,Mbar,by="PLOTNO")
  remove(Mbar)
  
  #CE is the complementary effect which is N * deltaRYbar * Mbar
  SECE$CE<-SECE$N*SECE$deltaRYbar*SECE$Mbar
  
  #SE is the selection effect which is deltaY - CE
  SECE$SE<-SECE$deltaY-SECE$CE
  df<-cbind(df,SECE)
  remove(SECE)
  return(df)
}

NESECE<-function(df,PLOTNO=""){
  PLOTNOcolno<-which(colnames(df)==PLOTNO)
  col1<-which(colnames(df)=="deltaY")
  col2<-which(colnames(df)=="CE")
  col3<-which(colnames(df)=="SE")
  col4<-which(colnames(df)=="Y.o")
  NESECE<-data.frame(df[,c(PLOTNOcolno,col1,col2,col3,col4)])
  names(NESECE)[1]<-paste("PLOTNO")
    
  mNESECE<-NESECE %>% group_by(PLOTNO) %>% 
      summarise(mean(Y.o),
                mean(deltaY),
                mean(CE),
                mean(SE))
  names(mNESECE)[2:5]<-paste(c("PROPTOTPLOT","deltaY","CE","SE"))
  return(mNESECE)
  remove(NESECE,mNESECE,col1,col2,col3,col4,PLOTNOcolno)
}

meanNESECE<-function(df,treatment=""){
  trtcolno<-which(colnames(df)==treatment)
  col1<-which(colnames(df)=="deltaY")
  col2<-which(colnames(df)=="CE")
  col3<-which(colnames(df)=="SE")
  col4<-which(colnames(df)=="PROPTOTPLOT")
  
  df2<-data.frame(df[,c(trtcolno,col1,col2,col3,col4)])
  names(df2)[1]<-paste("treatment")
  
  meandf<-data.frame(
    aggregate(df2$PROPTOTPLOT,by=list(df2$treatment),mean),
    tapply(df2$PROPTOTPLOT,df2$treatment,sd),
    tapply(df2$deltaY,df2$treatment,mean),
    tapply(df2$deltaY,df2$treatment,sd),
    tapply(df2$CE,df2$treatment,mean),
    tapply(df2$CE,df2$treatment,sd),
    tapply(df2$SE,df2$treatment,mean),
    tapply(df2$SE,df2$treatment,sd)
  )
  names(meandf)[1:9]<-paste(c(treatment,"Cover","sdCover","deltaY","sddeltaY","CE","sdCE","SE","sdSE"))
  return(meandf)
  remove(df,df2,trtcolno,col1,col2,col3,col4)
  }

#****************************************************************
#Phi function is a measure of community synchrony, and is described:

#Loreau, M., de Mazancourt, C., 2008. Species synchrony and its drivers: neutral and nonneutral community dynamics in fluctuating environments. Am. Nat. 172, E48–E66. 

#Husse, S., Huguenin-Elie, O., Buchmann, N. and Lüscher, A., 2016. Larger yields of mixtures than monocultures of cultivated grassland species match with asynchrony in shoot growth among species but not with increased light interception. Field Crops Research, 194, pp.1-11. 

PhiFUNCTION<-function(df,PLOTNO="PLOTNO",SppID="SPPID",TotYield="PLOTCOVER",SppYield="SPPCOVER"){
  col1<-which(colnames(df)==PLOTNO) 
  col2<-which(colnames(df)==TotYield) 
  col3<-which(colnames(df)==SppYield) 
  col4<-which(colnames(df)==SppID) 
  
  df$PLOTNOSPPID<-paste(df[,col1],df[,col4],sep=",")
  df2<-aggregate(df[,col2],by=list(df[,col1]),var)
  df3<-aggregate(df[,col3],by=list(df$PLOTNOSPPID),sd)
  df3<-separate(data = df3, col = Group.1, into = c("PLOTNO", "SPPID"), sep = ",")
  df3$PLOTNO<-as.numeric(df3$PLOTNO)
  df4<-aggregate(df3$x,by=list(df3$PLOTNO),sum)
  df4$x2<-df4$x^2
  df5<-left_join(df2,df4[,c(1,3)],by="Group.1")
  df5$phi<-df5$x/df5$x2
  df1<-df5[,c(1,4)]
  names(df1)[1]<-paste(PLOTNO)
  return(df1)
  remover(df1,df2,df3,df4,df5)
}

N<-function(x){length(x)}
se<-function(x){sd(x)/sqrt(length(x))}

sigfig <- function(vec, n=3){ 
  ### function to round values to N significant digits
  # input:   vec       vector of numeric
  #          n         integer is the required sigfig  
  # output:  outvec    vector of numeric rounded to N sigfig
  
  formatC(signif(vec,digits=n), digits=n,format="fg", flag="#") 
  
}      # end of function   sigfig

