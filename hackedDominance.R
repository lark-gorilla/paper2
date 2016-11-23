# hacked version of dominanace analyses script so that it works
# homemade function dm2 is basically the same but omits two r2 metrics
# that cause the current code to fail with my model


dominanceMatrix<-function(x,undefined.value=0.5) {
  vars<-colnames(x)
  m<-length(vars)
  ma<-matrix(undefined.value,m,m,dimnames=list(vars,vars))
  
  for(i in 1:(m-1)) {
    for(j in (i+1):m) {
      comps<-na.omit(cbind(x[,i,drop=F],x[,j,drop=F]))
      if(mean(comps[,1]>comps[,2])==1) 
      {
        ma[i,j]<-1
        ma[j,i]<-0
      }
      
      if(mean(comps[,1]<comps[,2])==1) 
      {
        ma[i,j]<-0
        ma[j,i]<-1
      }
      
    }
  }
  ma
}

daAverageContributionByLevel<-function(x) {
  ff<-x$fit.functions
  
  out<-list()
  for(i in ff) {
    
    res<-aggregate(x$fits[[i]],list(level=x$level),mean,na.rm=T)
    
    out[[i]]<-res[res$level<max(x$level),]
    
  }
  out
}

daCompleteDominance<-function(daRR) {
  lapply(daRR$fits,dominanceMatrix)
}

daConditionalDominance<-function(daRR) {
  daACBL<-daAverageContributionByLevel(daRR)
  analize<-function(x) {
    x<-x[,-1]
    dominanceMatrix(x)
  }
  lapply(daACBL,analize)
}

daGeneralDominance<-function(daRR) {
  daACBL<-daAverageContributionByLevel(daRR)
  analize<-function(x) {
    x<-x[,-1]
    gm<-matrix(colMeans(x),1,ncol(x),dimnames=list(1,colnames(x)))
    dominanceMatrix(gm)
  }
  lapply(daACBL,analize)
}


#slightly hacked version, had to remove 2 of the r2
# metrics (Nagelkerke and Cox and Snell) as throwing NAs

dm2<-function (x, constants = c(), fit.functions = "default", data = NULL, 
               null.model = NULL, ...) 
{
  daModels <- daSubmodels(x, constants)
  daRaw4r2 <- daRawResults(x, constants, fit.functions, data, 
                           null.model, ...)
  daRaw<-daRaw4r2
  daRaw$fit.functions<-daRaw$fit.functions[c(1,4)]
  daRaw$fits<-daRaw$fits[c(1,4)]
  daRaw$base.fits<-daRaw$base.fits[,c(1,4)]
  
  daAverageByLevel <- daAverageContributionByLevel(daRaw)
  daAverageGeneral <- lapply(daAverageByLevel, function(x) {
    colMeans(x[, -1])
  })
  list(predictors = daModels$predictors, constants = daModels$constants, 
       fit.functions = daRaw$fit.functions, fits = daRaw, contribution.by.level = daAverageByLevel, 
       contribution.average = daAverageGeneral, complete = daCompleteDominance(daRaw), 
       conditional = daConditionalDominance(daRaw), general = daGeneralDominance(daRaw))
}
