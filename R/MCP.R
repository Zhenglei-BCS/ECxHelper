
stepDownTrendTest.CATT <- function (x, g, n,test = c("catt_exact"), alternative = c("two.sided", "greater", "less"), continuity = FALSE, ...) 
{
  if (is.list(x)) {
    if (length(x) < 2L) 
      stop("'x' must be a list with at least 2 elements")
    DNAME <- deparse(substitute(x))
    x <- lapply(x, function(u) u <- u[complete.cases(u)])
    k <- length(x)
    l <- sapply(x, "length")
    if (any(l == 0)) 
      stop("all groups must contain data")
    g <- factor(rep(1:k, l))
    if (!is.null(x$alternative)) 
      alternative <- x$alternative
    x <- unlist(x)
  }
  else {
    if (length(x) != length(g)) 
      stop("'x' and 'g' must have the same length")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
    OK <- complete.cases(x, g)
    x <- x[OK]
    g <- g[OK]
    n <- n[OK]
    if (!all(is.finite(g))) 
      stop("all group levels must be finite")
    g <- factor(g)
    k <- nlevels(g)
    if (k < 2) 
      stop("all observations are in the same group")
  }
  alternative <- match.arg(alternative)
  test <- match.arg(test)
  lev <- levels(g)
  res <- lapply(2:k, function(i) {
    ok <- g %in% lev[1:i]
    browser()
    inp <- list( cases = x[ok], dose.ratings = g[ok], totals=n[ok])
    do.call(test, inp)
  })
  browser()
  METHOD <- paste0("Step down ", res[[1]]$method)
  pval <- sapply(1:(k - 1), function(j) res[[j]]$p.value)
  PVAL <- matrix(pval, nrow = k - 1, ncol = 1, byrow = TRUE)
  colnames(PVAL) <- lev[1]
  rownames(PVAL) <- lev[2:k]
  stat <- sapply(1:(k - 1), function(j) res[[j]]$statistic)
  STAT <- matrix(stat, nrow = k - 1, ncol = 1, byrow = TRUE)
  colnames(STAT) <- lev[1]
  rownames(STAT) <- lev[2:k]
  dist <- names(res[[1]]$statistic)
  ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL, 
              statistic = STAT, alternative = alternative, p.adjust.method = "none", 
              dist = dist)
  class(ans) <- "trendPMCMR"
  ans
}



#' Summary Williams Test results.
#'
#' @param object 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
summaryZG <-function (object, ...) 
{
  critVal <- as.numeric(object$crit.value)
  stat <- as.numeric(object$statistic)
  dist <- object$dist
  dec <- ifelse(stat > critVal, "reject", "accept")
  critDist <- paste0(dist, "-crit")
  if (!is.matrix(object$statistic)) {
    if (grepl(pattern = "Hayter's", x = object$method)) {
      H0 <- switch(object$alternative, greater = paste("Mean(xi) - Mean(xj) <= 0"), 
                   less = paste("Mean(xi) - Mean(xj) >= 0"))
    }
    else {
      H0 <- switch(object$alternative, greater = paste("Med(xi) - Med(xj) <= 0"), 
                   less = paste("Med(xi) - Med(xj) >= 0"))
    }
  }
  else {
    grp1 <- as.numeric(c(col(object$statistic)))
    grp2 <- as.numeric(c(row(object$statistic)))
    cnam <- colnames(object$statistic)
    rnam <- rownames(object$statistic)
    H0 <- switch(object$alternative, greater = paste(rnam[grp2], 
                                                     "-", cnam[grp1], "<=", "0"), less = paste(rnam[grp2], 
                                                                                               "-", cnam[grp1], ">=", "0"))
    ok <- !is.na(stat)
    stat <- stat[ok]
    H0 <- H0[ok]
    dec <- dec[ok]
  }
  dist <- paste0(dist, "-value")
  cat("\n\t", object$method, "\n\n")
  cat("data: ", object$data.name, "\n")
  if (!is.null(object$alternative)) {
    cat("alternative hypothesis: ", object$alternative, "\n")
  }
  paramName <- names(object$parameter)
  if (length(paramName) == 2) {
    suppressWarnings(expr = xdf <- data.frame(STATISTIC = round(stat, 
                                                                3), PARAM1 = object$parameter[1], PARAM2 = object$parameter[2], 
                                              CRITDIST = round(critVal, 3), DECISION = dec, ALPHA = 0.05))
    names(xdf) <- c(dist, paramName[1], paramName[2], critDist, 
                    "decision", "alpha")
  }
  else {
    suppressWarnings(expr = xdf <- data.frame(STATISTIC = round(stat, 
                                                                3), PARAM1 = object$parameter[1], CRITDIST = round(critVal, 
                                                                                                                   3), DECISION = dec, ALPHA = 0.05))
    names(xdf) <- c(dist, paramName[1], critDist, "decision", 
                    "alpha")
  }
  rownames(xdf) <- H0
  return(xdf) 
  invisible(object)
}

#' get from william res accept/reject
#'
#' @param william 
#'
#' @return
#' @export
#'
#' @examples
getwilliamRes <- function(william){
  if(class(william)=="try-error"){
    return(rep(NA,3))
  }else(return(as.character(summaryZG(william)$decision)))
}

#' getEndpoint according to a series of p-values
#'
#' @param pvals 
#' @param doses 
#' @param procedure 
#'
#' @return
#' @export
#'
#' @examples
getEndpoint <- function(pvals,doses=c("Control","B","C","D"),procedure="stepDown"){
  nd <-length(doses)
  if(all(is.na(pvals))) {
    pvals <- rep(0,nd-1) ## length of pvalues
  }
  np <- length(pvals) 
  if(np < (nd-1)){
    pvals <- c(pvals,rep(0,nd-1-np))
  }
  sig <- pvals<0.05
  
  
  NOEC <- doses[1]
  if(procedure=="stepDown")
  {
    for(i in (nd-1):1){
      if(!sig[i]){ NOEC <- doses[i+1]
      
      return(NOEC)}
    }
    if(sig[1]) return(NOEC)
    
  }else{
    for(i in 1:(nd-1)){
      if(sig[i]){ 
        NOEC <- doses[i]
        return(NOEC)
      }
      
    }
    if(!sig[nd-1]) return(doses[nd])
  }
}

#' get Endpoint for continuous data
#'
#' @param paov 
#' @param pks 
#' @param pnormal 
#' @param phomogeneity 
#' @param monotonicity 
#' @param william 
#' @param dunnett 
#' @param dunn 
#' @param jonckheere 
#' @param procedure 
#' @param doses 
#'
#' @return
#' @export
#'
#' @examples
contEndpoint <- function(paov,pks,pnormal,phomogeneity,monotonicity,william,dunnett,dunn,jonckheere,procedure="stepDown",doses=c("A","B","C","D")){
  ## if monotone not rejected. linear contrast not significant quadratic sig!
  test <- "NK"
  if(is.character(monotonicity)){
    monotonicity <- as.numeric(gsub("<","",monotonicity))
  }
  if(monotonicity[1] >0.05 & monotonicity[2] < 0.05){
    if(pnormal>0.05 & phomogeneity>0.05){
      ## check if PAVA Problem!
      pvals <- dunnett$p.value
      test <- "Dunnett"
    }else{
      pvals <- dunn$p.value
      test <- "Dunn"
    }
  }else{
    if(pnormal>0.05 & phomogeneity>0.05){
      ## check if PAVA Problem!
      if(class(william)!="try-error"){
        pvals <- as.character(summaryZG(william)$decision)
        pvals <- as.numeric(plyr::mapvalues(pvals,from=c("accept","reject"),to=c(0.2,0.01)))
        test <- "Williams"
      }else{
        pvals <- jonckheere$p.value
        test <- "Jonckheere"
      } 
    }else{
      pvals <- jonckheere$p.value
      test <- "Jonckheere"
    }
  }
  ############
  NOEC <- getEndpoint (pvals,doses=doses,procedure=procedure)
  attributes(NOEC) <- list(test=test)
  return(NOEC)
}


#' Testing Monocity
#'
#' @param Data 
#' @param Treatment 
#' @param Response 
#'
#' @return
#' @export
#'
#' @examples
monotonicityTest <-function(Data,Treatment,Response){
    #' @export
    #This is the test for monotonicity as done in the SAS Version
    #it sets up linear and Quadratic contrast for an ANOVA
    #Uses .stdEndEnv$WeightVar 
    Data <- as.data.frame(Data)
    #Transform data for rank response 
    Data<-rankTransform(Data,Response)
    
    
    #Form Contrasts
    LineContrast<-getLineContrast(Data,Treatment)
    QuadContrast<-getQuadContrast(Data,Treatment)
    Contrasts <- cbind(LineContrast,QuadContrast)
    colnames(Contrasts)<-c('Line','QuadContrast')
    
    #ANOVA
    Data[ ,Treatment]<-as.factor(Data[ ,Treatment])
    contrasts(Data[ ,Treatment])<-Contrasts
    AnovaTable<-aov(Data[ ,'TransformedResponse']~as.factor(Data[ ,Treatment]))
    # if (is.null(.stdEndEnv$WeightVar)==FALSE){  #if there is a weight and the data is not averaged 
    #   if (length(.stdEndEnv$WeightVar)==length(Data[ ,Response])){
    #     AnovaTable<-aov(Data[ ,'TransformedResponse']~as.factor(Data[ ,Treatment]), weight=.stdEndEnv$WeightVar )
    #   }
    # }
    # 
    #gather information and clean the table
    CAnova<-summary.lm(AnovaTable)
    MonocityTable<-as.data.frame(CAnova$coefficients[2:3,3:4])
    rownames(MonocityTable)<-NULL
    MonocityTable<-cbind(c('Linear','Quadratic'),MonocityTable,'.')
    colnames(MonocityTable)[1]<-'Test'
    colnames(MonocityTable)[4]<-'Significance'
    MonocityTable$Significance<-as.character(MonocityTable$Significance)
    
    MonocityTable$Significance[MonocityTable[ ,'Pr(>|t|)']<0.05]<-'*'
    MonocityTable$Significance[MonocityTable[ ,'Pr(>|t|)']<0.01]<-'**'
    MonocityTable$Significance[MonocityTable[ ,'Pr(>|t|)']<0.001]<-'***'
    
    MonocityTable[ ,'Pr(>|t|)']<-round(MonocityTable[ ,'Pr(>|t|)'],4)
    if (length(which(MonocityTable[ ,'Pr(>|t|)']<10^-4))>0){
      MonocityTable[which(MonocityTable[ ,'Pr(>|t|)']<10^-4),'Pr(>|t|)']<-'<0.0001'
    }
    MonocityTable[ ,'t value']<-round(MonocityTable[ ,'t value'],2)
    
    return(MonocityTable)
  }

getLineContrast <-
  function(Data,Treatment){
    #' @export
    #This attains the contrast used to test a linear relationship 
    K<-nlevels(as.factor(Data[ ,Treatment]))
    
    #Contrasts
    switch (K,
            return(0),
            return(c(-1,1)),
            return(c(-1,0,1)),
            return(c(-3,-1,1,3)),
            return(c(-2,-1,0,1, 2)),
            return(c(-5, -3, -1,1,3,5)),
            return(c(-3,-2, -1,0,1,2,3)),
            return(c(-7,-5,-3,-1,1,3,5,7)),
            return(c(-4,-3,-2,-1,0,1,2,3,4)),
            return(c(-9,-7,-5,-3,-1, 1,3, 5, 7,9)),
    )
    return()
  }

getQuadContrast <-
  function(Data,Treatment){
    #' @export
    #This attains the contrast used to test a quadratic relationship 
    
    K<-nlevels(as.factor(Data[ ,Treatment]))
    #Contrasts
    switch (K,
            return(0),
            return(c( 0,0)),
            return(c( 1,-2,1)),
            return(c(1,-1,-1,1)),
            return(c(2,-1,-2,-1,2)),
            return(c(5,-1,-4,-4,-1,5)),
            return(c(5,0,-3,-4,-3,0,5)),
            return(c(7,1,-3,-5,-5,-3,1,7)),
            return(c(28,7,-8,-17,-20,-17,-8,7,28)),
            return(c(6,2,-1,-3,-4,-4,-3,-1,2,6)),
    )
    return()
  }
rankTransform <-
  function(Data,VecName){
    #' @export
    #This function will rank transform the data the same way the SAS code below does
    
    #proc rank data=sorted out=assumtst normal=blom 
    #ties=mean;
    #var Counts;
    #by generation;   
    #ranks rank_EGGS;
    #run;
    
    #Data- data set to be modified 
    #VecName Name of Variable  to be modified
    #returns a data set with TransformedResponse as the rank transformed responce
    
    
    Vector<-Data[ ,VecName];
    
    #This handles the rank transform 
    n<-length(Vector)
    r.i<-rank(Vector, na.last = TRUE, ties.method = c("random"))
    y.i<-(r.i-3/8)/(n+1/4) 
    Data$y.i.norm<-qnorm(y.i)  #this will be for output
    Data<-Data[with(Data, order(Data$y.i.norm)), ]
    
    
    Repeats<-xtabs(~as.factor(Data[ ,VecName]))
    #This Block of code is for finding ties and will break if no ties exist 
    if(max(Repeats)>1){ #Check for ties
      Repeats<-Repeats[-which(Repeats==1)] #Remove non-ties
      Names<-dimnames(Repeats)
      for (e in Names[[1]]){  #For every number that has a tie
        e<-as.numeric(e)
        Data[which(Data[ ,VecName]==e),'y.i.norm']<-mean(Data[which(Data[ ,VecName]==e),'y.i.norm']) #Average
      }
    }
    colnames(Data)[which(colnames(Data)=='y.i.norm')]<-'TransformedResponse' #rename column
    
    return(Data)
  }


RSCABK <-  function(x.i.j,n.i.j,m.i,TestK,test.type){
    #This function is called for the detailed results functions 
    #This is just for 1 slice 
    #' @export
    K=1
    #Variable names are the same as in Rao and Scott 1992
    
    x.i<-colSums(x.i.j)  #Number of successful observations 
    n.i<-colSums(n.i.j)  #Total Number in a treatment 
    
    p.i.hat=x.i/n.i
    
    
    #Calculate r.ij.sum.sq
    #this uses the expanded form of sigma(r.ij) squared
    r.ij.sum.sq<-mat.or.vec(length(m.i)[1],K)
    
    p.i.hat.k<-apply(t(p.i.hat),2,rep,dim(n.i.j)[1])  #p.hat 
    
    if(test.type =='CA'){
      d.i<-1  #No adjustment
    }
    if(test.type =='RS'){
      #Calculate the Rao-Scott adjustment
      x.ij.sum.sq<-colSums(x.i.j^2)
      cross<- -2*colSums(x.i.j*n.i.j*p.i.hat.k)
      last.sum.sq<-colSums((n.i.j*p.i.hat.k)^2)
      r.ij.sum.sq<-x.ij.sum.sq+cross+last.sum.sq;
      
      v.i <- m.i/(m.i - 1)/n.i^2 *r.ij.sum.sq
      d.i <- n.i * v.i/(p.i.hat * (1 - p.i.hat))
      
      
      
      d.i[which(is.na(d.i))]<-1    #Ignores the correction when incidence is 0% or 100%
      d.i[which(is.nan(d.i))]<-1   #Ignores the correction when incidence is 0% or 100%
      d.i[which(d.i<1)]<-1         #Cf Rao and Scott Clustered Poisson paper;  
    }
    
    #apply adjustment
    x.i.new <- x.i/d.i
    n.i.new <- n.i/d.i
    
    scores<-matrix(1:length(x.i.new)-1)
    p.hat <- sum(x.i.new)/sum(n.i.new)
    mean.score <- sum(scores * n.i.new)/sum(n.i.new)
    var.scores <- sum(n.i.new * (scores - mean.score)^2)
    
    RS <- (sum(x.i.new * scores) - p.hat * sum(n.i.new * scores))/sqrt(p.hat * (1 - p.hat) * var.scores) #the statistic value 
    #P-Value
    p.val <- pnorm(abs(RS), lower.tail = TRUE)  
    p.val<-1-p.val
    if (is.finite(p.val)==FALSE){
      p.val<-1
    }
    #Adds marks for significance 
    Sig<-rep('.',length(p.val))
    Sig[which(p.val<=.05)]<-'*'
    Sig[which(p.val<=.01)]<-'**'
    Sig[which(p.val<=.001)]<-'***'
    
    return (list('Treatment'=max(scores)+1,'R-Score'=TestK,'Statistic' = RS, 'P-Value' = p.val,'Signif'=Sig))
  }

stepDownRSCABS <-
  function(TestK,x.i.j,n.i.j,m.i,Effect,test.type){
    Next<-0; #Step Down test
    #' @export
    #Replicates are rows, Treatment are columns  [Replicate,Treatment]
    Result.K<-as.data.frame({})
    x.i.j.K<-x.i.j[ , ,TestK]
    while (Next==0){
      Result<-RSCABK(x.i.j.K,n.i.j,m.i,TestK,test.type=test.type)
      
      #The step down
      Next<- Result['P-Value']>0.05 | dim(x.i.j.K)[2]<=2 #Stop conditions of too high of p.val or 2 treatments
      x.i.j.K<-x.i.j.K[ ,-dim(x.i.j.K)[2]]
      n.i.j<-n.i.j[ ,-dim(n.i.j)[2]]
      m.i<-m.i[-length(m.i)]
      Result.K<-rbind(Result.K, as.data.frame(Result))
    }
    
    
    Result.K<-rbind(Result.K, as.data.frame(Result))
    
    Effect<-paste(Effect,TestK,sep='')
    Result.K<-cbind(Effect,Result.K)
    return (Result.K)
  }


runRSCABS <-  function(Data,Treatment,Replicate='',Effects='',test.type='RS'){
    #This function will produce a table of step-down Cochran-Armitage trend tests with possible Rao-Scott adjustment by slices
    #It will Run the test on every effect in the Effect list 
    #' @export
    #Turn Replicate and Treatment into factors
    
    Data[ ,Treatment]<-as.factor(Data[ ,Treatment])
    if ( identical(Replicate,'') == FALSE && identical(test.type,'CA')==FALSE){
      Data[ ,Replicate]<-as.factor(Data[ ,Replicate])
    }else{
      test.type='CA'
      Replicate<-'FalseReplicateFill'
      Data$FalseReplicateFill<-as.factor(1:length(Data[ ,Treatment])) #Needed so R data structure do not change type. 
    }
    
    
    if (test.type != 'RS' & test.type != 'CA'){
      message('Error: Invalid test type')
      return()
    }
    #Remove all non-whole numbers
    
    
    #Default for effect every column name that is not a Treatment, Replicate, and has a 0 < K.max < 20
    #turn off warnings
    options(warn=-1)
    if (Effects ==''){
      Effects<-colnames(Data)
      Maxes<-apply(Data,2,max,na.rm=TRUE)
      Remove<-which(as.numeric(Maxes) <= 0 | as.numeric(Maxes) >  20 | is.na(as.numeric(Maxes))==TRUE)
      Remove<-c(Remove,which(Effects==Treatment),which(Effects==Replicate))
      Remove<-c(Remove,which(is.factor(Data[ ,Effects])))
      Remove<-unique(Remove)
      if (length(Remove)>0){
        Effects<-Effects[-Remove]
      }
    }
    
    #Changed to account for just 1 effect
    if (length(Effects)>1){
      Data[ ,Effects]<-apply(Data[ ,Effects],2,convert2Score)
    }else{
      Data[ ,Effects]<-convert2Score(Data[ ,Effects])
    }
    
    options(warn=0) #turn on warnings
    
    
    #Need to remove factors
    #Prep Data
    
    
    Data.Prep<-sapply(c('&Fill#',Effects),prepDataRSCABS,Data=Data,Treatment=Treatment,Replicate=Replicate)	
    Results.Raw<-sapply(c('&Fill#',Effects),stepKRSCABS,Data=Data.Prep,Treatment=Treatment,Replicate=Replicate,test.type=test.type)
    Results<-do.call("rbind", lapply(Results.Raw, data.frame, stringsAsFactors = FALSE))
    
    
    if (length(which(is.finite(Results[  ,'Statistic'])==FALSE))>0){
      Results<-Results[-which(is.finite(Results[  ,'Statistic'])==FALSE), ]
    }
    rownames(Results)<-{}
    
    return(Results)
}
convert2Score <-
  function(Dvec){
    #This Function will convert any object that is not zero or a positive number to NA 
    #' @export
    oldw <- getOption("warn")
    options(warn=-1)  #Turn off warnings
    Dvec<-as.numeric(as.character(Dvec))  #Will induce a a warning 
    Dvec[which(Dvec<0)]<-NA
    options(warn=oldw) #Turn on warnings
    return(Dvec)
  }

prepDataRSCABS <-
  function(Effect='',Data={},Treatment='',Replicate=''){
    #Data Transform from a list of individuals to matrix format 
    #This will take Clustered Data and convert it to by 
    #' @export
    
    if(Effect=='&Fill#'){ #ensures the correct data structure output
      return()
    }
    if(length(which(colnames(Data)==Effect))==0){
      print(paste(Effect,' is not in data set. Ending function.',sep=''))
      return()
    }
    K.max<-max(Data[ ,Effect],na.rm=TRUE) #max K score
    
    
    #Remove NA and negative numbers
    if (length(which(is.na(Data[ ,Effect])))>0){
      Data<-Data[-which(is.na(Data[ ,Effect])), ]
    }
    if (length(which(Data[ ,Effect]<0))>0){
      Data<-Data[-which(Data[ ,Effect]<0), ]
    }
    #Convert Factors to Numerics
    Data[ ,Effect]<-as.numeric(Data[ ,Effect])
    
    if (K.max==0){
      print(paste('There is no variation in ',Effect,'. Ending function.',sep=''))
      return()
    }
    #Replicates are rows, Treatment are columns  [Replicate,Treatment]
    
    n.i.j<-xtabs( ~Data[[Replicate]]+Data[[Treatment]])
    m.i<-apply(n.i.j,2,function(Vec){
      if (length(which(Vec==0))>0){
        Vec<-Vec[-which(Vec==0)]
      }
      return(length(Vec))
    })
    
    
    
    
    x.i.j<-array(dim=c(dim(n.i.j)[1],dim(n.i.j)[2],K.max)) #Declare x.i.j , frequency array of scores k or larger 
    #Each k level is on the 3rd dimension
    for (K in 1:K.max){
      x.i.j[ , ,K]<-xtabs( ~Data[[Replicate]]+Data[[Treatment]],subset=Data[ ,Effect]>=K)
    }
    RSCABS.Prep.Data<-list(x.i.j=x.i.j,n.i.j=n.i.j,m.i=m.i,K.max=K.max)
    return(RSCABS.Prep.Data)
  }

stepKRSCABS <-
  function(Effect,Data.Prep,Treatment,Replicate,test.type){
    #This iterates though each k-level
    #' @export
    if (is.null(Data.Prep[[Effect]])==TRUE){ #Skip effects not needed
      return()
    }
    
    x.i.j<-Data.Prep[[Effect]]$x.i.j
    n.i.j<-Data.Prep[[Effect]]$n.i.j
    m.i<-Data.Prep[[Effect]]$m.i
    K.max<-Data.Prep[[Effect]]$K.max
    
    Results.Effect<-lapply(1:K.max,stepDownRSCABS,x.i.j=x.i.j,n.i.j=n.i.j,m.i=m.i,Effect=Effect,test.type=test.type) #each k-level
    Results.Effect<-do.call("rbind", lapply(Results.Effect, data.frame, stringsAsFactors = FALSE))  #'fix' output table
    return(Results.Effect)
  }


###### CATTexact