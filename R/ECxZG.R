

checkEDres <- function(edResTab,...){
  ## ToDo!!!!!!!!!!!!!!!!!!!!! what do I want to do I forgot!!
}
#' adding ECx estimation and interval to the model output plot
#'
#' @param p model plot as ggplot2 output
#' @param object fitted drc object
#' @param EDres EDresults corresponding to respLev
#' @param trend "Decrease" or "Increase"
#' @param endpoint ErC or EbC or EyC or EC or LD
#' @param respLev reponse levels
#' @param textAjust.x  where to label ECx
#' @param textAjust.y where to label ECx
#' @param useObsCtr whether to use observed control mean
#' @param d0 to be used control mean. 
#' @param textsize label text size
#' @param lineheight errorbar height
#' @param xmin when confidence intervals very wide including even negative values.
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
addECxCI <- function(p=NULL,object,EDres=NULL,trend="Decrease",endpoint="ErC", respLev=c(10,20,50),
                     textAjust.x=0.1,textAjust.y=0.05,useObsCtr=FALSE,d0=NULL,textsize=2,lineheight=1,xmin=0.05,...){
  if(respLev[1]<1) respLev <- respLev*100
  
  if(!inherits(object,"try-error")){
    type <- object$type
    coefs <- coef(object)
    num <- length(coefs)
    if("c:(Intercept)" %in% names(coefs)) cVal <- coefs["c:(Intercept)"] else cVal <-0 ## switch(response,Increase=0,Decrease=1)
    if("d:(Intercept)" %in% names(coefs)) d <- coefs["d:(Intercept)"] else d <- 1
    # For 3 parameter models, if setting the 4th parameter being 0, then there is no need to change ECx estimation unless the observed control average should be used instead of the estimated control mean. 
    ## As shown in the equation below, c=0 as assumed, d is the estimated control mean. 
    ## EC_true x=EC_software {((d-(1-x)· mean(y_0 ))⁄((d-c)})
    if(d < cVal){
      d1 <- d
      d <- cVal
      cVal <- d1
    }
    if(type=="continuous" || type=="Poisson"){
      if(trend=="Increase"){
        ## for example, mortality can increase, however, 
        ## the x% effect should be compared to the maximal survival probabilities.
        ## which means, (1-cVal)=max survival probability.
        ## There should be only reduction or inhibition.
        x.relative <- (respLev-cVal*100)/(d-cVal)
        ### need to think about it!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        warning("If not inhibition or reduction, be careful! ")
      }else{
        if(!useObsCtr) d0 <- d
        x.relative <- d0*respLev/(d-cVal)
      }
    }else{
      if(type=="binomial"){
        if(num>2) {
          warning("need more study")
          x.relative <-respLev
        }else{
          x.relative <- respLev
        }
      }
    }
    
  }else{
    stop("No need to plot as error occured in the model fitting already ")
  }
  
  
  if(is.null(EDres)) EDres <- ED.ZG(object,trend=trend,respLev=respLev,...)
  ##browser()
  if(is.null(p)) p <- plot.modList(object,...)
  if(trend=="Decrease"){
    Effect=((1-respLev/100)*d)
  }else{
    warning("If not inhibition or reduction, be careful! ")
    Effect=cVal+respLev/100*d ## note, for survival, d=1, for continuous, it is a different topic!!!!
  }
  EDdat <- data.frame(EDres,Effect=Effect,
                           text=paste0(endpoint,respLev))
  EDdat$Lower <- ifelse(EDdat$Lower>xmin,EDdat$Lower,xmin)
  p <- p+geom_point(data=EDdat,aes(x=(Estimate),y=Effect),size=3,col="red")+
    geom_text(data=EDdat,aes(x=Estimate+textAjust.x,y=Effect+textAjust.y,label=text),size=textsize,col="darkgreen")+
    geom_errorbarh(data=EDdat,aes(x=Estimate,y=Effect,xmin=Lower,xmax=Upper,height=lineheight),lty=2)
  return(p)
}

#' Calculating ED following the regulatrory ED definition.
#'
#' @param object fitted model using drm
#' @param respLev response level, the x of ECx.
#' @param maxEff Not used for now. maximum effect.
#' @param trend Increase of Decrease, whether the dose response is decreasing or increasing compared control.
#' @param range using Percentage, not used for now.
#' @param interval methods to calculate the confidence intervals!
#' @param ... other parameters for ED function
#'
#' @return Back calculated regulatory ECx
#' @export ED.ZG
#'
#' @examples
ED.ZG <- function(object,respLev,maxEff=TRUE,trend="Increase",range="Percentage",CI=c("delta","inv"),...){
  ## Note that this might not be suitable for models with fixed c or d parameters, where
  CI <- match.arg(CI)
  if(grepl("log(ED50) as parameter",object$fct$text,fixed = TRUE)) CI <-"delta"
  if(!inherits(object,"try-error")){
    if(CI=="delta"){
      type <- object$type
      coefs <- coef(object)
      num <- length(coefs)
      if(grepl("log(ED50) as parameter",object$fct$text,fixed = TRUE)) {
        interval<- "fls" 
        ##browser()
      }else interval <- "delta"
      if("c:(Intercept)" %in% names(coefs)) cVal <- coefs["c:(Intercept)"] else cVal <-0 ## switch(response,Increase=0,Decrease=1)
      if("d:(Intercept)" %in% names(coefs)) d <- coefs["d:(Intercept)"] else d <- 1
      # For 3 parameter models, if setting the 4th parameter being 0, then there is no need to change ECx estimation unless the observed control average should be used instead of the estimated control mean. As shown in the equation below, c=0 as assumed, d is the estimated control mean. EC_true x=EC_software {((d-(1-x)· mean(y_0 ))⁄((d-c)})
      if(d < cVal){
        print("check this case: using decay model to fit growth or vice versa.") ## using decay model to fit growth or vice versa.
        d1 <- d
        d <- cVal
        cVal <- d1
        #browser()
      }
      if(type=="continuous" || type=="Poisson"){
        if(trend=="Increase"){
          x.relative <- (respLev-cVal*100)/(d-cVal)
          
        }else{
          
          ## x.relative <- d-(1-respLev/100)*y0/(d-c)
          x.relative <- d*respLev/(d-cVal)
        }
      }else{
        if(type=="binomial"){
          if(num>2) {
            warning("need more study, survival as response should be fine, but not mortality as response")
            x.relative <-respLev
          }else{
            x.relative <- respLev
          }
        }
      }
      ## Add a step to check if x.relative is out of range!
      id <- (x.relative>100)
      if(all(id)){
        res <- matrix(NA,length(respLev),4)#
        colnames(res) <- c("Estimate", "Std. Error","Lower","Upper")
      }else{
        res <- ED(object,respLev = x.relative[!id] ,display=FALSE,interval=interval,...)
        if(sum(id)>0){
          tmp <- matrix(NA,(sum(id)),ncol(res))
          ##browser()
          res <-rbind(res,tmp)
        }
      }
    }else{
      
      require(bmd)
      # tmp <- lapply(cross2(respLev,modList),function(l){
      #   x <- l[[1]]
      #   y <- l[[2]]
      #   res <- try(bmd(y,x/100,backgType="modelBased",def="relative",interval="inv",display=F))
      #   if(!inherits(res,"try-error")){
      #     return(c(res$Results[1],res$SE,res$interval))
      #   }else{
      #     return(rep(NA,4))
      #     }
      #   
      # })
      tmp <- lapply(respLev,function(l){
        x <- l[[1]]
        y <- object
        if(object$type=="binomial"){
         ## browser()
          res <- try(bmd(y,x/100,backgType="modelBased",def="excess",interval="delta",display=F))
        } else{
          res <- try(bmd(y,x/100,backgType="modelBased",def="relative",interval="inv",display=F))
        }
        if(!inherits(res,"try-error")){
          return(c(res$Results[1],res$SE,res$interval))
        }else{
          return(rep(NA,4))
        }
        
      })
      res <- plyr::ldply(tmp)
      colnames(res) <- c("Estimate", "Std. Error","Lower","Upper")
      
    }


  }else{
    res <- matrix(NA,length(respLev),4)#
    colnames(res) <- c("Estimate", "Std. Error","Lower","Upper")
  }

  rownames(res) <- paste("EC",respLev)
  return(res)

}
#' Title
#'
#' @param object
#' @param fctList
#' @param nested
#' @param sorted
#' @param linreg
#' @param icfct
#' @param respCol
#' @param doseCol
#' @param data
#' @param type
#' @param EFSA
#'
#' @return mselect.ZG
#' @export
#'
#' @examples
mselect.ZG <- function (object=NULL, fctList = NULL, nested = FALSE, sorted = c("IC",
                                                                                "Res var", "Lack of fit", "no"), linreg = FALSE, icfct = AIC,respCol="effect",doseCol="dose",data=NULL,type="continuous",EFSA=NULL)
{
  sorted <- match.arg(sorted)
  if (!is.logical(nested)) {
    stop("'nested' argument takes only the values: FALSE, TRUE")
  }
  contData <- identical(object$type, "continuous")

  mc <- match.call()

  ## If continuos data, we need to calculate ResVar, if nested, we need to do an ANOVA step.
  fctNames <- sapply(fctList,function(f)f$name)

  if(is.null(object)) {
    if(is.null(data)){stop("No Data Input, please check your input parameters")}
    modList <- lapply(fctList, function(fct){
      try(drm( data[["response"]] ~ data[["dose"]], fct = fct, type=type, control=drmc(errorm = T)), silent=T)
    })
  }else{
    modList <- lapply(fctList, function(fct){
      try(update(object, fct = fct),
          silent = TRUE)})
    modList[[length(fctList)+1]] <- object
    ## fctList[[length(fctList)+1]] <- mod$fct possbile bug
    fctList[[length(fctList)+1]] <- object$fct  ## not necessary!!
  }

  Comparison <- plyr::laply(modList, function(tempObj){
    if (!inherits(tempObj, "try-error")) {
      criteria <- c(logLik(tempObj),icfct(tempObj),modelFit(tempObj)[2, 5])
      if (contData) {
        tryRV2 <- try(summary(tempObj)$resVar, silent = TRUE)
        if (inherits("tryRV2", "try-error")) {
          tryRV2 <- NA
        }

        criteria <- c(criteria,tryRV2)
      }
    }else{
      print("Model not Converged, Please consult a statistician.")
      if(contData) criteria <- rep(NA,4) else criteria <- rep(NA,3)
    }
    return(criteria)
  })


  rownames(Comparison) <- sapply(fctList,function(f)f$name)
  cnames <- c("logLik", "IC", "Lack of fit")
  if (contData) {
    cnames <- c(cnames, "Res var")
  }

  colnames(Comparison) <- cnames
  if (linreg) {
    if(is.null(object)){object <- modList[[1]]}
    drcData <- as.data.frame(object$data[, c(2, 1)])
    names(drcData) <- c("yVec", "xVec")
    linFitList <- list(lm(yVec ~ xVec, data = drcData), lm(yVec ~
                                                             xVec + I(xVec * xVec), data = drcData), lm(yVec ~
                                                                                                          xVec + I(xVec * xVec) + I(xVec * xVec * xVec), data = drcData))
    linModMat <- matrix(unlist(lapply(linFitList, function(listObj) {
      c(logLik(listObj), icfct(listObj), NA, (summary(listObj)$sigma)^2)
    })), 3, 4, byrow = TRUE)
    rownames(linModMat) <- c("Lin", "Quad", "Cubic")
    colnames(linModMat) <- cnames[1:4]

    Comparison <- rbind(Comparison, linModMat)
  }
  names(modList) <- rownames(Comparison)
  if (sorted != "no") {
    res <- list(Comparison=Comparison[order(Comparison[, sorted]), ],modList=modList[order(Comparison[, sorted])],EFSA=EFSA)

  }
  else {
    res <- (list(Comparison=Comparison,modList=invisible(modList),EFSA=EFSA))
  }

  class(res) <- "drcComp"
  return(res)
}
print.drcComp <- function(x,..){
  x=x[c("Comparison","EFSA")]
  NextMethod()
}
#' Plot a list of models together. 
#'
#' @param modList
#' @param respLev
#' @param data
#' @param xmin
#' @param xmax
#' @param scale
#' @param npts
#' @param ... user input when the use does not want to use the data provided in the modList
#'
#' @return
#' @export plot.modList
#'
#' @examples
plot.modList <- function(modList,respLev=NULL,data=NULL,xmin,xmax,scale=c("logx","logy","logxy","orig"),
                         npts=100,plot_respLev=FALSE,xbreaks=NULL,ymin=NULL,ymax=NULL,...){
  if(is.null(data)){
    mod1 <- modList[[1]]
    i <- 1
    while(inherits(mod1,"try-error")){
      mod1 <- modList[[i+1]]
      i <- i+1
    }
    data <- mod1$origData
    responseName <- mod1$dataList$names$orName
    doseName <- mod1$dataList$names$`dName`
    dose <- mod1$dataList$`dose`
    response <- mod1$dataList$origResp
  }else{
    print("need user input of responseName and doseName")
    if(missing(doseName)) stop("need user input of responseName and doseName!")
  }
  ## fctNames <- names(modList)
  ## Calculate Prediction and Confidence Bands
  ex1 <- data.frame(data)
  ex1$conc0 <- pull(data,doseName)
  if(missing(xmin)) xmin <- min(ex1$conc0[ex1$conc0>0])/10
  ex1$conc0[ex1$conc0 == 0] <- xmin
  ex1$response <- pull(data,responseName)
  if(missing(xmax)) xmax <- max(ex1[,doseName])*1.1
  if(is.null(xbreaks)){
    DoseData <- data.frame(Dose = exp(seq(log(xmin), log(xmax),length = npts)))
  }else{
    DoseData <- data.frame(Dose=xbreaks)
  }
  names(DoseData) <- doseName
  predData <- plyr::ldply(modList,function(tempObj){
    if (!inherits(tempObj, "try-error")) {
      pm <- predict(tempObj, newdata = DoseData,
                    interval = "confidence",...)
      newdata <- DoseData
      newdata$p <- pm[, 1]
      newdata$pmin<- pm[, 2]
      newdata$pmax<- pm[, 3]
      newdata$Model <- tempObj$fct$name
      if(!is.null(ymin)) newdata$pmin[newdata$pmin<ymin] <- ymin
      if(!is.null(ymax)) newdata$pmax[newdata$pmax>ymax] <- ymax
    }else{
      newdata <- NULL
    }
    return(newdata)
  })
  predData$Model <- factor( predData$Model,levels=unique(predData$Model))
 
  p <- ggplot(ex1, aes_(x = ~conc0, y = as.name(responseName))) +
    geom_ribbon(data = predData, aes_(x = as.name(doseName), y = ~ p, ymin = ~pmin,ymax = ~pmax,col=~Model,fill=~Model),alpha = 0.2) +
    xlab(doseName) + ylab(responseName) +geom_line(data = predData,aes_(x = as.name(doseName), y = ~ p,col=~Model)) + geom_point()

  if(!is.null(respLev) & plot_respLev){
    y0 <- mean(response[dose==0])## calculate the control mean
    p <- p + geom_hline(yintercept=y0*(100-respLev)/100,lty=2)
  }
  scale <- match.arg(scale)
  if(scale=="logx") p <- p+scale_x_log10()
  if(scale=="logy") p <- p+scale_y_log10()
  if(scale=="logxy") p <- p+scale_x_log10()+scale_y_log10()
  return(p)
}

#' Plot the ECx estimation and confidence intervals from the list of models.
#'
#' @param edList
#' @param fctNames
#' @param ...
#'
#' @return
#' @export plot.edList
#'
#' @examples
plot.edList <- function(edList,fctNames,...){
  if(!is.data.frame(edList)){
    edResTab <- plyr::ldply(lapply(edList,function(x) {
      x <- as.data.frame(x)
      if(ncol(x)==3) {
        x$`Std. Error`<- rep(NA,nrow(x))
        x <- x[,c(1,4,2,3)]
      }
      x$ED <- rownames(x)
      rownames(x) <- NULL
      x
    }))
    if(is.null(names(edList))){
      edResTab$Model <- fctNames
    }else edResTab$Model <- names(edList)
  }else{
    ## Need to check if edRes Tab
    edResTab <- edList
    edResTab$Model <- edResTab$.id
  }
  edResTab$Model <- factor(edResTab$Model,levels=unique(edResTab$Model))
  if(is.null(edResTab$Rating)){
    plotED <- ggplot(edResTab, aes(x = ED, y = Estimate, col = Model)) +
      geom_pointrange(aes(ymin = Lower, ymax = Upper),
                      position = position_dodge(width = 0.4))#+scale_y_log10()
  }else{
    edResTab$Rating <- factor(edResTab$Rating,levels=c("Excellent","Good","Fair","Poor","Bad"))

    plotED <- ggplot(edResTab, aes(x = EC, y = Estimate, col = Model,shape=Rating)) +
      geom_pointrange(aes(ymin = Lower, ymax = Upper),
                      position = position_dodge(width = 0.4))#+scale_y_log10()
  }
  return(plotED)
}
## Steepness and overlap
#' Title
#'
#' @param mod
#' @param obj
#' @param trend
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
calcSteepnessOverlap <- function(mod=NULL,obj=NULL,trend="Decrease",...){
  # |Overlapping  Conditions                  |Certainty of the Protection Level |
  # |-----------------------------------------|----------------------------------|
  # |EC$_{10}$ < EC$_{20,low}$                | High                             |
  # |EC$_{20,low}$ < EC$_{10}$ < EC$_{50,low}$| Medium                           |
  # |EC$_{10}$ > EC$_{50,low}$                | Fair                             |
  if(!is.null(mod)){
    if(class(mod)=="drc") obj <- ED.ZG(mod,c(10,20,50),trend=trend,...)
  }
  res <- rep(NA,2)
  obj <- as.data.frame(obj)
  steep <- obj$Estimate[1]/obj$Estimate[3] # the Ratio between EC 10 and EC 50
  if(!is.na(steep)){
    if(steep <0.33) res[2] <- "Shallow" else{
      if(steep <0.66) res[2] <- "Medium" else res[2] <- "Steep"
    }
  }else steep <- "Not Defined"
  if(nrow(obj)!=3){
    print("Without EC50, not able to calculate Steepness")
  }else{


    if(!is.na(obj$Lower[2]) & !is.na(obj$Estimate[1])){
      if(obj$Estimate[1] < obj$Lower[2]){
        res[1] <- "High"
      }else{
        if(!is.na(obj$Lower[3])){
          if(obj$Estimate[1] < obj$Lower[3])res[1] <- "Medium" else res[1] <-"Fair"
        }else{
          res[1] <- "Not Low"
        }
      }
    }else{
      res[1] <- "Not Defined"
    }

  }
  res
}

#' Select ECx estimation from Models with mselect.ZG and EFSA SO criteria
#'
#' @param modList
#' @param respLev
#' @param trend
#' @param ...
#'
#' @return
#' @export mselect.ED
#'
#' @examples
mselect.ED <- function(modList,respLev=c(10,20,50),trend="Decrease",...){
  edList <- lapply(modList,function(obj) ED.ZG(obj,respLev,trend = trend,...))
  ##map2(modList,edList,function(x,y){predict(x,newdata=data.frame(Concentration=y[,1]))}) ## Check if Prediction Make sense

  NW <- lapply(edList, calcNW)
  edRes <- map2(edList,NW,cbind)


  edResTab <- plyr::ldply(lapply(edRes,function(x) {
    x <- as.data.frame(x)
    if(ncol(x)==3) {
      x$`Std. Error`<- rep(NA,nrow(x))
      x <- x[,c(1,4,2,3)]
    }
    x$EC <- rownames(x)
    rownames(x) <- NULL
    x
  }))
  return(edResTab)
}

## Normalised width
#' Title
#'
#' @param x
#' @param ED
#'
#' @return
#' @export
#'
#' @examples
calcNW <- function(x,ED="ZG"){
  x <- as.data.frame(x)  ## x is the ED object from ED function.
  out <- as.data.frame((x$Upper-x$Lower)/x$Estimate)
  if(ED=="ZG") newRowNames <- rownames(x)
  if(ED=="drc") newRowNames <- paste0("EC_",do.call(c, lapply(strsplit(row.names(x), ":"), function(y) y[3])))
  out <- cbind(out,NA)
  out[,2] <- ECx_rating(out[,1])
  row.names(out) <- newRowNames
  names(out) <- c("NW","Rating")
  return(out)
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
ECx_rating <- function(x){
  ## normalised width (NW) of CI
  ## NW   |  Rating
  ## ----------------
  ## <0.2 | Excellent
  ## <0.5 | Good
  ## <1   | Fair
  ## <2   | Poor
  ## >2   | Bad
  outRating <- NULL
  for (i in x){
    if (is.nan(i) || is.na(i)){
      outRating <- c(outRating,"not defined")
    }else{
      if (i<0.2){
        outRating <- c(outRating,"Excellent")
      }else if (i <0.5){
        outRating <- c(outRating,"Good")
      }else if (i <1){
        outRating <- c(outRating,"Fair")
      }else if (i <2){
        outRating <- c(outRating,"Poor")
      }else if (i >2){
        outRating <- c(outRating,"Bad")
      }else{
        outRating <- c(outRating,"not defined")
      }
    }
  }
  return(outRating)
}

#' Compare models with multiple criteria
#'
#' @param modRes
#' @param modList
#' @param trend
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
drcCompare <- function(modRes=NULL,modList=NULL,trend = "Decrease",CI=c("delta","inv"),...){
  CI <- match.arg(CI)
  if(!is.null(modRes)){
    ## take the comarison matrix from modRes
    modList <- modRes$modList
    comparison <- modRes$Comparison
  }else{
    if(is.null(modList)) stop("Need the model output list from previous step!")
  }
  #A significance test is provided for the comparison of the dose-response model considered and the simple linear regression model with slope 0 (a horizontal regression line corresponding to no dose effect)
  # The likelihood ratio test statistic and the corresponding degrees of freedom and p-value are reported.
  a <- plyr::ldply(modList,function(obj) if(inherits(obj,"try-error")) c(`Chi-square test` = NA, Df = NA, `p-value` = NA) else noEffect(obj))
  Overlap <- plyr::ldply(modList,calcSteepnessOverlap,CI=CI)
  names(Overlap) <- c(".id","Certainty_Protection","Steepness")
  CompRes <- cbind(comparison,Overlap[,-1])
  CompRes <- cbind(CompRes,a[,"p-value"])
  colnames(CompRes)[ncol(CompRes)] <- "No Effect p-val"
  return(CompRes)
}

#' for scaling x-axis
#'
#' @param x 
#' @param a 
#'
#' @return
#' @export
#'
#' @examples
logxp <- function(x,a){log(a+x)}
#' inverse transformation of logxp
#'
#' @param x 
#' @param a 
#'
#' @return
#' @export
#'
#' @examples
invlogxp <- function(x,a){exp(x)-a}

## logxp_trans <- function(a=0.001) scales::trans_new("logxp",transform=function(x)logxp(x,a=a),inverse=function(x)invlogxp(x,a=a))
