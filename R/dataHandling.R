

#' Reading raw ToxRat/.. data according to Data type.
#'
#' @param path path to the file (excel sheets)
#' @param sheet sheet of relevant endpoint
#' @param DataType  either "Algae", "Myriophyllum", "Lemna","Earthworm", or ""
#'
#' @return a data table of long format, suitable for drm function
#' @export
#'
#' @examples
Raw2Long <- function( path,sheet="RawData Fresh Weight",DataType="Myriophyllum"){
  rawdat <- readxl::read_excel(path.expand(path),sheet=sheet)
  ## remove all NA columns.
  ind <- sapply(1:ncol(rawdat),function(i){
    if(all(is.na(rawdat[i]))) FALSE else TRUE
  })
  # cols <- 1:ncol(rawdat)
  # cols<- cols[ind]
  rawdat <- rawdat[,ind]

  ## Note in the raw data, the first row are units

  trts <- as.vector(t(rawdat[2,3:ncol(rawdat)]))
  rawdat <- rawdat %>% fill(...2)
  rawdat <- rawdat[-(1:2),]
  names(rawdat) <- c("Time","Replicates",trts)
  id <- apply(rawdat[,trts],1,function(x)!(all(is.na(x))))
  rawdat <- rawdat[id,]
  ## if Replicate is specified, then there is no need to do 
  if(all(is.na(as.numeric(rawdat$Replicates)))){
    ## specific data type in the case of "~/ECx/Algae/Data2/Copy of NewOECD201-M-002353-02-1-Bowers-Skele.xlsx"
    dat.wide <- rawdat %>%  filter(is.na(Replicates))
    
  }else{
    dat.wide <- rawdat %>% mutate(Replicates=as.numeric(Replicates))%>% filter(!is.na(Replicates))
    
  }
  
 if(DataType=="Myriophyllum"){ ## For this we need to average over pseudo-replicates
   nRep <- nrow(dat.wide)/length(unique(dat.wide$Time))/3
   nTime <- length(unique(dat.wide$Time))
   dat.wide<- dat.wide %>% mutate(Replicates=rep(rep(1:nRep,each=3),nTime))
   dat <- dat.wide %>% gather(Treatment,Response,-c(Time,Replicates)) %>% 
     mutate(Concentration=recode(as.numeric(Treatment),.missing=0),Response=as.numeric(Response),Time=as.numeric(Time)) %>% 
     filter(!is.na(Response)) %>% group_by(Time,Replicates,Treatment,Concentration) %>% 
     summarise(Response=mean(Response)) %>% ungroup
   
 } else{
 
 
   dat <- dat.wide %>% gather(Treatment,Response,-c(Time,Replicates)) %>% 
     mutate(Concentration=dplyr::recode(as.numeric(Treatment),.missing=0),Response=as.numeric(Response),Time=as.numeric(Time)) %>% 
     filter(!is.na(Response)) 
 }
  time0 <- dat[dat$Time==0,] %>% dplyr::select(Replicates,Treatment,Concentration,Response) %>% rename(Initial=Response)
  if(DataType %in% c("NTA")){
    dat <- left_join(dat,time0) %>% mutate(Survival=Response/Initial) %>% mutate(Treatment=factor(Treatment,levels=unique(Treatment))) 
  }else{
    if(DataType %in% c("Earthworm","Bees")){
      dat <- left_join(dat,time0)%>% mutate(Treatment=factor(Treatment,levels=unique(Treatment))) 
      ## Note, initial could be NA!!
    }else{
      dat <- left_join(dat,time0) %>% mutate(Yield=Response-Initial) %>% mutate(Treatment=factor(Treatment,levels=unique(Treatment))) %>% 
        mutate(GR=(log(Response)-log(Initial))/(Time))
    }
    
  }
  
  return(dat)
}
