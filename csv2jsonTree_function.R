###############################################################################
# csv2jsonTree_function.R, collection of function to convert csv table to json structure
# Author: Martial M Sankar
# copyright: 2015-2016, SIB Swiss Institute of Bioinformatics
###############################################################################


library(jsonlite)


# csv2jsonTree, function that return json file and export json file
#             @fn, String, path and file name of source input file
#             @outfn, String, path and file name of export JSON file
#             @sub, INT, whether to use subset of the source csv file
csv2jsonTree  <- function(fn, outfn = NULL , sub=NULL) {
  
  # import source csv 
  
  if(!is.null(sub)){    
    tab <- read.csv(fn,sep =";", as.is=T)[1:sub,]
  }else {    
    tab <- read.csv(fn,sep =";", as.is=T)
  }
  
  # format table
  
  naix <- which(is.na(tab$level)) # NA rm
  if(length(naix>0)){tab <- tab[-naix,]}
  
  
  options(stringsAsFactors = FALSE)
  tab <- data.frame(cbind(cbind(tab, rep(FALSE, nrow(tab)), rep("undefined", nrow(tab)) ))) # add done and count column 
  names(tab) <- c("name", "level", "label", "annotationtype", "mode", "figureground", "instructions", "done", "count")
  
  # iterate to create data structure
  maxV <- max(tab$level)
  #maxl<- 2
  cat("# init max level", maxV, "\n")
  
  
  
  # select subset max level then iterate until level 0    
  tabCur <- tab[which(tab$level==maxV),]
  
  # split subset by name to get children
  tabCurL <- split(tabCur, gsub("\\.\\d+$", "", tabCur$name))
  
  # convert
  
  for (j in 1:maxV){
    
    
    # go to max level -1
    maxl <- maxV - j
    
    cat("# cur max level", maxl, "\n")
    
    
    
    tmpCur <- tab[which(tab$level == maxl),]
    if(maxl==0){
      if(length(tmpCurL)==1){
        tmpCurL <- list(unbox(tmpCur))
      }else {
        tmpCurL <- lapply(1:nrow(tmpCur), function(i){list(unbox(tmpCur[i,]))})
      }
    }else {
      
      #tmpCurL <- split( tmpCur, gsub("\\.\\d$", "", as.character(tmpCur$name)))
      tmpCurL <- split( tmpCur, as.character(tmpCur$name))
      
      checkDuplicates <- unique(unlist(lapply(tmpCurL,nrow)))
      if(length(checkDuplicates)>1){
        stop("# ERROR : Duplicates found at lines ", paste(as.character(which(unlist(lapply(tmpCurL,nrow))!=1)), collapse=" , "))
      }
      tmpCurL <- lapply( tmpCurL, unbox)
      
      
      
    }
    
    
    
    ll <- list()
    
    
    # add children if annotation type expend  
    expix <- which(tmpCur$annotationtype == "expand")
    
    
    if(length(expix)>0 ) {
      
      # gather children in list
      for ( xx in 1:length(tabCurL)){
        
        dfCur <- tabCurL[[xx]]
        nm <- tabCurL[[xx]]$name[1]
        if(is.null(nm)){ next } # if nm = NULL continue = no children for this element
        
        if(maxl==0) {
          gs  <- gsub("\\d+$", "", as.character(nm))
          
        }else {
          gs <- gsub("\\.\\d+$", "", as.character(nm))
        }
        
        #cat("# gather children for parent : ", gs, "\n")
        ll[[gs]][[nm]] <- dfCur
      }
      
      
      
      
      
      # add children to parent node 
      if(maxl==0) {
        #cat("## Added Children to ",, "\n")
        if(length(tmpCurL)==1){
          names(ll[[1]]) <- NULL
          tmpCurL[[1]][[1]]$children <- ll
          
        }else {
          for (i in 1:length(tmpCurL)){
            
            orepl <- ll[[tmpCurL[[i]][[1]]$name]]
            
            # reorder ll
            if (length(orepl)>9){
              ix <- 1:length(orepl)
                            
              valV <- unlist(lapply(names(orepl), function(nml){ 
                if(length(grep("\\.", nml))>0){
                  snml <- strsplit(nml, "\\.")[[1]]
                }else {
                  snml <- gsub("^[A-Z]", "", nml)
                }
                val <- as.numeric(snml[length(snml)])
                return(val)
              }))
              
              ix <- ix[order(valV)]                
              orepl <- orepl[ix]
            }
            
            names(orepl) <- NULL
            tmpCurL[[i]][[1]]$children <- list(orepl)
          }
          
          
          
        }
        
        
      }else {
        
        
       
        for (xl in names(ll)){ 
          cat("## Added Children to ",xl, "\n")
          
          # order names(ll) if length ll>9
          
          if(length(ll[[xl]])>9){
            ix <- 1:length(ll[[xl]])
            
            
            valV <- unlist(lapply(names(ll[[xl]]), function(nml){ 
              if(length(grep("\\.", nml))>0){
                snml <- strsplit(nml, "\\.")[[1]]
              }else {
                snml <- gsub("^[A-Z]", "", nml)
              }
              val <- as.numeric(snml[length(snml)])
              return(val)
            }))
            
            ix <- ix[order(valV)]  
            
            orepl <- ll[[xl]][ix]
            
          }else {
            orepl <- ll[[xl]]
          }
          
          names(orepl) <- NULL
          
          
          if(maxl==(maxV-1)){
            #tmpCurL[[xl]]$children <- lapply(1:nrow(orepl), function(i){list(unbox(orepl[i,]))})
            tmpCurL[[xl]]$children <- orepl
            
            
          }else {
            tmpCurL[[xl]]$children <- list(orepl)
          }
        }
      } 
    }
    
    names(tmpCurL)  <- NULL
    tabCurL <- tmpCurL
  }
  
  
  # export 
  names(tmpCurL)  <- NULL # remove object name
  objJSON <- toJSON(lapply(tmpCurL, function(oo){ unbox(oo[[1]])}), pretty=TRUE)
  
  
  #   if(length(tmpCurL)>1){
  #    objJSON <- toJSON(lapply(tmpCurL, function(oo){ unbox(oo[[1]])}), pretty=TRUE)
  #   }else {
  #     objJSON <- toJSON(tmpCurL, pretty=TRUE)   
  #   }
  
  if(!is.null(outfn)){
    write(objJSON, outfn)
    cat(outfn, "\n")
    
  }
  return(objJSON)
}