#'Convert camel-case string to underscore separated
#'
#'@param x Character vector to be converted to camel-case
#'@return A Character vector with underscore separated words
#'@export
fillUnderscore<-function(x){
  under<-function(x){
    whichCaps<-which(unlist(strsplit(x,""))!=unlist(strsplit(tolower(x),"")))
    if(1 %in% whichCaps){whichCaps<-whichCaps[-1]}
    xSplit<-substring(x,c(1,whichCaps),c(whichCaps-1,nchar(x)))
    filled<-paste(tolower(xSplit),collapse="_")
    return(filled)
  }
  return(unlist(lapply(x,under)))
}

#'Convert string of separated words to camel-case
#'
#'@param x Character vector to be converted to camel-case
#'@param sep Character separating words in the original string
#'@return Camel-case string(s)
#'@export
camelCase<-function(x,sep="_"){
  camel<-function(toCamel,sep){
    xSplit<-unlist(strsplit(toCamel,sep))
    if(length(xSplit)==1){
      xSplit<-paste0(tolower(substr(xSplit,1,1)),substr(xSplit,2,nchar(xSplit)))
      return(xSplit)}
    xSplit[1]<-tolower(xSplit[1])
    for(n in 2:length(xSplit)){
      xSplit[n]<-paste0(toupper(substr(xSplit[n],1,1)),
                        substr(xSplit[n],2,nchar(xSplit[n])))
    }
    camelled<-paste(xSplit,collapse="")
    return(camelled)
  }
  return(unlist(lapply(x,camel,sep=sep)))
}
