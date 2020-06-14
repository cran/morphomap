#' morphomapWriteMorphologika
#' 
#' Export an array in the morphologika format file
#' @param array an array
#' @param groups a vector containing a classifier
#' @param variables list containing further classifiers
#' @param file path of the file to be saved
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @export

morphomapWriteMorphologika<-function (array, groups = NULL, variables = NULL, file) 
{
  out<-list("array"=NULL,"groups"=NULL,"variables"=NULL)
  nind <- dim(array)[3]
  nlan <- dim(array)[1]
  ndim <- dim(array)[2]
  if(is.na(nind)==TRUE){
    nind<- 1
    array<-array(array,dim=c(nlan,ndim,nind))
  } 
  
  speclabels <- dimnames(array)[[3]]
  if (is.null(speclabels) == TRUE) {
    speclabels <- paste("specimen", 1:nind)
  }
  if (is.null(groups) == FALSE) {
    levgroups <- levels(as.factor(groups))
    sizgroups <- as.numeric(as.factor(groups))
  }
  if (is.null(groups) == FALSE) {
    newgroupv <- NULL
    newspecsp <- NULL
    sizegroup <- NULL
    for (i in 1:length(levgroups)) {
      group_i <- (which(as.character(groups == levgroups[i]) == 
                          TRUE))
      l1 <- length(group_i)
      newgroupv <- c(newgroupv, rep(levgroups[i], l1))
      newspecsp <- c(newspecsp, group_i)
      sizegroup <- c(sizegroup, l1)
    }
    array <- array[, , newspecsp]
    speclabels <- speclabels[newspecsp]
    if (is.null(variables) == FALSE) 
      for (i in 1:length(variables)) {
        variables[[i]] <- variables[[i]][newspecsp]
      }
  }
  out$array <- array
  out$labels <- speclabels
  if (is.null(groups) == FALSE) {
    out$groups <- newgroupv
  }
  if (is.null(variables) == FALSE) {
    out$variables <- variables
  }
  cat("[individuals]\n", file = file, sep = "")
  cat(paste(nind, "\n", sep = ""), file = file, 
      append = TRUE, sep = "")
  cat("[landmarks]\n", file = file, append = TRUE, sep = "")
  cat(paste(nlan, "\n", sep = ""), file = file, 
      append = TRUE, sep = "")
  cat("[dimensions]\n", file = file, append = TRUE, sep = "")
  cat(paste(ndim, "\n", sep = ""), file = file, 
      append = TRUE, sep = "")
  if (is.null(groups) == FALSE) {
    cat("[groups]\n", file = file, append = TRUE, sep = "")
    cat(paste(levgroups, sizegroup), file = file, append = TRUE, 
        collapse = "")
    cat(paste("", "\n", sep = ""), file = file, 
        append = TRUE, sep = "")
  }
  cat("[names]\n", file = file, append = TRUE, sep = "")
  cat(paste(speclabels, "\n", sep = ""), file = file, 
      append = TRUE, sep = "")
  if (is.null(variables) == FALSE) {
    cat("[labels]\n", file = file, append = TRUE, sep = "")
    cat(paste(names(variables), collapse = " "), file = file, 
        append = TRUE, sep = "")
    cat(paste("", "\n", sep = ""), file = file, 
        append = TRUE, sep = "")
    cat("[labelvalues]\n", file = file, append = TRUE, 
        sep = "")
    for (i in 1:nind) {
      cat(paste((unlist(do.call(Map, c(c, variables))[[i]])), 
                collapse = " "), file = file, append = TRUE, 
          sep = "")
      cat(paste("", "\n", sep = ""), 
          file = file, append = TRUE, sep = "")
    }
  }
  cat("[rawpoints]\n", file = file, append = TRUE, sep = "")
  for (i in 1:nind) {
    cat(paste("' ", i, "\n", sep = ""), 
        file = file, append = TRUE, sep = "")
    write.table(format(array[, , i], scientific = F, trim = T), 
                file = file, sep = " ", append = TRUE, quote = FALSE, 
                row.names = FALSE, col.names = FALSE, na = "", 
                eol = "\n")
  }
}