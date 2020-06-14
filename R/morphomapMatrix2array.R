#' morphomapMatrix2array
#' 
#' Convert a matrix into an array 
#' @param matrix an array
#' @param nsects number of cross sections
#' @return array an array
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @export

morphomapMatrix2array<-function(matrix,nsects){
  nland<-dim(matrix)[1]/(nsects)
  k<-dim(matrix)[2]
  start<-seq(1,nsects*nland,nland)
  end<-seq(nland,nsects*nland,nland)
  array<-array(NA,dim=c(nland,k,nsects))
  for(i in 1:nsects){
    array[,,i]<-matrix[c(start[i]:end[i]),]
  }
  return(array)
}