#' morphomapArray2matrix
#' 
#' Convert an array into a matrix
#' @param array an array
#' @return mat a matrix
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @export

morphomapArray2matrix<-function(array){
  mat <- matrix(aperm(array, c(1, 3, 2)), nrow=dim(array)[1]* dim(array)[3], ncol=dim(array)[2] )
  return(mat)
}