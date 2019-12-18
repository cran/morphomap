#' morphomapCircle
#' 
#' Define a circular outline
#' @param r numeric: radius of the outline 
#' @param n numeric: number of points along the outline
#' @return mat matrix with coordinates
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' extsec<-morphomapCircle(10,100)
#' intsec<-morphomapCircle(8,100)
#' plot(extsec,asp=1,type="l")
#' points(intsec,type="l",col=2)
#' @export

morphomapCircle<-function(r=1,n=1000){
  X<-NULL
  Y<-NULL
  seqs<-seq(0,2*pi,length.out = n)
  for(i in 1:length(seqs)){
    X<-c(X,r*sin(seqs[i]))
    Y<-c(Y,r*cos(seqs[i]))
  }
  mat<-cbind(X,Y)
}

