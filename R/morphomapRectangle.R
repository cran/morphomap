#' morphomapRectangle
#' 
#' Define a rectangular outline
#' @param l numeric: length of the rectangle
#' @param h numeric: height of the rectangle
#' @param n numeric: number of points along the outline
#' @return mat matrix with coordinates
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' extsec<-morphomapRectangle(10,6,100)
#' intsec<-morphomapRectangle(8,4,100)
#' plot(extsec,asp=1,type="l")
#' points(intsec,type="l",col=2)
#' @export

morphomapRectangle<-function(l=1,h=1,n=1000){
  
  xL<-seq(-l,l,length.out = n)
  yL1<-rep(h,n)
  yL2<-rep(-h,n)
  
  yH<-seq(-h,h,length.out = n)
  xH1<-rep(-l,n)
  xH2<-rep(l,n)
  
  mat<-rbind(cbind(xL,yL1),
             cbind(xH2,rev(yH)),
             cbind(rev(xL),yL2),
             cbind(xH1,yH))
}

