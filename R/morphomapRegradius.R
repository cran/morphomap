#' morphomapRegradius
#' 
#' Wrapper of the function regularradius written by Julien Claude (Morphometrics with R)
#' @param mat a kx2 matrix 
#' @param center coordinates of the center from which the calculation of regular radius started 
#' @param n number of points
#' @return V2 position of landmarks equi angular spaced
#' @author Julien Claude, Antonio Profico
#' @references Claude, J. (2008). Morphometrics with R. Springer Science & Business Media.
#' @examples
#' extsec<-morphomapCircle(10,1000)
#' sel<-morphomapRegradius(extsec,center = c(0,0),n=11)
#' selcoo<-extsec[sel,]
#' plot(extsec,type="l",asp=1)
#' points(selcoo,col="red",pch=19)
#' @export

morphomapRegradius<-function(mat, center, n){
  Rx<-mat[,1]
  Ry<-mat[,2]
  le<-length(Rx)
  M<-matrix(c(Rx, Ry), le, 2)
  M1<-matrix(c(Rx-center[1], Ry-center[2]), le, 2)
  V1<-complex(real=M1[,1], imaginary=M1[,2])
  M2<-matrix(c(Arg(V1), Mod(V1)), le, 2)
  V2<-NA
  for (i in 0:(n-1))
  {V2[i+1]<-which.max((cos(M2[,1]-2*i*pi/n)))}
  return(V2)
}
