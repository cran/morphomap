#' morphomapCentroid
#' 
#' Calculate the barycenter of the cortical area
#' @param cp matrix: coordinates of the external outline of the section 
#' @param mp matrix: coordinates of the internal outline of the section 
#' @param delta numeric: picture elements of adjustable side length
#' @return centroid numeric vector: coordinates of the cortical area
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' extsec<-morphomapCircle(10,100)
#' intsec<-morphomapCircle(8,100)
#' plot(extsec,asp=1,type="l")
#' points(intsec,col=2,type="l")
#' cent<-morphomapCentroid(extsec,intsec,delta = 0.1)  
#' points(cent[1],cent[2],pch=19,col=3)
#' @export

morphomapCentroid<-function(cp,mp,delta=0.1){
  maxx<-max(cp[,1],na.rm=T)
  minx<-min(cp[,1],na.rm=T)
  maxy<-max(cp[,2],na.rm=T)
  miny<-min(cp[,2],na.rm=T)
  X<-seq(minx+delta/2,maxx-delta/2,delta)
  Y<-seq(miny+delta/2,maxy-delta/2,delta)
  grid_sect<-as.matrix(expand.grid(X, Y))
  A<-point.in.polygon(grid_sect[,1], grid_sect[,2], cp[,1], cp[,2], mode.checked=FALSE)
  B<-point.in.polygon(grid_sect[,1], grid_sect[,2], mp[,1], mp[,2], mode.checked=FALSE)
  sel<-which(A==1 & B==0)
  centroid<-colMeans(grid_sect[sel,])
  return(centroid)
}