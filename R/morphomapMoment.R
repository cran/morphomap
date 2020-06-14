#' morphomapMoment
#' 
#' Calculate the moment of inertia around the x and y axes and the product of inertia
#' @param cp matrix: coordinates of the external outline
#' @param mp matrix: coordinates of the internal outline 
#' @param delta numeric: picture elements of adjustable side length
#' @return Ix numeric: moment of inertia around the x axis
#' @return Iy numeric: moment of inertia around the y axis
#' @return Ixy numeric: product of inertia around the x and y axis
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' #create a section
#' extsec<-morphomapCircle(10,1000)
#' intsec<-morphomapCircle(8,1000)
#' InMs<-morphomapMoment(extsec,intsec,delta=0.1)
#' @export

morphomapMoment<-function(cp,mp,delta=0.1){
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
  delta2=delta*delta
  Ix<-sum((grid_sect[sel,2]^2)*delta2)
  Iy<-sum((grid_sect[sel,1]^2)*delta2)  
  Ixy<-sum((grid_sect[sel,1]*grid_sect[sel,2])*delta2)
  names(Ix)<-"Ix"
  names(Iy)<-"Iy"
  names(Ixy)<-"Ixy"
  return(c(Ix,Iy,Ixy))
}
