#' morphomapZmoment
#' 
#' Calculate the polar moment of inertia around the x and y axes and the polar section module
#' @param cp matrix: coordinates of the external outline of the section 
#' @param mp matrix: coordinates of the internal outline of the section 
#' @param Cx numeric: x coordinate of the section center
#' @param Cy numeric: y coordinate of the section center
#' @param delta numeric: picture elements of adjustable side length 
#' @return Zx numeric: moment of inertia around the x axis
#' @return Zy numeric: moment of inertia around the y axis
#' @return dx numeric: maximum chord length from y axis
#' @return dy numeric: maximum chord length from x axis
#' @return Zpol numeric: polar moment of inertia 
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' extsec<-morphomapCircle(10,1000)
#' intsec<-morphomapCircle(8,1000)
#' ZMs<-morphomapZmoment(extsec,intsec,delta=0.1)
#' @export

morphomapZmoment<-function(cp,mp,Cx=0,Cy=0,delta=0.1){
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
  
  Zpol_t<-(sqrt(rowSums((grid_sect[sel,]-c(Cx,Cy))^2)))
  maxr<-max(Zpol_t)
  Zpol<-sum((Zpol_t^2)*delta2)/maxr
  
  sumx<-sum((grid_sect[sel,1]*grid_sect[sel,1])*delta2)
  sumy<-sum((grid_sect[sel,2]*grid_sect[sel,2])*delta2)
  
  dx<-max(abs(range(cp[,1])))
  dy<-max(abs(range(cp[,2])))
  
  Zx<-Ix/dy
  Zy<-Iy/dx
  
  Zms<-c(Zx,Zy,dx,dy,Zpol)
  names(Zms)<-c("Zx","Zy","dx","dy","Zpol")
  return(Zms)
}
