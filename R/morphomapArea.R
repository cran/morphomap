#' morphomapArea
#' 
#' Shoelace formula to calculate the area of a closed outline
#' @param p matrix: kx2 matrix 
#' @param delta numeric: picture elements of adjustable side length 
#' @param method character: the user can choice to calculate the area applying the "shoelace" formula or discretizing the cross sections in dA areas (method = "delta") 
#' @return ar numeric: area
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' extsec<-morphomapCircle(10,100)
#' area<-morphomapArea(extsec, method="shoelace")
#' @export

morphomapArea<-function (p,delta=0.1,method="shoelace") 
{
  if (method=="shoelace"){
    k <- length(p[, 1])
    p <- rbind(p, p[1, ])
    ar <- 0
    for (i in 1:k) {
      ar <- ar + (p[i, 1] * p[i + 1, 2]) - (p[i + 1, 1] * p[i, 
                                                            2])
    }
    ar <- abs(ar/2)
  }
  
  if (method=="delta"){
    maxx<-max(p[,1],na.rm=T)
    minx<-min(p[,1],na.rm=T)
    maxy<-max(p[,2],na.rm=T)
    miny<-min(p[,2],na.rm=T)
    X<-seq(minx+delta/2,maxx-delta/2,delta)
    Y<-seq(miny+delta/2,maxy-delta/2,delta)
    grid_sect<-as.matrix(expand.grid(X, Y))
    A<-point.in.polygon(grid_sect[,1], grid_sect[,2], p[,1], p[,2], mode.checked=FALSE)
    sel<-which(A==1)
    delta2=delta*delta  
    ar<-length(sel)*delta2
  }
  return(ar)
}