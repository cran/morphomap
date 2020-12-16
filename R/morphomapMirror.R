#' morphomapMirror
#'
#' Mirror a long bone mesh along the yz plane
#' @param mesh object of class mesh3d
#' @return mesh: object of class mesh3d
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' if(interactive()){
#' #a left human femur bone
#' require(rgl)
#' data(HomFem38023)
#' lfem<-HomFem38023
#' rfem<-morphomapMirror(lfem)
#' rgl::open3d()
#' rgl::wire3d(lfem,col="green")
#' rgl::ire3d(rfem,col="red")
#' }
#' @export

morphomapMirror<-function(mesh){
  p1 <- c(0, 0, 0)
  p2 <- c(0, 0, 100)
  p3 <- c(0, 100, 100)
  normal <- crossProduct(p2-p1,p3-p1)
  mesh<-mirror2plane(mesh, p1, normal = normal, v2 = p2, v3 = p3)
  return(mesh)
}
