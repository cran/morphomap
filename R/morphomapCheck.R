#' morphomapCheck
#'
#' Plot the long bone mesh to check the orientation of the long bone
#' @param mesh 3D mesh: long bone 3D model
#' @param col character: color mesh
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' \donttest{
#' library(morphomap)
#' data(HomFem38023)
#' morphomapCheck(HomFem38023)
#' }
#' @export

morphomapCheck<-function(mesh,col="white"){
  open3d()
  wire3d(mesh,col=col)
  axis3d('x', pos = c(NA, 0, 0),lwd=5,col="red")
  axis3d('y', pos = c(0, NA, 0),lwd=5,col="green")
  axis3d('z', pos = c(0, 0, NA),lwd=5,col="blue")
  title3d(main=NULL,xlab="X axis",zlab="Z (Biomechanical length)")
  bbox3d()
  view3d(userMatrix=rotationMatrix(pi, 1, 0, 1),fov=0)
  U <- par3d("userMatrix")
  par3d(userMatrix = rotate3d(U, pi, 0,0,1))

  p1 <- c(0, 0, 0)
  p2 <- c(0, 0, 100)
  p3 <- c(100, 0, 100)
  normal <- crossProduct(p2-p1,p3-p1)
  planes3d(normal[1],normal[2],normal[3],d=0,alpha=0.5,col="violet")

  p1 <- c(0, 0, 0)
  p2 <- c(100, 0, 0)
  p3 <- c(0, 100, 0)
  normal <- crossProduct(p2-p1,p3-p1)
  planes3d(normal[1],normal[2],normal[3],d=0,alpha=0.5,col="yellow")

  p1 <- c(0, 0, 0)
  p2 <- c(0, 0, 100)
  p3 <- c(0, 100, 100)
  normal <- crossProduct(p2-p1,p3-p1)
  planes3d(normal[1],normal[2],normal[3],d=0,alpha=0.5,col="lightblue")

}
