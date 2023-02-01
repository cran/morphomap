#' morphomapVariations
#'
#' Calculate cortical map variation from PCA
#' @param PCA list: list containing morphomapShape objects
#' @param scores list: list containing morphomapShape objects
#' @param PC list: list containing morphomapShape objects
#' @param pal list: list containing morphomapShape objects
#' @param asp numeric: aspect ratio of the morphometric map
#' @return mapvar: matrix containing values of cortical thickness
#' @author Antonio Profico
#' @examples
#' \donttest{
#' data(Ex_mpShapeList)
#' PCA<-morphomapPCA(Ex_mpShapeList)
#' plot(PCA$PCscores)
#' barplot(PCA$Variance[,2])
#' morphomapVariations(PCA,min(PCA$PCscores[,1]),PCA$PCs[,1])
#' morphomapVariations(PCA,max(PCA$PCscores[,1]),PCA$PCs[,1])
#' }
#' @export

morphomapVariations<-function(PCA,scores, PC,pal=blue2green2red(101),asp=2){
  mapvar<-matrix(restoreFromPCA(scores,PC,PCA$meanMap),nrow=dim(PCA$CorMaps)[2],ncol=dim(PCA$CorMaps)[1])
  map<-levelplot(mapvar,asp=asp,col.regions =pal)
  graphics::plot(map)
  return(mapvar)
}
