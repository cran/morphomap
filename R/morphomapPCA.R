#' morphomapPCA
#'
#' Calculate maps of cortical thickness and perform a Principal Component Analysis
#' @param mpShapeList list: list containing morphomapShape objects
#' @param gamMap list: list containing morphomapShape objects
#' @param nrow list: list containing morphomapShape objects
#' @param ncol list: list containing morphomapShape objects
#' @param rem.out list: list containing morphomapShape objects
#' @param scaleThick list: list containing morphomapShape objects
#' @param fac.out list: list containing morphomapShape objects
#' @param method list: list containing morphomapShape objects
#' @param scalePCA list: list containing morphomapShape objects
#' @param unwrap list: list containing morphomapShape objects
#' @return PCscores PC scores
#' @return PCs loadings
#' @return Variance Table of the explained Variance by the PCs
#' @return meanMap mean morphometric map
#' @return CorMaps morphometric maps
#' @author Antonio Profico
#' @examples
#' \donttest{
#' data(Ex_mpShapeList)
#' PCA<-morphomapPCA(Ex_mpShapeList)
#' plot(PCA$PCscores)
#' barplot(PCA$Variance[,2])
#' }
#' @export

morphomapPCA<-function(mpShapeList,gamMap = TRUE,nrow = 61,ncol = 24,
                       rem.out = TRUE,scaleThick = FALSE,fac.out = 1.5,
                       method="equiangular",scalePCA=TRUE,unwrap = "A"){
  
  dims<-dim(mpShapeList[[1]]$`3D_out`)
  nland<-dims[1]
  nrows<-dims[3]
  if(length(unique(unlist(lapply(mpShapeList, function(x) dim(x$`3D_out`)[1])))) != 1){
    stop("the number of landmarks is different")
  }
  if(length(unique(unlist(lapply(mpShapeList, function(x) dim(x$`3D_out`)[1]))))!=1){
    stop("the number of cross-sections is different")
  }
  if(isFALSE(gamMap)){nrow<-nrows; ncol<-nland}
  ThickMat<-array(NA,dim=c(nrow,ncol,length(mpShapeList)))
  for(i in 1:dim(ThickMat)[3]){
    map2Di<-morphomap2Dmap(mpShapeList[[i]],rem.out = rem.out,scale = scaleThick,fac.out = fac.out,gamMap = gamMap,nrow = nrow,ncol = ncol,
                           plot=FALSE,unwrap = unwrap )
    ThickMat[,,i]<-matrix(map2Di$data[,3],nrow = nrow,ncol=ncol)
  }
  
  CortThick_mat<-vecx(ThickMat)
  Thick_PCA<-prcomp(CortThick_mat,scale. = scalePCA)
  values <- 0
  eigv <- Thick_PCA$sdev^2
  values <- eigv[which(eigv > 1e-16)]
  lv <- length(values)
  PCs <- Thick_PCA$rotation[, 1:lv]
  PCscores <- as.matrix(Thick_PCA$x[, 1:lv])
  Variance <- cbind(sqrt(eigv), eigv/sum(eigv), cumsum(eigv)/sum(eigv)) * 100
  Variance <- Variance[1:lv, ]
  
  out<-list("PCscores"=PCscores,"PCs"=PCs,"Variance"=Variance,"meanMap"=Thick_PCA$center,"CorMaps"=ThickMat)
  return(out)
}