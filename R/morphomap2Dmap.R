#' morphomap2Dmap
#' 
#' Create a 2D cortical thickness map
#' @param morphomap.shape list: output from morphomapShape function
#' @param rem.out logical: if TRUE the outlier will be removed 
#' @param fac.out numeric: parameter to set the threshold in outliers detection
#' @param smooth logical: if TRUE a smooth filter is applied
#' @param scale logical: if TRUE the thichkness matrix is scaled from 0 to 1 
#' @param smooth.iter numeric: number of smoothing iterations 
#' @param gamMap logical: if TRUE gam smoothing is applied
#' @param nrow numeric: number of rows for gam smoothing
#' @param ncol numeric: number of columns for gam smoothing
#' @param gdl numeric: number of degree of freedom for gam smoothing
#' @param method character: if set on "equiangular" the cortical thickness is meant as the distance of the segment intersecting the external and internal outline starting from the centroid of the section. If set on "closest" the cortical thickness is calculated at each point as the closest distance between external and internal outlines
#' @param unwrap character: starting qaudrant to unwrap the diaphysis ("A"=anterior, "L"=lateral, "P"=posterior, "M"=mesial)
#' @param plot logical: if TRUE the 2D morphometric map is plotted
#' @param pal character vector: colors to be used in the map production
#' @param aspect numeric: axis ratio for 2D morphometric map
#' @return dataframe dataframe for colormap production
#' @return 2Dmap thickness color map
#' @return gamoutput output from GAM
#' @return data input used to build the GAM map
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' if (interactive()){
#' library(colorRamps)
#' #morphomap on a human femur bone
#' data(HomFem38023)
#' meshes<-morphomapSegm(HomFem38023, param1=4)
#' perMesh<-meshes$external
#' endMesh<-meshes$internal
#' mech_length<-380.23
#' rawSections<-morphomapCore(out.sur=perMesh,
#' inn.sur=endMesh,num.sect=61,mech.len = mech_length, start = 0.2,end=0.8)
#' shapeSections<-morphomapShape(rawSections,24,sects_vector=NULL,cent.out="CCA",
#' delta=0.1, side="left")
#' 
#' #built 2D morphometric map without GAM smoothing
#' bone2Dmap<-morphomap2Dmap(morphomap.shape=shapeSections,
#'                           plot = TRUE, rem.out = TRUE,fac.out = 1.0, pal = blue2green2red(101),
#'                           aspect=2)
#' #built 2D morphometric map with GAM smoothing
#' bone2Dmap<-morphomap2Dmap(morphomap.shape=shapeSections,gam=TRUE,
#'                           plot = TRUE, rem.out = TRUE,fac.out = 1.0, pal = blue2green2red(101),
#'                           aspect=2)
#' 
#' #morphomap on a chimpanzee femur bone
#' data(PanFem27713)
#' meshes<-morphomapSegm(PanFem27713, param1=3)
#' perMesh<-meshes$external
#' endMesh<-meshes$internal
#' mech_length<-277.13
#' rawSections<-morphomapCore(out.sur=perMesh,
#'                            inn.sur=endMesh,num.sect=61,mech.len = mech_length, start = 0.2,end=0.8)
#' shapeSections<-morphomapShape(rawSections,24,sects_vector=NULL,cent.out="CCA",
#' delta=0.1, side="left")
#' #built 2D morphometric map without GAM smoothing
#' bone2Dmap<-morphomap2Dmap(morphomap.shape=shapeSections,plot = TRUE, 
#' rem.out = TRUE,fac.out = 1.0,pal = blue2green2red(101),aspect=2)
#' #built 2D morphometric map with GAM smoothing
#' bone2Dmap<-morphomap2Dmap(morphomap.shape=shapeSections,gam=TRUE,
#'                           plot = TRUE, rem.out = TRUE,fac.out = 1.0,pal = blue2green2red(101),
#'                           aspect=2)
#' }
#' @export


morphomap2Dmap<-function (morphomap.shape,rem.out = FALSE,
                          fac.out = 0.5,smooth = FALSE,scale = TRUE,
                          smooth.iter = 5,gamMap = FALSE,nrow = 90,
                          ncol = 100,gdl = 250,method = "equiangular",
                          unwrap = "A",plot = TRUE,pal = blue2green2red(101),
                          aspect = 2) 
  {
  thickarray<-morphomapThickness(morphomap.shape)
  data_t<-morphomapDF(thickarray,rem.out = rem.out, fac.out = fac.out, smooth = smooth, scale = scale, 
                    smooth.iter = smooth.iter, method = method, unwrap = unwrap)
  data <- data_t$XYZ
  labels <- data_t$labels
  if (gamMap == FALSE) {
    minx <- min(data$X, na.rm = T)
    miny <- min(data$Y, na.rm = T)
    maxx <- max(data$X, na.rm = T)
    maxy <- max(data$Y, na.rm = T)
    maxz <- max(data$Z, na.rm = T)
    minz <- min(data$Z, na.rm = T)
    atvalues <- seq(minz, maxz, length.out = 100)
    xat <- seq(1, 100, length.out = 5)
    ylabels<-round(seq(morphomap.shape$start,morphomap.shape$end,length.out = 5)*100,2)
    yat<-round(seq(miny,maxy,length.out = 5),2)
    map <- levelplot(data[, 3] ~ data[, 1] + data[, 2], xlab = "", 
                     ylab = "Biomechanical length", aspect = aspect, main = "2D thickness map", 
                     at = atvalues, col.regions = pal, scales = list(x = list(at = xat,labels = labels, rot = 90, alternating = 1),
                                                                     y = list(at = yat,labels = ylabels, rot = 90, alternating = 1)))
    if (plot == TRUE) {
      graphics::plot(map)
    }
    mat <- data
  }
  if (gamMap == TRUE) {
    m1 <- gam(Z ~ s(X, Y, k = gdl), data = data)
    minx <- min(data$X, na.rm = T)
    maxx <- max(data$X, na.rm = T)
    miny <- min(data$Y, na.rm = T)
    maxy <- max(data$Y, na.rm = T)
    wx <- maxx - minx
    wy <- maxy - miny
    sx <- wx/ncol
    sy <- wy/nrow
    xx <- seq(minx, maxx, length.out = ncol)
    yy <- seq(miny, maxy, length.out = nrow)
    xp <- vector()
    yp <- vector()
    xc <- vector()
    yc <- vector()
    for (i in seq(1, nrow)) for (j in seq(1, ncol)) {
      xp[(i - 1) * ncol + j] <- xx[j]
      yp[(i - 1) * ncol + j] <- yy[i]
      yc[(i - 1) * ncol + j] <- yy[i]
      xc[(i - 1) * ncol + j] <- xx[j]
    }
    fit <- predict.gam(m1, list(X = xp, Y = yp))
    data <- data.frame(X = xc, Y = yc, Z = fit)
    if (scale == TRUE) {
      x <- data$Z
      thickvector <- (x - min(x))/(max(x) - min(x))
      data$Z <- thickvector
    }
    minx <- min(data$X, na.rm = T)
    miny <- min(data$Y, na.rm = T)
    maxx <- max(data$X, na.rm = T)
    maxy <- max(data$Y, na.rm = T)
    maxz <- max(data$Z, na.rm = T)
    minz <- min(data$Z, na.rm = T)
    atvalues <- seq(minz, maxz, length.out = 100)
    xat <- seq(1, 100, length.out = 5)
    ylabels<-round(seq(morphomap.shape$start,morphomap.shape$end,length.out = 5)*100,2)
    yat<-round(seq(miny,maxy,length.out = 5),2)
    map <- levelplot(data[, 3] ~ data[, 1] + data[, 2], xlab = "", 
                     ylab = "Biomechanical length", aspect = aspect, main = "2D thickness map", 
                     at = atvalues, col.regions = pal, scales = list(x = list(at = xat,labels = labels, rot = 90, alternating = 1),
                                                                     y = list(at = yat,labels = ylabels, rot = 90, alternating = 1)))
    if (plot == TRUE) {
      graphics::plot(map)
    }
    mat <- data
  }
  if(gamMap==TRUE){
  out <- list(dataframe = mat, `2Dmap` = map,gamoutput=m1,data=data)}
  if(gamMap==FALSE){
  out <- list(dataframe = mat, `2Dmap` = map,gamoutput=NULL,data=data)}

  return(out)
  
}
