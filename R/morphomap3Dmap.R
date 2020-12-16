#' morphomap3Dmap
#' 
#' Plot a 3D thickness map in four different anatomical views
#' @param morphomap.shape list: output from morphomapShape function
#' @param out.sur 3D mesh: 3D mesh of the long bone
#' @param method character: if set on "equiangular" the cortical thickness is meant as the distance of the segment intersecting the external and internal outline starting from the centroid of the section. If set on "closest" the cortical thickness is calculated at each point as the closest distance between external and internal outlines
#' @param scale logical: if TRUE the cortical thickness matrix will be scaled from 0 to 1
#' @param k integer: neighbourhood of kd-tree to search the nearest semilandmarks to each vertex
#' @param plot logical: if TRUE the 3D map is plotted
#' @param rem.out logical: if TRUE outliers are identified and removed from thickness matrix
#' @param fac.out numeric: parameter to set the threshold in outliers detection
#' @param smooth logical: if TRUE the smoothing filter is applied on the thickness matrix
#' @param smooth.iter numeric: number of smoothing iterations
#' @param pal character vector: colors to be used in the map production
#' @return cols color associated at each vertex of 3D mesh
#' @return thickmat thickness matrix after smoothing and outliers removal
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' if(interactive()){
#' #morphomap on a human femur bone
#' data(HomFem38023)
#' meshes<-morphomapSegm(HomFem38023, param1=4)
#' perMesh<-meshes$external
#' endMesh<-meshes$internal
#' mech_length<-380.23
#' rawSections<-morphomapCore(out.sur=perMesh,
#' inn.sur=endMesh,num.sect=61,
#' mech.len = mech_length,param1 = 0.5,
#' radius.fact = 2.5,npovs = 100,clean_int_out = TRUE,
#' num.points = 500, start = 0.2,end=0.8)
#' shapeSections<-morphomapShape(rawSections,24,sects_vector=NULL,cent.out="CCA",
#' delta=0.1, side="left")
#' 
#' #built 3D morphometric map 
#' bone3Dmap<-morphomap3Dmap(shapeSections, out.sur=perMesh,
#'                           plot = TRUE,rem.out=TRUE,
#'                           fac.out=1.5,smooth=TRUE,
#'                           smooth.iter=5)
#' #or
#' require(rgl)
#' rgl::open3d()
#' rgl::shade3d(perMesh,col=bone3Dmap$cols,specular="black")
#' 
#' #morphomap on a chimpanzee femur bone
#' data(PanFem27713)
#' meshes<-morphomapSegm(PanFem27713, param1=3)
#' perMesh<-meshes$external
#' endMesh<-meshes$internal
#' mech_length<-277.13
#' rawSections<-morphomapCore(out.sur=perMesh,
#'                            inn.sur=endMesh,num.sect=61,mech.len = mech_length, 
#'                            start = 0.2,end=0.8)
#' shapeSections<-morphomapShape(rawSections,24,sects_vector=NULL,cent.out="CCA",
#' delta=0.1, side="left")
#' #built 3D morphometric map 
#' bone3Dmap<-morphomap3Dmap(shapeSections, out.sur=perMesh,
#'                           plot = TRUE,rem.out=TRUE,
#'                           fac.out=1.5,smooth=TRUE,
#'                           smooth.iter=5)
#' #or
#' require(rgl)
#' rgl::open3d()
#' rgl::shade3d(perMesh,col=bone3Dmap$cols,specular="black")
#' }
#' @export

morphomap3Dmap<-function (morphomap.shape,out.sur,method = "equiangular",
                          scale = TRUE,rem.out = FALSE,fac.out = 0.5,
                          smooth = FALSE,smooth.iter = 5, k=5,
                          plot = TRUE,pal = blue2green2red(101)) 
{
  morphomap.thickness<-morphomapThickness(morphomap.shape)
  morphomap.shape<-morphomap.thickness$morphomap.shape
  
  if (method == "equiangular") {
    thicks <- morphomap.thickness$sect_thickness_eq
  }
  if (method == "closest") {
    thicks <- morphomap.thickness$sect_thickness_cp
  }
  X <- NULL
  Y <- NULL
  Z <- NULL
  summary <- dim(thicks)
  num.sect <- summary[3]
  num.land <- summary[1]
  
  if (scale == TRUE) {
    x <- as.vector(thicks)
    thickvector <- (x - min(x))/(max(x) - min(x))
    Thick_arr <- array(thickvector, dim = c(num.land, 1, num.sect))
    thicks <- Thick_arr
  }
  
  if (rem.out == TRUE) {
    x <- as.vector(thicks)
    qnt <- stats::quantile(x, probs = c(0.25, 0.75))
    H <- fac.out * IQR(x)
    y <- x
    y[x < (qnt[1] - H)] <- min(y[x > (qnt[1] - H)])
    y[x > (qnt[2] + H)] <- max(y[x < (qnt[2] + H)])
    thick_values <- y
    Thick_arr <- array(thick_values, dim = c(num.land, 1, num.sect))
    thicks <- Thick_arr
  }
  
  if (smooth == TRUE) {
    thick_smoo <- apply(thicks, 2, function(x) (matrixSmooth(x, 
                                                             passes = smooth.iter)))
    thick_out <- array(thick_smoo, dim = c(num.land, 1, num.sect))
    thicks <- thick_out
  }
  
  for (i in 1:num.sect) {
    X_i <- seq(1, 100, length.out = num.land)
    Y_i <- morphomap.shape$`3D_out`[, 3, i]
    Z_i <- thicks[, , i]
    X <- c(X, X_i)
    Y <- c(Y, Y_i)
    Z <- c(Z, Z_i)
  }  
  
  mech.len <- morphomap.shape$mech.len
  mat <- morphomap.shape$`3D_out`[, , 1]
  for (i in 2:dim(morphomap.shape$`3D_out`)[3]) {
    mat <- rbind(mat, morphomap.shape$`3D_out`[, , i])
  }
  pmeshdist <- vcgClostKD(mat, out.sur)
  
  
  thick_values <- as.vector(unlist(thicks))
  cols <- pal
  intervals <- seq(min(thick_values), max(thick_values), length = length(cols))
  pos_cols <- NULL
  for (i in 1:length(thick_values)) {
    pos_cols[i] <- which.min(abs(thick_values[i] - intervals))
  }
  clostInd <- mcNNindex(mat, t(out.sur$vb)[, 1:3], k = k)
  distInd <- clostInd
  for (i in 1:ncol(clostInd)) {
    distInd[, i] <- sqrt(rowSums((t(out.sur$vb)[, 1:3] - 
                                    mat[clostInd[, k], ])^2))
  }
  colsInd <- NULL
  for (i in 1:dim(clostInd)[1]) {
    colsInd[i] <- round(weighted.mean(pos_cols[clostInd[i, 
                                                        ]], distInd[i, ]/sum(distInd[i, ])))
  }
  white_col <- c(which(t(out.sur$vb)[, 3] < mech.len * morphomap.shape$start), 
                 which(t(out.sur$vb)[, 3] > mech.len * morphomap.shape$end))
  cols <- c(cols, "#FFFFFF")
  colsInd[white_col] <- length(cols)
  col <- cols[colsInd]
  it <- as.vector(out.sur$it)
  colls <- col[it]
  if (plot == TRUE) {
    mat <- matrix(1:4, ncol = 2, byrow = TRUE)
    open3d()
    layout3d(mat, sharedMouse = TRUE)
    triangles3d(t(out.sur$vb[1:3, out.sur$it]), col = colls, 
                specular = "black")
    next3d()
    out.sur_r <- rgl::rotate3d(out.sur, pi/2, 0, 0, 1)
    triangles3d(t(out.sur_r$vb[1:3, out.sur_r$it]), col = colls, 
                specular = "black")
    next3d()
    out.sur_r <- rgl::rotate3d(out.sur, pi, 0, 0, 1)
    triangles3d(t(out.sur_r$vb[1:3, out.sur_r$it]), col = colls, 
                specular = "black")
    next3d()
    out.sur_r <- rgl::rotate3d(out.sur, pi * 3/2, 0, 0, 1)
    triangles3d(t(out.sur_r$vb[1:3, out.sur_r$it]), col = colls, 
                specular = "black")
  }
  out <- list(cols = col, thickmat = thick_values)
  return(out)
}