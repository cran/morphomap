#' morphomapPlotShape
#' 
#' Visualize 2D and 3D cross sections 
#' @param Shape list: output from morphomapShape function
#' @param dims numeric:  2 = bi-dimensional cross sections, 3 = three-dimensional cross sections
#' @param col1 color of the external outline
#' @param col2 color of the internal outline
#' @param colc color of the centroid of the cross section 
#' @param colr color of the radii
#' @param coll1 color of the lines on the enternal outline
#' @param coll2 color of the lines on the internal outline
#' @param size numeric: points and spheres size
#' @param lwd numeric: line width in pixels
#' @param colmesh1 color of the periosteal mesh
#' @param colmesh2 color of the endosteal mesh
#' @param alpha numeric: alpha value between 0(fully transparent) and 1 (opaque)
#' @param tri logical: if TRUE the semilandmarks configuration is triangulated
#' @param outlines logical: if TRUE the 2D and 3D outlines are plotted
#' @param points logical: if TRUE points (2D) and spheres (3D) are plotted
#' @param lines logical: if TRUE 2D and 3D lines are plotted
#' @param centroid logical: if TRUE 2D and 3D centroids are plotted
#' @param cent.out how to define the center of each section. The method allowed are "CCA" (center of cortical area), "E" (barycenter of the external outline) and "I" (barycenter of the internal outline)
#' @param delta pixel size used to calculate the CCA
#' @param vecs numeric: which sections will be plotted. If dims is set on 2 only the first element of the vector vecs is considered
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' if(interactive()){
#' #morphomap on a human femur bone
#' data(HomFem38023)
#' meshes<-morphomapSegm(HomFem38023)
#' perMesh<-meshes$external
#' endMesh<-meshes$internal
#' mech_length<-380.23
#' rawSections<-morphomapCore(out.sur=perMesh,
#' inn.sur=endMesh,num.sect=61,mech.len = mech_length, start = 0.2,end=0.8)
#' shapeSections<-morphomapShape(rawSections,21,sects_vector=NULL,cent.out="CCA",
#' delta=0.1, side="left")
#' #Plot the object morphomapShape in three dimensions
#' morphomapPlotShape(shapeSections,dims=3, size=0.5)
#' #Plot a 2D cross-section 
#' morphomapPlotShape(shapeSections,dims=2,lines=TRUE,vecs=31)
#' }
#' @export

morphomapPlotShape<-function (Shape, dims = 3, col1 = "red", col2 = "green", 
          colc = "orange", colr = "violet", coll1 = "darkred", 
          coll2 = "darkgreen", size = 1.5, lwd = 0.7, colmesh1 = "red", 
          colmesh2 = "green", alpha = 0.3, tri = TRUE, outlines = TRUE, 
          points = TRUE, lines = FALSE, centroid = FALSE, cent.out = "CCA", 
          delta = 0.1, vecs = NULL) 
{
  if (dims == 3 & length(vecs) == 1) {
    tri <- FALSE
  }
  if (dims == 2 & any(c(outlines, lines, centroid))) {
    points <- TRUE
  }
  if (is.null(vecs) == TRUE) {
    vecs <- 1:dim(Shape$`3D_out`)[3]
  }
  if (lines == TRUE | centroid == TRUE) {
    centroids <- matrix(NA, ncol = 3, nrow = dim(Shape$`3D_out`)[3])
    if (cent.out == "CCA") {
      for (i in 1:dim(Shape$`3D_out`)[3]) {
        cent_x_y <- morphomapCentroid(Shape$`3D_out`[, 
                                                     , i], Shape$`3D_inn`[, , i], delta = delta)
        cent <- c(cent_x_y, mean(Shape$`3D_out`[, 
                                                3, i]))
        centroids[i, ] <- cent
      }
    }
    if (cent.out == "E") {
      for (i in 1:dim(Shape$`3D_out`)[3]) {
        centroids[i, ] <- colMeans(Shape$`3D_out`[, 
                                                  , i])
      }
    }
    if (cent.out == "I") {
      for (i in 1:dim(Shape$`3D_out`)[3]) {
        centroids[i, ] <- colMeans(Shape$`3D_out`[, 
                                                  , i])
      }
    }
  }
  if (dims == 3) {
    if (length(vecs) > 1) {
      all_coo_ext <- morphomapArray2matrix(Shape$`3D_out`[, 
                                                          , vecs])
      all_coo_int <- morphomapArray2matrix(Shape$`3D_inn`[, 
                                                          , vecs])
    }
    else {
      all_coo_ext <- Shape$`3D_out`[, , vecs]
      all_coo_int <- Shape$`3D_inn`[, , vecs]
    }
  }
  if (dims == 3) {
    open3d()
    if (centroid == TRUE) {
      spheres3d(centroids[vecs, ], col = colc, radius = size)
    }
    if (points == TRUE) {
      spheres3d(all_coo_ext, col = col1, radius = size)
      spheres3d(all_coo_int, col = col2, radius = size)
    }
    if (lines == TRUE) {
      for (i in 1:length(vecs)) {
        for (j in 1:dim(Shape$`3D_out`)[1]) {
          lines3d(rbind(Shape$`3D_out`[j, , i], 
                        centroids[vecs[i], ]), col = colr, lwd = lwd)
        }
      }
    }
    if (tri == TRUE) {
      mesh_ext <- morphomapTriangulate(all_coo_ext, length(vecs))
      mesh_int <- morphomapTriangulate(all_coo_int, length(vecs))
      triangles3d(t(mesh_ext$vb[, mesh_ext$it]), col = colmesh1, 
                  alpha = alpha, lit = T, specular = "black")
      triangles3d(t(mesh_int$vb[, mesh_ext$it]), col = colmesh2, 
                  alpha = alpha, lit = T, specular = "black")
    }
    if (outlines == TRUE) {
      for (i in 1:length(vecs)) {
        set1 <- Shape$`3D_out`[, , vecs[i]]
        set1 <- rbind(set1, set1[1, ])
        lines3d(set1[, 1], set1[, 2], set1[, 3], lwd = lwd, 
                col = 2)
        set2 <- Shape$`3D_inn`[, , vecs[i]]
        set2 <- rbind(set2, set2[1, ])
        lines3d(set2[, 1], set2[, 2], set2[, 3], lwd = lwd, 
                col = 3)
      }
    }
  }
  if (dims == 2) {
    all_coo_ext <- Shape$`2D_out`[, , vecs[1]]
    all_coo_int <- Shape$`2D_inn`[, , vecs[1]]
  }
  if (dims == 2) {
    if (points == TRUE) {
      plot(all_coo_ext, col = col1, cex = size, asp = 1, 
           xlab = "x", ylab = "y", main = paste("cross section number", 
                                                vecs[1]), pch = 19)
      points(all_coo_int, col = col2, cex = size, pch = 19)
    }
    if (outlines == TRUE) {
      set1 <- rbind(all_coo_ext, all_coo_ext[1, ])
      set2 <- rbind(all_coo_int, all_coo_int[1, ])
      points(set1, type = "l", col = coll1, lwd = lwd, 
             pch = 19)
      points(set2, type = "l", col = coll2, lwd = lwd, 
             pch = 19)
    }
    if (lines == TRUE) {
      for (i in 1:dim(all_coo_ext)[1]) {
        points(rbind(all_coo_ext[i, ], centroids[vecs[1], 
                                                 c(1, 2)]), type = "l", col = colr, lty = 3, 
               lwd = lwd, pch = 19)
      }
    }
    if (centroid == TRUE) {
      points(centroids[vecs[1], 1], centroids[vecs[1], 
                                              2], col = colc, cex = size, pch = 19)
    }
  }
}