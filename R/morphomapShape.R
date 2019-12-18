#' morphomapShape
#' 
#' Tool for the extraction of equiangular landmarks on the entire diaphysis 
#' @param morphomap.core list: morphomap.core object 
#' @param num.land numeric: number of landmarks defining each section
#' @param sects_vector numeric: number of sections
#' @param cent.out how to define the center of each section. The method allowed are "CCA" (center of cortical area), "E" (barycenter of the external outline) and "I" (barycenter of the internal outline)
#' @param delta pixel size used to calculate the CCA
#' @param side character: specify if the long bone is "left" or "right" side 
#' @return 3D_out num.pointsx3xnum.sect array in which the external outlines are stored
#' @return 3D_inn num.pointsx3xnum.sect array in which the internal outlines are stored
#' @return 2D_out num.pointsx2xnum.sect array in which the external outlines are stored
#' @return 2D_inn num.pointsx2xnum.sect array in which the interal outlines are stored
#' @return ALPM_inn array with the coordinates of ALPM coordinates on the external outline
#' @return ALPM_out array with the coordinates of ALPM coordinates on the internal outline
#' @return mech_length mechanical length of the long bone
#' @return start percentage of the mechanical length from which the first section is defined
#' @return end percentage of the mechanical length from which the last section is defined
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' \donttest{
#' data(HomFem38023)
#' meshes<-morphomapSegm(HomFem38023)
#' perMesh<-meshes$external
#' endMesh<-meshes$internal
#' mech_length<-380.23
#' rawSections<-morphomapCore(out.sur=perMesh,
#' inn.sur=endMesh,num.sect=61,
#' mech.len = mech_length, start = 0.2,end=0.8,num.points = 500)
#' # Shape coordinates defining as center the barycenter of the cortical area
#' shapeSections_CCA<-morphomapShape(rawSections,21,sects_vector=NULL, cent.out="CCA",delta=0.1)
#' sect1_ext<-shapeSections_CCA$`2D_out`[,,1]
#' sect1_int<-shapeSections_CCA$`2D_inn`[,,1]
#' centroid_CCA<-morphomapCentroid(rawSections$`2D_out`[,,1],rawSections$`2D_inn`[,,1], delta=0.1)
#' plot(sect1_ext,type="b",asp=1,xlab="x",ylab="y",main="Section 1 - CCA")
#' points(sect1_int,type="b",asp=1)
#' points(centroid_CCA[1],centroid_CCA[2],pch=19)
#' start<-c(centroid_CCA[1],centroid_CCA[2])
#' for(i in 1:21){
#'   to_e<-sect1_ext[i,]
#'   points(rbind(start,to_e),type="l")
#' }
#' 
#' # Shape coordinates defining as center the barycenter of the external perimeter
#' shapeSections_E<-morphomapShape(rawSections,21,sects_vector=NULL, cent.out="E",
#' delta=0.1, side="left")
#' sect1_ext<-shapeSections_E$`2D_out`[,,1]
#' sect1_int<-shapeSections_E$`2D_inn`[,,1]
#' centroid_E<-colMeans(rawSections$`2D_out`[,,1])
#' plot(sect1_ext,type="b",asp=1,xlab="x",ylab="y",main="Section 1 - E")
#' points(sect1_int,type="b",asp=1)
#' points(centroid_E[1],centroid_E[2],pch=19)
#' start<-c(centroid_E[1],centroid_E[2])
#' for(i in 1:21){
#'   to_e<-sect1_ext[i,]
#'   points(rbind(start,to_e),type="l")
#' }
#' 
#' # Shape coordinates defining as center the barycenter of the internal perimetershape
#' shapeSections_I<-morphomapShape(rawSections,21,sects_vector=NULL, cent.out="I",
#' delta=0.1, side="left")
#' sect1_ext<-shapeSections_I$`2D_out`[,,1]
#' sect1_int<-shapeSections_I$`2D_inn`[,,1]
#' centroid_I<-colMeans(rawSections$`2D_inn`[,,1])
#' plot(sect1_ext,type="b",asp=1,xlab="x",ylab="y",main="Section 1 - I")
#' points(sect1_int,type="b",asp=1)
#' points(centroid_I[1],centroid_I[2],pch=19)
#' start<-c(centroid_I[1],centroid_I[2])
#' for(i in 1:21){
#' to_e<-sect1_ext[i,]
#' points(rbind(start,to_e),type="l")
#' }
#' }
#' @export

morphomapShape<-function (morphomap.core, num.land, sects_vector, cent.out = "CCA", 
          delta = 0.1,side="left") 
{
mech.len <- morphomap.core$mech_length
start <- morphomap.core$start
end <- morphomap.core$end
if (is.null(sects_vector) == TRUE) {
  sects_vector <- 1:dim(morphomap.core$`3D_out`)[3]
}
num.sect <- length(sects_vector)
sect_poi <- seq(mech.len * start, mech.len * end, length = num.sect)
out_coo_3D <- array(NA, dim = c(num.land, 3, num.sect))
inn_coo_3D <- array(NA, dim = c(num.land, 3, num.sect))
out_coo_2D <- array(NA, dim = c(num.land, 2, num.sect))
inn_coo_2D <- array(NA, dim = c(num.land, 2, num.sect))
ALPM_i_tot <- array(NA, dim = c(4, 2, num.sect))
ALPM_o_tot <- array(NA, dim = c(4, 2, num.sect))

for (m in sects_vector) {
  out_outline <- morphomap.core$`2D_out`[, , m]
  inn_outline <- morphomap.core$`2D_inn`[, , m]
  if (cent.out == "CCA") {
    centroid <- morphomapCentroid(out_outline, inn_outline, 
                                  delta)
    ids_o <- morphomapRegradius(out_outline, centroid, 
                                num.land)
    ALPM_o <- morphomapRegradius(out_outline, centroid, 4)
    ids_i <- morphomapRegradius(inn_outline, centroid, 
                                num.land)
    ALPM_i <- morphomapRegradius(inn_outline, centroid,4)

    if(side=="right"){
      fho<-ids_o[which.min(abs(ids_o-ALPM_o[1])):(which.min(abs(ids_o-ALPM_o[3]))-1)]
      sho<-ids_o[which.min(abs(ids_o-ALPM_o[3])):length(ids_o)]
      out_coo_2D[, , m] <- out_outline[c(sho,fho)[c(1,num.land:2)], ]
      out_coo_3D[, , m] <- cbind(out_coo_2D[, , m], sect_poi[m])
      ALPM_o_tot[, , m] <- out_outline[ALPM_o[c(3,2,1,4)], ]
    }else{
      ALPM_o_tot[, , m] <- out_outline[ALPM_o, ]
      out_coo_2D[, , m] <- out_outline[ids_o, ]
      out_coo_3D[, , m] <- cbind(out_coo_2D[, , m], sect_poi[m])
    }  
    
    if(side=="right"){
      fhi<-ids_i[which.min(abs(ids_i-ALPM_i[1])):(which.min(abs(ids_i-ALPM_i[3]))-1)]
      shi<-ids_i[which.min(abs(ids_i-ALPM_i[3])):length(ids_i)]
      inn_coo_2D[, , m] <- inn_outline[c(shi,fhi)[c(1,num.land:2)], ]
      inn_coo_3D[, , m] <- cbind(inn_coo_2D[, , m], sect_poi[m])
      ALPM_i_tot[, , m] <- inn_outline[ALPM_o[c(3,2,1,4)], ]
    }else{
      ALPM_i_tot[, , m] <- inn_outline[ALPM_i, ]
      inn_coo_2D[, , m] <- inn_outline[ids_i, ]
      inn_coo_3D[, , m] <- cbind(inn_coo_2D[, , m], sect_poi[m])
    }    
    
  }
  if (cent.out == "E")  {
    ids_o <- morphomapRegradius(out_outline, colMeans(out_outline), 
                                num.land)
    ALPM_o <- morphomapRegradius(out_outline, colMeans(out_outline), 
                                 4)
    
    
    if(side=="right"){
      fho<-ids_o[which.min(abs(ids_o-ALPM_o[1])):(which.min(abs(ids_o-ALPM_o[3]))-1)]
      sho<-ids_o[which.min(abs(ids_o-ALPM_o[3])):length(ids_o)]
      out_coo_2D[, , m] <- out_outline[c(sho,fho)[c(1,num.land:2)], ]
      out_coo_3D[, , m] <- cbind(out_coo_2D[, , m], sect_poi[m])
      ALPM_o_tot[, , m] <- out_outline[ALPM_o[c(3,2,1,4)], ]
    }else{
      ALPM_o_tot[, , m] <- out_outline[ALPM_o, ]
      out_coo_2D[, , m] <- out_outline[ids_o, ]
      out_coo_3D[, , m] <- cbind(out_coo_2D[, , m], sect_poi[m])
    }
    
    ids_i <- morphomapRegradius(inn_outline, colMeans(out_outline), 
                                num.land)
    ALPM_i <- morphomapRegradius(inn_outline, colMeans(out_outline), 
                                 4)
    if(side=="right"){
      fhi<-ids_i[which.min(abs(ids_i-ALPM_i[1])):(which.min(abs(ids_i-ALPM_i[3]))-1)]
      shi<-ids_i[which.min(abs(ids_i-ALPM_i[3])):length(ids_i)]
      inn_coo_2D[, , m] <- inn_outline[c(shi,fhi)[c(1,num.land:2)], ]
      inn_coo_3D[, , m] <- cbind(inn_coo_2D[, , m], sect_poi[m])
      ALPM_i_tot[, , m] <- inn_outline[ALPM_i[c(3,2,1,4)], ]
    }else{
      ALPM_i_tot[, , m] <- inn_outline[ALPM_i, ]
      inn_coo_2D[, , m] <- inn_outline[ids_i, ]
      inn_coo_3D[, , m] <- cbind(inn_coo_2D[, , m], sect_poi[m])
    }
  }
  if (cent.out == "I")  {
    ids_o <- morphomapRegradius(out_outline, colMeans(inn_outline), 
                                num.land)
    ALPM_o <- morphomapRegradius(out_outline, colMeans(inn_outline), 
                                 4)
    if(side=="right"){
      fho<-ids_o[which.min(abs(ids_o-ALPM_o[1])):(which.min(abs(ids_o-ALPM_o[3]))-1)]
      sho<-ids_o[which.min(abs(ids_o-ALPM_o[3])):length(ids_o)]
      out_coo_2D[, , m] <- out_outline[c(sho,fho)[c(1,num.land:2)], ]
      out_coo_3D[, , m] <- cbind(out_coo_2D[, , m], sect_poi[m])
      ALPM_o_tot[, , m] <- out_outline[ALPM_o[c(3,2,1,4)], ]
    }else{
      ALPM_o_tot[, , m] <- out_outline[ALPM_o, ]
      out_coo_2D[, , m] <- out_outline[ids_o, ]
      out_coo_3D[, , m] <- cbind(out_coo_2D[, , m], sect_poi[m])
    }
    
    
    ids_i <- morphomapRegradius(inn_outline, colMeans(inn_outline), 
                                num.land)
    ALPM_i <- morphomapRegradius(inn_outline, colMeans(inn_outline), 
                                 4)
    
    if(side=="right"){
      fhi<-ids_i[which.min(abs(ids_i-ALPM_i[1])):(which.min(abs(ids_i-ALPM_i[3]))-1)]
      shi<-ids_i[which.min(abs(ids_i-ALPM_i[3])):length(ids_i)]
      inn_coo_2D[, , m] <- inn_outline[c(shi,fhi)[c(1,num.land:2)], ]
      inn_coo_3D[, , m] <- cbind(inn_coo_2D[, , m], sect_poi[m])
      ALPM_i_tot[, , m] <- inn_outline[ALPM_i[c(3,2,1,4)], ]
    }else{
      ALPM_i_tot[, , m] <- inn_outline[ALPM_i, ]
      inn_coo_2D[, , m] <- inn_outline[ids_i, ]
      inn_coo_3D[, , m] <- cbind(inn_coo_2D[, , m], sect_poi[m])
    }
    
    
  }
}

out <- list(`3D_out` = round(out_coo_3D, 2), 
            `3D_inn` = round(inn_coo_3D, 2), 
            `2D_out` = round(out_coo_2D, 2), 
            `2D_inn` = round(inn_coo_2D, 2), 
            ALPM_inn = ALPM_i_tot, 
            ALPM_out = ALPM_o_tot, 
            start = start, end = end, 
            mech.len = mech.len)
return(out)
}