#' morphomapThickness
#' 
#' Tool for the extraction of equiangular landmarks on the entire diaphysis 
#' @param morphomap.shape list: morphomap.shape object 
#' @return sect_thickness cortical thickness at each pair of landmarks on the external and internal outlines
#' @return ALPM_thickness cortical thickness at ALPM quadrants
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' \donttest{
#' #morphomap on a human femur bone
#' data(HomFem38023)
#' meshes<-morphomapSegm(HomFem38023, param1=4)
#' perMesh<-meshes$external
#' endMesh<-meshes$internal
#' mech_length<-380.23
#' rawSections<-morphomapCore(out.sur=perMesh,
#'                            inn.sur=endMesh,num.sect=61,mech.len = mech_length, 
#'                            start = 0.2,end=0.8)
#' shapeSections<-morphomapShape(rawSections,21,sects_vector=NULL,cent.out="CCA",delta=0.1)
#' femthick<-morphomapThickness(shapeSections)
#' plot(femthick$ALPM_thickness[1,,],type="l",
#'      main="LAMP thickness",xlab="section",ylab="thickness")
#' points(femthick$ALPM_thickness[2,,],type="l",col=2)
#' points(femthick$ALPM_thickness[3,,],type="l",col=3)
#' points(femthick$ALPM_thickness[4,,],type="l",col=4)
#' }
#' @export

morphomapThickness<-function(morphomap.shape) 
{
  num.land <- dim(morphomap.shape$`3D_out`)[1]
  num.sect <- dim(morphomap.shape$`3D_out`)[3]
  dist_sections_eq <- array(NA, dim = c(num.land, 1, num.sect))
  dist_sections_cp <- array(NA, dim = c(num.land, 1, num.sect))
  dist_ALPM <- array(NA, dim = c(4, 1, num.sect))
  for (m in 1:num.sect) {
    dist_section <- sqrt(rowSums((morphomap.shape$`2D_out`[, 
                                                           , m] - morphomap.shape$`2D_inn`[, , m])^2))
    closedPs <- vcgKDtree(morphomap.shape$`2D_inn`[, , m], 
                          morphomap.shape$`2D_out`[, , m], k = 1)$distance
    dist_sections_cp[, , m] <- closedPs
    dist_sections_eq[, , m] <- dist_section
  }
  for (m in 1:num.sect) {
    dist_quadr <- sqrt(rowSums((morphomap.shape$ALPM_out[, 
                                                         , m] - morphomap.shape$ALPM_inn[, , m])^2))
    dist_ALPM[, , m] <- dist_quadr
  }
  out <- list(sect_thickness_eq = dist_sections_eq, sect_thickness_cp = dist_sections_cp, 
              ALPM_thickness = dist_ALPM,morphomap.shape=morphomap.shape)
  return(out)
}
