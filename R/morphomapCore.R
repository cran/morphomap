#' morphomapCore
#'
#' Tool to build 3D and 2D cross sections 
#' @param out.sur object of class mesh3d
#' @param inn.sur object of class mesh3d
#' @param num.sect number of sections
#' @param mech.len mechanical length of the long bone
#' @param num.points number of equiengular points to be defined on each section
#' @param start percentage of the mechanical length from which the first section is defined
#' @param end percentage of the mechanical length from which the last section is defined
#' @param clean_int_out logical if TRUE the inner section will be cleaned by using spherical flipping
#' @param param1 numeric parameter for spherical flipping operator (how much the section will be deformed)
#' @param radius.fact numeric parameter for spherical flipping operator (distance from the center of the outline at which the povs are defined)
#' @param npovs numeric: number of points of view defined around the section
#' @param print.progress logical: if TRUE a progress bar is printed to the screen
#' @return 3D_out num.pointsx3xnum.sect array of the external outlines 
#' @return 3D_inn num.pointsx3xnum.sect array of the internal outlines
#' @return 2D_out num.pointsx2xnum.sect array of the external outlines
#' @return 2D_inn num.pointsx2xnum.sect array of the internal outlines
#' @return mech_length mechanical length of the long bone
#' @return start percentage of the mechanical length from which the first section is defined
#' @return end percentage of the mechanical length from which the last section is defined
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' \donttest{
#' #raw section on a human femur bone
#' data(HomFem38023)
#' meshes<-morphomapSegm(HomFem38023, param1=4)
#' perMesh<-meshes$external
#' endMesh<-meshes$internal
#' mech_length<-380.23
#' rawSections<-morphomapCore(out.sur=perMesh,
#'                            inn.sur=endMesh,num.sect=61,mech.len = mech_length,
#'                            start = 0.2,end=0.8)
#' #2D plot of the first section
#' plot(rawSections$`2D_out`[,,1],col="grey",asp=1,xlab="x",ylab="y",type="l")
#' points(rawSections$`2D_inn`[,,1],col="red",type="l")
#' #3D plot of the first section
#' library(rgl)
#' open3d()
#' plot3d(rawSections$`3D_out`[,,1],aspect=FALSE,col="grey",type="l",lwd=5,xlab="x",ylab="y",zlab="z")
#' plot3d(rawSections$`3D_inn`[,,1],aspect=FALSE,col="red",type="l",lwd=5,add=TRUE)
#'
#' #raw section on a chimpanzee femur bone
#' data(PanFem27713)
#' meshes<-morphomapSegm(PanFem27713, param1=3)
#' perMesh<-meshes$external
#' endMesh<-meshes$internal
#' mech_length<-277.13
#' rawSections<-morphomapCore(out.sur=perMesh,
#'                            inn.sur=endMesh,num.sect=61,mech.len = mech_length,
#'                            start = 0.2,end=0.8)
#' #2D plot of the first section
#' plot(rawSections$`2D_out`[,,1],col="grey",asp=1,xlab="x",ylab="y",type="l")
#' points(rawSections$`2D_inn`[,,1],col="red",type="l")
#' #3D plot of the first section
#' library(rgl)
#' open3d()
#' plot3d(rawSections$`3D_out`[,,1],aspect=FALSE,col="grey",type="l",lwd=5,xlab="x",ylab="y",zlab="z")
#' plot3d(rawSections$`3D_inn`[,,1],aspect=FALSE,col="red",type="l",lwd=5,add=TRUE)
#' }
#' @export

morphomapCore<-function (out.sur = out.sur, inn.sur = inn.sur, num.sect = 61,
                         mech.len, clean_int_out = TRUE, param1 = 0.5, radius.fact = 2.5,
                         npovs = 100, num.points = 500, start = 0.2, end = 0.8, print.progress = TRUE)
{
  ext_raw_sects_2D <- array(NA, dim = c(num.points, 2, num.sect))
  inn_raw_sects_2D <- array(NA, dim = c(num.points, 2, num.sect))
  ext_raw_sects_3D <- array(NA, dim = c(num.points, 3, num.sect))
  inn_raw_sects_3D <- array(NA, dim = c(num.points, 3, num.sect))
  sect_poi <- seq(mech.len * start, mech.len * end, length = num.sect)

  if(print.progress==TRUE){
    pb <- txtProgressBar(min=0,max=num.sect-1,initial=0,style=3)
    step<-0
  }

  for (m in 1:num.sect) {
    p1 <- c(0, 0, sect_poi[m])
    p2 <- c(100, 0, sect_poi[m])
    p3 <- c(0, 100, sect_poi[m])
    inter_out <- NULL
    inter_inn <- NULL
    inter_out <- meshPlaneIntersect(out.sur, p1, p2, p3)[,
                                                         c(1, 2)]
    inter_inn <- meshPlaneIntersect(inn.sur, p1, p2, p3)[,
                                                         c(1, 2)]
    inters <- inter_inn
    ordered_out_temp <- morphomapSort(inter_out)
    ordered_out_temp <- rbind(ordered_out_temp, ordered_out_temp[1,
                                                                 ])
    ordered_out <- ordered_out_temp
    if (clean_int_out == TRUE) {
      inter_inn <- morphomapFlip(inter_inn, param1 = param1,
                                 radius.fact = radius.fact, npovs = npovs)
    }
    ordered_inn_temp <- morphomapSort(inter_inn)
    ordered_inn_temp <- rbind(ordered_inn_temp, ordered_inn_temp[1,
                                                                 ])
    ordered_inn <- ordered_inn_temp
    ev_out <- equidistantCurve(ordered_out, n = num.points,
                               iterations = 1, increment = 0)
    ev_inn <- equidistantCurve(ordered_inn, n = num.points,
                               iterations = 1, increment = 0)
    ext_raw_sects_2D[, , m] <- ev_out
    inn_raw_sects_2D[, , m] <- ev_inn
    ext_raw_sects_3D[, , m] <- cbind(ev_out, sect_poi[m])
    inn_raw_sects_3D[, , m] <- cbind(ev_inn, sect_poi[m])
    if(print.progress==TRUE){
      setTxtProgressBar(pb,step)
      step <- step+1
    }
  }
  if(print.progress==TRUE){
    close(pb)
  }

  out <- list(`3D_out` = ext_raw_sects_3D, `3D_inn` = inn_raw_sects_3D,
              `2D_out` = ext_raw_sects_2D, `2D_inn` = inn_raw_sects_2D,
              mech_length = mech.len, start = start, end = end)
  return(out)
}
