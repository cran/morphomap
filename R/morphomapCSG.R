#' morphomapCSG
#' 
#' Tool for Cross-sectional geometry
#' @param cp matrix: coordinates of the external outline
#' @param mp matrix: coordinates of the internal outline
#' @param translate logical: if TRUE the section will be centered
#' @param center how to define the center of each section. The method allowed are "CCA" (center of cortical area), "E" (barycenter of the external outline) and "I" (barycenter of the internal outline)
#' @param delta numeric: picture elements of adjustable side length
#' @param Cx numeric: new x center coordinate
#' @param Cy numeric: new y center coordinate
#' @param I_xy logical: if TRUE the product of inertia around the x and y axis is calculated
#' @param I_minmax logical: if TRUE the Imin and Imax will be calculated
#' @param Zxy logical: if TRUE the polar moment of inertia will be calculated
#' @return Cx x coordinate of the centered section
#' @return Cy y coordinate of the centered section
#' @return T_area total area
#' @return M_area medullar area
#' @return CA cortical area
#' @return Ext_perim external perimeter
#' @return Med_perim medullar perimiter
#' @return Mean_thick mean thickness of the section
#' @return Sd_thick thickness standard deviation
#' @return Min_thick minimum thickness
#' @return Max_thick maximum thickness
#' @return Ix numeric: moment of inertia around the x axis
#' @return Iy numeric: moment of inertia around the y axis
#' @return Zx numeric: moment of inertia around the x axis
#' @return Zy numeric: moment of inertia around the y axis
#' @return Zpol numeric: polar moment of inertia 
#' @return dx new centered coordinates of the internal outline
#' @return dy new centered coordinates of the internal outline
#' @return Imin numeric: minimum moment of inertia
#' @return Imax numeric: maximum moment of inertia
#' @return J numeric: polar moment of inertia
#' @return Zmax numeric: the maximum polar section
#' @return Zmin numeric: the minimum polar section
#' @return theta numeric: theta angle
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' \donttest{
#' #calculation of csg parameter on a human femur cross section
#' data(HomFem38023)
#' meshes<-morphomapSegm(HomFem38023)
#' perMesh<-meshes$external
#' endMesh<-meshes$internal
#' mech_length<-380.23
#' rawSections<-morphomapCore(out.sur=perMesh,
#'                            inn.sur=endMesh,num.sect=61,mech.len = mech_length, 
#'                            start = 0.2,end=0.8)
#' shapeSections<-morphomapShape(rawSections,250,sects_vector=NULL,cent.out="CCA",
#' delta=0.1, side="left")
#' csgSect31<-morphomapCSG(cp = shapeSections$`2D_out`[,,31],
#'                         mp=shapeSections$`2D_inn`[,,31],
#'                         translate = FALSE,center="CCA")
#' 
#' #Cross sectional geometry along the entire femur bone
#' results<-matrix(NA,ncol=24,nrow=61)
#' rownames(results)<-paste("section",c(1:61))
#' colnames(results)<-c("Cx","Cy","T_area","M_area","CA",
#'                      "Ext_perim","Med_perim","Mean_thick","Sd_thick" ,
#'                      "Min_thick","Max_thick","Ix","Iy","Zx" ,"Zy","Zpol" ,
#'                      "dx","dy","Imin","Imax","J","Zmax","Zmin","theta")
#' 
#' for(i in 1:61){
#'   results[i,]<-unlist(morphomapCSG(cp = shapeSections$`2D_out`[,,i],
#'                                      mp=shapeSections$`2D_inn`[,,i],
#'                                    translate = FALSE,center="CCA",delta = 0.5))
#' }
#' 
#' plot(c(1:61),results[,24],type="b",main="Theta",cex=1,
#'      xlab="section",ylab="radians")
#' 
#' 
#' #calculation of csg parameter on a chimpanzee femur cross section
#' data(PanFem27713)
#' meshes<-morphomapSegm(PanFem27713)
#' perMesh<-meshes$external
#' endMesh<-meshes$internal
#' mech_length<-277.13
#' rawSections<-morphomapCore(out.sur=perMesh,
#'                            inn.sur=endMesh,num.sect=61,mech.len = mech_length, 
#'                            start = 0.2,end=0.8)
#' shapeSections<-morphomapShape(rawSections,250,sects_vector=NULL,cent.out="CCA",
#' delta=0.1, side="left")
#' csgSect31<-morphomapCSG(cp = shapeSections$`2D_out`[,,31],
#'                         mp=shapeSections$`2D_inn`[,,31],
#'                         translate = FALSE,center="CCA")
#' 
#' #Cross sectional geometry along the entire femur bone
#' results<-matrix(NA,ncol=24,nrow=61)
#' rownames(results)<-paste("section",c(1:61))
#' colnames(results)<-c("Cx","Cy","T_area","M_area","CA",
#'                      "Ext_perim","Med_perim","Mean_thick","Sd_thick" ,
#'                      "Min_thick","Max_thick","Ix","Iy","Zx" ,"Zy","Zpol" ,
#'                      "dx","dy","Imin","Imax","J","Zmax","Zmin","theta")
#' 
#' for(i in 1:61){
#'   results[i,]<-unlist(morphomapCSG(cp = shapeSections$`2D_out`[,,i],
#'                                    mp=shapeSections$`2D_inn`[,,i],
#'                                    translate = FALSE,center="CCA",delta = 0.5))
#'   }
#' 
#' plot(c(1:61),results[,24],type="b",main="Theta",cex=1,
#'      xlab="section",ylab="radians")
#' }
#'@export

morphomapCSG<-function (cp, mp, translate = FALSE, center = c("I", "E", "CCA"), 
                        delta = 0.1, Cx = NULL, Cy = NULL, I_xy = TRUE, I_minmax = TRUE, 
                        Zxy = TRUE) 
{
  out <- list(Cx = NULL, Cy = NULL, T_area = NULL, M_area = NULL, 
              CA = NULL, Ext_perim = NULL, Med_perim = NULL, Mean_thick = NULL, 
              Sd_thick = NULL, Min_thick = NULL, Max_thick = NULL, 
              Ix = NULL, Iy = NULL, Zx = NULL, Zy = NULL, Zpol = NULL, 
              dx = NULL, dy = NULL, Imin = NULL, Imax = NULL, J = NULL, 
              Zmax = NULL, Zmin = NULL, theta = NULL)
  cp <- rbind(cp, cp[1, ])
  mp <- rbind(mp, mp[1, ])
  if (translate == TRUE) {
    tra <- morphomapTranslate(cp, mp, Cx, Cy)
    cp <- tra$cortical
    mp <- tra$medullar
  }
  if (translate == FALSE) {
    if (center == "I") {
      origin <- colMeans(mp)
      Cx <- origin[1]
      Cy <- origin[2]
      tra <- morphomapTranslate(cp, mp, Cx, Cy)
      cp <- tra$cortical
      mp <- tra$medullar
    }
    if (center == "E") {
      origin <- colMeans(cp)
      Cx <- origin[1]
      Cy <- origin[2]
      tra <- morphomapTranslate(cp, mp, Cx, Cy)
      cp <- tra$cortical
      mp <- tra$medullar
    }
    if (center == "CCA") {
      origin <- morphomapCentroid(cp, mp, delta = delta)
      Cx <- origin[1]
      Cy <- origin[2]
      tra <- morphomapTranslate(cp, mp, Cx, Cy)
      cp <- tra$cortical
      mp <- tra$medullar
    }
  }
  out$Cx <- Cx
  out$Cy <- Cy
  area_c <- morphomapArea(cp)
  area_m <- morphomapArea(mp)
  CA <- area_c - area_m
  out$T_area <- area_c
  out$M_area <- area_m
  out$CA <- CA
  perim_c <- NULL
  perim_m <- NULL
  for (i in 1:dim(cp)[1]) {
    if (i == dim(cp)[1]) {
      perim_c[i] <- sqrt(sum((cp[i, ] - cp[1, ])^2))
    }
    else {
      perim_c[i] <- sqrt(sum((cp[i, ] - cp[i + 1, ])^2))
    }
  }
  for (i in 1:dim(mp)[1]) {
    if (i == dim(mp)[1]) {
      perim_m[i] <- sqrt(sum((mp[i, ] - mp[1, ])^2))
    }
    else {
      perim_m[i] <- sqrt(sum((mp[i, ] - mp[i + 1, ])^2))
    }
  }
  cor_perim <- sum(perim_c)
  med_perim <- sum(perim_m)
  dist_section <- sqrt(rowSums((cp - mp)^2))
  mean_thick <- mean(dist_section)
  sd_thick <- sd(dist_section)
  min_thick <- min(dist_section)
  max_thick <- max(dist_section)
  out$Ext_perim <- cor_perim
  out$Med_perim <- med_perim
  out$Mean_thick <- mean_thick
  out$Sd_thick <- sd_thick
  out$Min_thick <- min_thick
  out$Max_thick <- max_thick
  if (I_xy == TRUE) {
    I_xy <- morphomapMoment(cp, mp, delta = delta)
    Ix <- I_xy[1]
    Iy <- I_xy[2]
    out$Ix <- Ix
    out$Iy <- Iy
    out$J <- Ix + Iy
  }
  if (Zxy == TRUE) {
    I_xy_t <- morphomapZmoment(cp, mp, Cx = 0, Cy = 0, delta = delta)
    Zx <- I_xy_t[1]
    Zy <- I_xy_t[2]
    out$Zx <- Zx
    out$Zy <- Zy
    out$dx <- I_xy_t[3]
    out$dy <- I_xy_t[4]
  }
  if (I_minmax == TRUE) {
    I_xy <- morphomapMoment(cp, mp, delta = delta)
    theta <- atan(2 * I_xy[3]/(I_xy[1] - I_xy[2]))/2
    rotcp <- Rotate(cp[, 1], y = cp[, 2], mx = 0, my = 0, 
                    theta = theta, asp = 1)
    rotmp <- Rotate(mp[, 1], y = mp[, 2], mx = 0, my = 0, 
                    theta = theta, asp = 1)
    moms <- morphomapMoment(cbind(rotcp$x, rotcp$y), cbind(rotmp$x, 
                                                           rotmp$y), delta = delta)
    I_1 <- ((moms[1] + moms[2])/2) + sqrt((((moms[1] - moms[2])/2)^2) + 
                                            moms[3]^2)
    I_2 <- ((moms[1] + moms[2])/2) - sqrt((((moms[1] - moms[2])/2)^2) + 
                                            moms[3]^2)
    out$theta <- theta
    out$Imin <- min(c(I_1,I_2))
    out$Imax <- max(c(I_1,I_2))
    if (Zxy == TRUE) {
      Z_xy_t <- morphomapZmoment(cbind(rotcp$x, rotcp$y), 
                                 cbind(rotmp$x, rotmp$y), Cx = 0, Cy = 0, delta = delta)
      ddx <- Z_xy_t[3]
      ddy <- Z_xy_t[4]
      Zmax <- max(Z_xy_t[c(1, 2)])
      Zmin <- min(Z_xy_t[c(1, 2)])
      Zpol <- Z_xy_t[5]
      out$Zpol <- Zpol
      out$Zmax <- Zmax
      out$Zmin <- Zmin
    }
  }
  return(out)
}