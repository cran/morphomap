#' morphomapRaster
#' 
#' Convert a section in a raster image. It is useful to save cross section at the real size
#' @param cp numeric: radius of the outline 
#' @param mp numeric: number of points along the outline
#' @param pixel numeric: desired ratio pixel/mm
#' @param filename character: path of the file to be saved
#' @param save logical: if TRUE the raster image will be saved
#' @return rimg raster image of the cross section
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' library(raster)
#' #rectangular section
#' extsec<-morphomapRectangle(10,6,100)
#' intsec<-morphomapRectangle(8,4,100)
#' rimg<-morphomapRaster(extsec,intsec,pixel=0.1,save=FALSE)
#' plot(rimg,col=gray(0:2/2))
#' #circular section
#' extsec<-morphomapCircle(10,100)
#' intsec<-morphomapCircle(8,100)
#' rimg<-morphomapRaster(extsec,intsec,pixel=0.1,save=FALSE)
#' plot(rimg,col=gray(0:2/2))
#' @export

morphomapRaster<-function(cp,mp,pixel=1,filename,save=FALSE){
XX <- extendrange(cp[, 1])
YY <- extendrange(cp[, 2])
maxx <- max(XX)
minx <- min(XX)
maxy <- max(YY)
miny <- min(YY)
X <- seq(minx + pixel/2, maxx - pixel/2, pixel)
Y <- seq(miny + pixel/2, maxy - pixel/2, pixel)
M <- matrix(0, length(X), length(Y))
grid_sect <- as.matrix(expand.grid(X, Y))
A <- point.in.polygon(grid_sect[, 1], grid_sect[, 2], cp[, 
                                                         1], cp[, 2], mode.checked = FALSE)
B <- point.in.polygon(grid_sect[, 1], grid_sect[, 2], mp[, 
                                                         1], mp[, 2], mode.checked = FALSE)
sel <- which(A == 1 & B == 0)
selt <- rep(0, length(X) * length(Y))
selt[sel] <- 1
img <- list()
img$x <- (X)
img$y <- (Y)
img$z <- (matrix(t(selt), length(X), length(Y), byrow = F))

check_dx<-dim(img$z)[1]
check_dy<-dim(img$z)[2]
if(check_dx != check_dy){
if(check_dx>check_dy){
 diff<-check_dx-check_dy 
for(i in 1:diff){
  img$z<-cbind(img$z,0)
}
 if(check_dx<check_dy){
   diff<-check_dy-check_dx
   for(i in 1:diff){
     img$z<-rbind(img$z,0)
   }
 } 
}
}

r <- raster(t(img$z), xmn = 0, xmx = dim(img$z)[2], 
            ymn = 0, ymx = dim(img$z)[1], crs = CRS("+proj=utm +zone=11 +datum=NAD83"))
rimg <- flip(r, 2)
if (save == TRUE) {
  writeRaster(rimg, filename, "GTiff",overwrite=TRUE)
}
return(rimg)
}