#' morphomapFlip
#' 
#' Spherical flipping operator for bi-dimensional configuration
#' @param mat numeric matrix: coordinates of the bi-dimensional configuration
#' @param param1 numeric: first parameter for spherical flipping 
#' @param param2 numeric: second parameter for spherical flipping 
#' @param radius.fact mechanical length of the long bone
#' @param npovs number of evenly spaced points to be defined on each section 
#' @return mat matrix after spherical flipping
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' #create a section
#' extsec<-morphomapCircle(10,1000)
#' intsec<-morphomapCircle(8,1000)
#' #simulate noise 
#' noiseX<-rnorm(1000,mean = 0,sd = 0.2)
#' noiseY<-rnorm(1000,mean = 0,sd = 0.2)
#' noise<-cbind(noiseX,noiseY)
#' noisect<-intsec+noise
#' #spherical flipping
#' flipsect<-morphomapFlip(noisect,param1 = 2,radius.fact = 2)
#' sortsect<-morphomapSort(flipsect)
#' #original section
#' plot(extsec,asp=1,type="l",xlim=c(-15,15),ylim=c(-15,15))
#' points(intsec,asp=1,type="l",xlim=c(-15,15),ylim=c(-15,15))
#' #noise
#' points(noisect,col=2)
#' #new section after spherical flipping
#' points(sortsect,type="l",col=3,asp=1,lwd=2)
#' @export

morphomapFlip<-function(
  mat,param1=0.8,param2=10,
  radius.fact=1.5,npovs=100){ 
  
  numPoi<-dim(mat)[1]
  P <- mat
  center<-colMeans(P)
  matcent<-repmat(matrix(center,ncol=2),dim(P)[1],1)
  h<-center[1]
  k<-center[2]
  max_radius<-max(sqrt(rowSums((P-matcent)^2)))
  radius<-max_radius*radius.fact
  theta_seq<-seq(0,2*pi,length.out = 100)
  mat<-matrix(NA,ncol=2,nrow=npovs)
  for(i in 1:npovs){
    x = h + radius*cos(theta_seq[i])
    y = k + radius*sin(theta_seq[i]) 
    mat[i,1]<-x
    mat[i,2]<-y
  }
  sel<-NULL
  for(i in 1:npovs){
    C <- matrix(mat[i,],ncol=2,nrow=1)
    P2 <- P - (repmat(C, numPoi, 1))
    normp <- rowSums(P2^2)
    normp <- sqrt(normp)
    param <- param1
    R2 <- matrix(repmat(max(normp) * (param2^param1), 2, 1))
    SF <- P2 + 2 * repmat(R2[1] - cbind(normp), 1, 2) * P2/(repmat(cbind(normp),1, 2))
    cloud <- rbind(C, SF, P)
    tri <- t(convhulln(cloud, options = ""))
    subset <- unique(as.vector(tri))
    subset <- subset[which(subset > 1 & subset <= (dim(P)[1] + 1))]
    subset <- subset - 1
    sel<-c(sel,subset)
  }
  sel<-unique(sel)
  mat<-P[sel,]
  return(mat)
}