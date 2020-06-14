#' morphomapTri2sects
#' 
#' Triangulate the external and internal outlines of a 3D cross section
#' @param cp matrix: coordinates of the external outline of the section
#' @param mp matrix: coordinates of the internal outline of the section
#' @return matrix coordinates of the triangulated mesh
#' @return tri triangulations of the triangulated mesh
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @export


morphomapTri2sects<-function(cp,mp){
  mat<-matrix(NA,ncol=3,nrow=dim(cp)[1]*2)
  totsec<-rbind(cp,mp)
  extpoi<-1:dim(cp)[1]
  intpoi<-(dim(mp)[1]+1):dim(totsec)[1]
  counter<-0
  for(i in 1:length(extpoi)){
    if(i!=length(extpoi)){
      mati<-c(extpoi[i],intpoi[i],extpoi[i+1])
      mati2<-c(intpoi[i],intpoi[i+1],extpoi[i+1])
      
    }
    if(i==length(extpoi)){
      mati<-c(extpoi[i],intpoi[i],extpoi[1])
      mati2<-c(intpoi[i],intpoi[1],extpoi[1])
      
    }
    
    counter<-counter+1
    mat[counter,]<-mati
    counter<-counter+1
    mat[counter,]<-mati2
  }
  
  out<-list("matrix"=totsec,"tri"=mat)
  return(out)
}