#' morphomapTriangulate
#' 
#' Build a mesh starting from the coordinates of the diaphysis
#' @param set matrix: coordinates of the cross sections to be triangulated
#' @param n numeric: number of cross sections
#' @param close logical: if TRUE the two surfaces are closed
#' @return mesh a mesh of the triangulated semilandark configuration
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @export


morphomapTriangulate<-function(set,n,close=FALSE){
  start<-seq(1,dim(set)[1],dim(set)[1]/n)
  end<-seq(dim(set)[1]/n,dim(set)[1],dim(set)[1]/n)
  it<-NULL
  vb<-NULL
  nland<-(dim(set)[1]/n)
  maxit<-0
  for(i in 1:(n-1)){
    pro_z<-morphomapTri2sects(set[start[i]:end[i],],set[start[i+1]:end[i+1],]) 
    it<-rbind(it,pro_z$tri+maxit)
    vb<-rbind(vb,pro_z$matrix)
    maxit<-maxit+nland*2
  }
  mesh_d<-list("vb"=t(cbind(vb,1)),"it"=t(it))
  class(mesh_d)<-"mesh3d"
  
  if(close==TRUE){
  center1_e<-colMeans(set[1:nland,])
  cup1_d_e<-morphomapTri2sects(repmat(center1_e,nland,3),set[1:nland,])
  cup1_m_e<-list("vb"=t(cbind(cup1_d_e$matrix,1)),"it"=t(cup1_d_e$tri))
  class(cup1_m_e)<-"mesh3d"
  
  center2_e<-colMeans(set[(dim(set)[1]-nland+1):dim(set)[1],])
  cup2_d_e<-morphomapTri2sects(set[(dim(set)[1]-nland+1):dim(set)[1],],repmat(center2_e,nland,3))
  cup2_m_e<-list("vb"=t(cbind(cup2_d_e$matrix,1)),"it"=t(cup2_d_e$tri))
  class(cup2_m_e)<-"mesh3d"
  
  mesh<-mergeMeshes(cup1_m_e,mesh_d,cup2_m_e)
  }else {mesh<-mesh_d}
  return(mesh)
}