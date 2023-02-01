#' morphomapImport
#'
#' Import a morphomapShape object exported with morphomapExport 
#' @param file character: name of input file
#' @return 3D_out num.pointsx3xnum.sect array in which the external outlines are stored
#' @return 3D_inn num.pointsx3xnum.sect array in which the internal outlines are stored
#' @return 2D_out num.pointsx2xnum.sect array in which the external outlines are stored
#' @return 2D_inn num.pointsx2xnum.sect array in which the interal outlines are stored
#' @return ALPM_inn array with the coordinates of ALPM coordinates on the external outline
#' @return ALPM_out array with the coordinates of ALPM coordinates on the internal outline
#' @return mech_length mechanical length of the long bone
#' @return start percentage of the mechanical length from which the first section is defined
#' @return end percentage of the mechanical length from which the last section is defined
#' @author Antonio Profico
#' @export

morphomapImport<-function(file){
  A <- readLines(file)
  e1<-which(A == "@1");e2<-which(A == "@2");e3<-which(A == "@3")
  e4<-which(A == "@4");e5<-which(A == "@5");e6<-which(A == "@6")
  e7<-which(A == "@7");e8<-which(A == "@8");e9<-which(A == "@9")
  e10<-which(A == "@10"); e11<-which(A == "@11")
  
  num.lands<-as.numeric(read.table(file, skip = e11, nrows =(e11-e10-2)))
  num.sects<-as.numeric(read.table(file, skip = e10, nrows =(e10-e9-2)))
  se3D_out<-morphomapMatrix2array(as.matrix(read.table(file, skip = e1, nrows =(e2-e1-2))),nsects =num.sects)
  se3D_inn<-morphomapMatrix2array(as.matrix(read.table(file, skip = e2, nrows =(e3-e2-2))),nsects =num.sects)
  se2D_out<-morphomapMatrix2array(as.matrix(read.table(file, skip = e3, nrows =(e4-e3-2))),nsects =num.sects)
  se2D_inn<-morphomapMatrix2array(as.matrix(read.table(file, skip = e4, nrows =(e5-e4-2))),nsects =num.sects)
  ALPM_inn<-morphomapMatrix2array(as.matrix(read.table(file, skip = e5, nrows =(e6-e5-2))),nsects =num.sects)
  ALPM_out<-morphomapMatrix2array(as.matrix(read.table(file, skip = e6, nrows =(e7-e6-2))),nsects =num.sects)
  start<-as.numeric(read.table(file, skip = e7, nrows =(e8-e7-2)))
  end<-as.numeric(read.table(file, skip = e8, nrows =(e9-e8-2)))
  mech.len<-as.numeric(read.table(file, skip = e9, nrows =(e11-e10-2)))
  
  out<-list("3D_out"= se3D_out,"3D_inn"=se3D_inn,"2D_out"=se2D_out,
             "2D_inn"=se2D_inn,"ALPM_inn"=ALPM_inn,"ALPM_out"=ALPM_out,
            "start"=start,"end"=end,"mech.len"=mech.len)
  
  return(out)
}
  

