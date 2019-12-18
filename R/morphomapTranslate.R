#' morphomapTranslate
#' 
#' Translate a section to a new center defined by the user
#' @param corA matrix: coordinates of the external outline
#' @param medA matrix: coordinates of the internal outline
#' @param Cx numeric: new x center coordinate
#' @param Cy numeric: new y center coordinate
#' @return cortical new centered coordinates of the external outline
#' @return medullar new centered coordinates of the internal outline
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' extsec<-morphomapCircle(10,1000)
#' intsec<-morphomapCircle(8,1000)
#' plot(extsec,asp=1,type="l",xlim=c(-11,11),ylim=c(-11,11))
#' points(intsec,type="l")
#' traSect<-morphomapTranslate(extsec,intsec,1,1)
#' points(traSect$cortical,type="l",col="red")
#' points(traSect$medullar,type="l",col="red")
#' @export

morphomapTranslate<-function(corA,medA,Cx,Cy){
  Bar_s<-c(Cx,Cy)
  corAl<-corA
  medAl<-medA
  corAl[,1]<- corAl[,1]-Bar_s[1]
  corAl[,2]<- corAl[,2]-Bar_s[2]
  medAl[,1]<- medAl[,1]-Bar_s[1]
  medAl[,2]<- medAl[,2]-Bar_s[2]
  out<-list("cortical"=corAl,"medullar"=medAl)
  return(out)
}