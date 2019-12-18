#' morphomapSort
#' 
#' Sort a series of points stored as a 2D matrix 
#' @param mat numeric matrix: a kx2 matrix 
#' @return mat sorted kx2 matrix
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' rand<-sample(100)
#' extsec<-morphomapCircle(10,100)[rand,]
#' plot(extsec,type="l",asp=1)
#' sorted<-morphomapSort(extsec)
#' plot(sorted,type="l",asp=1)
#' @export

morphomapSort<-function(mat){
  center_x<-mean(mat[,1])
  center_y<-mean(mat[,2])
  
  diff_x<-mat[,1]-center_x
  diff_y<-mat[,2]-center_y
  
  angles <- atan2(diff_y, diff_x)
  order_1<-order(angles)
  
  diff_A<-mat[,1]-center_x
  diff_pos_A<-which(diff_A>=0)
  start_A<-diff_pos_A[which.min(diff_A[diff_pos_A])]
  
  start_A2<-which(order_1==start_A)
  new_order<-c(order_1[start_A2:length(order_1)],order_1[1:(start_A2-1)])
  
  mat<-mat[new_order,]
  return(mat)
}