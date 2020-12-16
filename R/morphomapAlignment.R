#' morphomapAlignment
#'
#' Align a femur bone following the protocol proposed by Ruff (2002)
#' @param mesh 3D mesh: femur long bone mesh
#' @param set matrix: 7 landmarks acquired on the mesh (see details)
#' @param side character: specify if the femur bone is "left" or "right" side
#' @param param1 numeric: parameter for spherical flipping (usually ranged between 3 and 4)
#' @param iter1 numeric: number of iterations first alignment
#' @param iter2 numeric: number of iterations second alignment
#' @param iter3 numeric: number of iterations third alignment
#' @param from1 numeric: inferior range of the allowed rotation in the first alignment  
#' @param to1 numeric: superior range of the allowed rotation in the first alignment 
#' @param from2 numeric: inferior range of the allowed rotation in the second alignment 
#' @param to2 numeric: superior range of the allowed rotation in the second alignment 
#' @param from3 numeric: inferior range of the allowed rotation in the third alignment 
#' @param to3 numeric: superior range of the allowed rotation in the third alignment 
#' @param tol numeric: maximum allowed error in the alignment expressed in mm 
#' @return sur: mesh of the aligned femur bone
#' @return coo: coordinates of the landmark used in the alignment (plus two added automatically)
#' @return mech_length: mechanical length of the aligned femur bone
#' @details The function 'morphomapAlignment' is designed to align a femur bone. I did not tested on other long bones. 
#' The function requires 7 anatomical landmarks samples as follow: 1-the point at the center of the diaphysis in posterior view after the less trochanter, 
#' 2- the most posterior point on the lateral epicondyle, 3-the most posterior point on the medial epicondyle, 4- the most inferior point on the intercondilar fossa,
#' 5- neck of the femur, 6- the most inferior point on the medial epicondyle and 7-the most inferior point on the lateral epicondyle. 
#' If the function in a short time does not complete the alignement, please stop the R session, check your landmark configuration or try to increase the value of the argument 'tol'.  
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @export
morphomapAlignment<-function(mesh,set,side=c("left","right"),
                                   param1=4,iter1=2000,iter2=2000,
                                   iter3=2000,from1=180,to1=360,
                                   from2=-5,to2=5,from3=-5,to3=5,
                                   tol=0.5){
morphomapAlOne<-function(mesh,set,iter,from,to,tol){
    
    pos_1<-aro.clo.points(t(mesh$vb)[,-4],set)$position
    set_1<-t(mesh$vb)[pos_1, -4] 
    sur_2 <- morphomapAlMesh(rbind(set_1[2, ], colMeans(set_1[c(1,2),]),set_1[1,]),mesh)
    set_2<-t(sur_2$mesh$vb)[pos_1, -4]
    sur_2<-sur_2$mesh
    degrees <- seq(from, to, length = iter)
    radians <- (degrees * pi)/180
    for (i in 1:length(degrees)) {
      set_3_t <- rotaxis3d(set_2, set_2[2,], c(500,0,0), radians[i])
      diff_align_z <- (set_3_t[3, 3] - set_3_t[4, 3])
      if (abs(diff_align_z) <= tol) {
        break
      }
    }
    sur_3 <- rotaxis3d(sur_2, set_2[2,], c(500,0,0), radians[i])
    set_3 <- t(sur_3$vb)[pos_1, -4]
    if((set_3[3,3] >0 & set_3[4,3] >0)==TRUE){
      sur_3 <- rotaxis3d(sur_3, set_3[2,], c(500,0,0), (180 * pi)/180)
      set_3 <- t(sur_3$vb)[pos_1, -4]
    }
    
    dist_1_z<-colMeans(set_3[c(3,4),])[3]
    dist_2_z<-set_3[1,3]
    sur_4<-sur_3
    sur_4$vb[3,]<-sur_4$vb[3,]+abs(dist_1_z)
    set_4 <- t(sur_4$vb)[pos_1, -4]
    out<-list("sur"=sur_4,"coo"=set_4)
  }
morphomapLanDia<-function(sur,side=c("left","right"),param1=4){
    sur<-morphomapSegm(sur)$external
    xaxis<-as.vector(range(sur$vb[1,]))
    seqs<-seq(xaxis[1],xaxis[2],length.out = 20)[c(4,15)]
    
    p1 <- c(seqs[1], 0, 0)
    p2 <- c(seqs[1], 100, 0)
    p3 <- c(seqs[1], 0, 100)
    normal <- crossProduct(p2 - p1, p3 - p1)
    zeroPro <- points2plane(rep(0,3),p1,normal)
    sig <- sign(crossprod(-zeroPro,normal))
    d <- sig*norm(zeroPro,"2")
    sect_t1<-meshPlaneIntersect(sur, p1, p2, p3)
    sect_t2 <- morphomapSort(sect_t1[,c(2,3)])
    sect_tp<-cbind(sect_t1[,1],sect_t2)
    points_1<-morphomapRegradius(sect_tp[,c(2,3)],n = 4,center=colMeans(sect_tp[,c(2,3)]))
    p1 <- c(seqs[2], 0, 0)
    p2 <- c(seqs[2], 100, 0)
    p3 <- c(seqs[2], 0, 100)
    normal <- crossProduct(p2 - p1, p3 - p1)
    zeroPro <- points2plane(rep(0,3),p1,normal)
    sig <- sign(crossprod(-zeroPro,normal))
    d <- sig*norm(zeroPro,"2")
    sect_t1<-meshPlaneIntersect(sur, p1, p2, p3)
    sect_t2 <- morphomapSort(sect_t1[,c(2,3)])
    sect_td<-cbind(sect_t1[,1],sect_t2)
    points_2<-morphomapRegradius(sect_td[,c(2,3)],n = 4,center=colMeans(sect_td[,c(2,3)]))
    if(side=="right"){
      sets<-rbind(sect_td[points_2[c(1,2)],],
                  sect_tp[points_1[c(1)],])
    }
    if(side=="left"){
      sets<-rbind(sect_td[points_2[c(3,2)],],
                  sect_tp[points_1[c(3)],])
    }
    return(sets)
  }
morphomapAlMesh<-function(set,mesh){
    eucl<-dist(set,method="euclidean")
    newP1<-c(0,0,0)
    newP2<-c(eucl[1],0,0)
    newP3<-c(((eucl[1]^2)+(eucl[2]^2)-(eucl[3]^2))/(2*eucl[1]), sqrt((eucl[2]^2)-(((eucl[1]^2)+(eucl[2]^2)- (eucl[3]^2))/(2*eucl[1]))^2), 0)
    newP3[which(is.na(newP3))]<-0
    tar<-rbind(newP1,newP2,newP3)
    rot_mesh<-rotmesh.onto(mesh,as.matrix(set),as.matrix(tar))
    return(rot_mesh)
  }
morphomapAlTwo<-function(mesh,set,iter,from,to,
                         tol, iter2, from2, 
                         to2){
    pos<-aro.clo.points(vert2points(mesh),set)$position
    set_2<-set
    for(j in 1:iter2){ 
      if((j>=2)==TRUE){from<-from2}
      if((j>=2)==TRUE){to<-to2}
      if((j>=2)==TRUE){tol<-tol}
      degrees <- seq(from, to, length = iter)
      degrees<-c(0,degrees)
      radians <- (degrees * pi)/180
      for (i in 1:length(degrees)) {
        set_2_t <- rotaxis3d(set_2, set_2[2,], c(0,500,500), radians[i])
        diff_align_z <- set_2_t[2,3]-set_2_t[5, 3] 
        if (abs(diff_align_z) <= tol) {
          break
        }
      }
      sur_3<-rotaxis3d(mesh, set_2[2,], c(0,500,500), radians[i])
      set_3<-t(sur_3$vb)[pos, -4]
      
      sur_4<-sur_3
      set_4<-t(sur_4$vb)[pos, -4]
      degrees <- seq(from, to, length = iter)
      degrees<-c(0,degrees)
      radians <- (degrees * pi)/180
      for (i in 1:length(degrees)) {
        set_4_t <- rotaxis3d(set_4, set_4[2,], c(0,0,500), radians[i])
        diff_align_z <- set_4_t[2,2]-set_4_t[6, 2] 
        if (abs(diff_align_z) <= tol) {
          break
        }
      }
      sur_5<-rotaxis3d(sur_4, set_4[2,], c(0,0,500), radians[i])
      set_5<-t(sur_5$vb)[pos, -4]
      
      degrees <- seq(from, to, length = iter)
      degrees<-c(0,degrees)
      radians <- (degrees * pi)/180
      for (i in 1:length(degrees)) {
        set_5_t <- rotaxis3d(set_5, set_5[2,], c(0,500,0), radians[i])
        diff_align_z <- set_5_t[5,3]-set_5_t[7,3] 
        if (abs(diff_align_z) <= 0.05) {
          break
        }
      }
      sur_6<-rotaxis3d(sur_5, set_5[2,], c(0,500,0), radians[i])
      set_6<-t(sur_6$vb)[pos, -4]
      
      degrees <- seq(from, to, length = iter)
      degrees<-c(0,degrees)
      radians <- (degrees * pi)/180
      for (i in 1:length(degrees)) {
        set_6_t <- rotaxis3d(set_6, set_6[2,], c(500,0,0), radians[i])
        diff_align_z <- set_6_t[6,2]-set_6_t[2,2] 
        if (abs(diff_align_z) <= 0.05) {
          break
        }
      }
      sur_7<-rotaxis3d(sur_6, set_6[2,], c(500,0,0), radians[i])
      set_7<-t(sur_7$vb)[pos, -4]
      
      sur_8<-sur_7
      sur_8$vb[1,]<-sur_7$vb[1,]-set_7[2,1]
      sur_8$vb[2,]<-sur_7$vb[2,]-set_7[2,2]
      sur_8$vb[3,]<-sur_7$vb[3,]-set_7[2,3]
      set_8<-t(sur_8$vb)[pos, -4]
      
      sur_9<-sur_8
      sur_9$vb[1,]<-sur_8$vb[1,]-mean(c(set_8[9,1],set_8[10,1]))
      set_9<-t(sur_9$vb)[pos,-4]
      
      ax1<-(abs(set_9[3,3]-set_9[4,3]))
      ax2<-(abs(set_9[7,3]-set_9[5,3]))
      ax3<-(abs(set_9[2,2]-set_9[6,2]))
      
      
      if((ax1<=tol & ax2 <= tol &ax3<= tol)==TRUE){
        break
      }
      
      if((ax1<=tol & ax2 <= tol &ax3<= tol)==FALSE){
        mesh<-sur_9
        set_2<-set_9
      }
    }
    
    mech_lengh<-set_9[8,1]-mean(c(set_9[c(9,10),1]))
    sur_9<-rotaxis3d(sur_9,c(0,0,0),c(0,500,0),(-90 * pi)/180)
    sur_9 <- rotaxis3d(sur_9, c(0,0,0), c(0,0,500), (270 * pi)/180)
    set_9 <- t(sur_9$vb)[pos, -4]
    out<-list("sur"=sur_9,"coo"=set_9,"mech_length"=mech_lengh)
    return(out)
  }
  
set0<-set
set1<-morphomapAlOne(mesh,set0[c(1,4,2,3),],iter=iter1,from=from1,to=to1,tol=tol)
set2<-morphomapLanDia(set1$sur,side,param1=param1)
posold<-aro.clo.points(vert2points(mesh),set0[c(1,4,2,3),])$position
set_new<-vert2points(set1$sur)[posold,]
posold1<-aro.clo.points(vert2points(mesh),set0[c(5:7),])$position
set_new1<-vert2points(set1$sur)[posold1,]
set_new1<-rbind(set2,set_new1)
setf<-morphomapAlTwo(set1$sur,set=rbind(set_new,set_new1),
                      iter=iter2,from=from2,to=to2,tol=tol,iter2=iter3,
                      from2=from3,to2=to3)
out<-list("sur"=setf$sur,"mech_length"=setf$mech_length)
return(out)
}

