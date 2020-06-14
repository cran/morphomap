#' morphomapPic
#' 
#' Save the sections defined via morphomapShape or morphomapCore 
#' @param morphomap.core list: morphomap.core object
#' @param morphomap.shape list: morphomap.shape object
#' @param vector numeric: define which sections will be saved
#' @param full logical: if TRUE the thickness at ALPM is reported
#' @param width numeric: width of the picture
#' @param height numeric: height of the picture
#' @param pointsize numeric: pointsize of plotted text
#' @param res numeric: the nominal resolution in ppi which will be recorded 
#' @param colthk specify the color for the numbers
#' @param collbs specify the color for the labels
#' @param dirpath character: path of the directory where the pictures will be saved
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @examples
#' \donttest{
#' #export picture from a human femur bone
#' data(HomFem38023)
#' meshes<-morphomapSegm(HomFem38023, param1=4)
#' perMesh<-meshes$external
#' endMesh<-meshes$internal
#' mech_length<-380.23
#' rawSections<-morphomapCore(out.sur=perMesh,
#'                            inn.sur=endMesh,num.sect=11,mech.len = mech_length, 
#'                            start = 0.2,end=0.8)
#' shapeSections<-morphomapShape(rawSections,250,sects_vector=NULL,cent.out="CCA",
#' delta=0.5, side="left")
#' morphomapPic(rawSections,shapeSections,full=TRUE,dirpath=tempdir(),
#'             width=2500,height=2500)
#' 
#' #export picture from a chimpanzee femur bone
#' data(PanFem27713)
#' meshes<-morphomapSegm(PanFem27713, param1=3)
#' perMesh<-meshes$external
#' endMesh<-meshes$internal
#' mech_length<-277.13
#' rawSections<-morphomapCore(out.sur=perMesh,
#'                            inn.sur=endMesh,num.sect=11,mech.len = mech_length, 
#'                            start = 0.2,end=0.8)
#' shapeSections<-morphomapShape(rawSections,250,sects_vector=NULL,cent.out="CCA",delta=0.5,
#'  side="left")
#' morphomapPic(rawSections,shapeSections,full=TRUE,dirpath=tempdir(),
#'              width=2500,height=2500)
#' }
#' @export

morphomapPic<-function (morphomap.core,morphomap.shape, vector = NULL, full = TRUE,
                        width=1500,height=1500,
                        pointsize = 12, res=300,colthk="red",collbs="blue", dirpath = tempdir()) 
{
  out.lines <- morphomap.core$`2D_out`
  in.lines <- morphomap.core$`2D_inn`
  alpm_o <- morphomap.shape$ALPM_out
  alpm_i <- morphomap.shape$ALPM_inn
  if (is.null(vector)) {
    vector <- c(1:dim(out.lines)[3])
  }
  if (full == FALSE) {
    for (i in vector) {
      xlim <- range(out.lines[, 1, vector])
      ylim <- range(out.lines[, 2, vector])
      ylim[1] <- ylim[1] - 2
      xlim[1] <- xlim[1] - 2
      tiff(paste(dirpath,"/" ,paste("section number ", vector[i], ".tiff", 
                            sep = ""), sep = "/"), width = width, height = height, 
           pointsize = pointsize, res = res, compression = "lzw")
      graphics::plot(out.lines[, , i], asp = 1, type = "n", pch = 19, 
           cex = 3, main = paste("section number ", i, sep = ""), 
           xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
           xlim = xlim, ylim = ylim)
      polygon(out.lines[, , i], col = "black")
      polygon(in.lines[, , i], col = "white")
      xy_p <- rbind(c((xlim[1] + abs(xlim[1] * 0.04)), 
                      (ylim[1] + abs(ylim[1] * 0.025))), c((xlim[1] + 
                                                              (abs(xlim[1] * 0.04) + 10)), (ylim[1] + abs(ylim[1] * 
                                                                                                            0.025))))
      points(xy_p, type = "l", lwd = 3)
      # text(mean(xy_p[, 1]), mean(xy_p[, 2]), labels = "1 cm", 
      #      pos = 3)
      dev.off()
    }
  }
  if (full == TRUE) {
    for (i in vector) {
      diffs <- round(sqrt(rowSums(alpm_o[, , i] - alpm_i[, 
                                                         , i])^2), 2)
      xlim <- range(out.lines[, 1, vector])
      ylim <- range(out.lines[, 2, vector])
      ylim[1] <- ylim[1] - 2
      xlim[1] <- xlim[1] - 2
      tiff(paste(dirpath,"/" ,paste("section number ", i, ".tiff", 
                            sep = ""), sep = "/"), width = width, height = height, 
           pointsize = pointsize, res = res, compression = "lzw")
      graphics::plot(out.lines[, , i], asp = 1, type = "n", pch = 19, 
           cex = 3, main = paste("section number ", i, sep = ""), 
           xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
           xlim = xlim, ylim = ylim)
      polygon(out.lines[, , i], col = "black")
      polygon(in.lines[, , i], col = "white")
      points(rbind(alpm_o[1, , i], alpm_i[1, , i]), type = "l", 
             col = colthk, lwd = 2)
      text(alpm_o[1, 1, i], alpm_o[1, 2, i], labels = diffs[1], 
           col = colthk, pos = 3)
      text(mean(alpm_o[1, 1, i], alpm_i[1, 1, i]), mean(alpm_o[1, 
                                                               2, i], alpm_i[1, 2, i]) * 0.95, labels = "L", 
           col = collbs, pos = 2, lwd = 4, cex = 1.5)
      points(rbind(alpm_o[2, , i], alpm_i[2, , i]), type = "l", 
             col = colthk, lwd = 2)
      text(alpm_o[2, 1, i], alpm_o[2, 2, i], labels = diffs[2], 
           col = colthk, pos = 3)
      text(mean(alpm_o[2, 1, i], alpm_i[2, 1, i]), mean(alpm_o[2, 
                                                               2, i], alpm_i[2, 2, i]) * 0.95, labels = "A", 
           col = collbs, pos = 2, lwd = 4, cex = 1.5)
      points(rbind(alpm_o[3, , i], alpm_i[3, , i]), type = "l", 
             col = colthk, lwd = 2)
      text(alpm_o[3, 1, i], alpm_o[3, 2, i], labels = diffs[3], 
           col = colthk, pos = 3)
      text(mean(alpm_o[3, 1, i], alpm_i[3, 1, i]), mean(alpm_o[3, 
                                                               2, i], alpm_i[3, 2, i]) * 0.95, labels = "M", 
           col = collbs, pos = 4, lwd = 4, cex = 1.5)
      points(rbind(alpm_o[4, , i], alpm_i[4, , i]), type = "l", 
             col = colthk, lwd = 2)
      text(alpm_o[4, 1, i], alpm_o[4, 2, i], labels = diffs[4], 
           col = colthk, pos = 1)
      text(mean(alpm_o[4, 1, i], alpm_i[4, 1, i]), mean(alpm_o[4, 
                                                               2, i], alpm_i[4, 2, i]) * 1.05, labels = "P", 
           col = collbs, pos = 4, lwd = 4, cex = 1.5)
      xy_p <- rbind(c((xlim[1] + abs(xlim[1] * 0.04)), 
                      (ylim[1] + abs(ylim[1] * 0.025))), c((xlim[1] + 
                                                              (abs(xlim[1] * 0.04) + 10)), (ylim[1] + abs(ylim[1] * 
                                                                                                            0.025))))
      points(xy_p, type = "l", lwd = 3)
      text(mean(xy_p[, 1]), mean(xy_p[, 2]), labels = "1 cm", 
           pos = 3)
      dev.off()
    }
  }
}