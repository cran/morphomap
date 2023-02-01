#' morphomapExport
#'
#' Export the output from ToothAlignement 
#' @param mpShapeObject list: list containing morphomapShape objects
#' @param id character: label name
#' @param file character: name the output file
#' @author Antonio Profico
#' @export

morphomapExport<-function(mpShapeObject,id,file){
  cat(paste("# morphomap ","version ",packageVersion("morphomap"), "\n", 
            sep = ""), file = file, append = FALSE, sep = "")
  cat(paste("# ",id, "\n", 
            sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("# ","List of contents", "\r\n", 
            sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("3D_out ", "{ 3D coordinates}", " @1", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("3D_inn ", "{ 3D coordinates}", " @2", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("2D_out ", "{ 2D coordinates}", " @3", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("2D_inn ", "{ 2D coordinates}", " @4", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("ALPM_inn ", "{ 2D coordinates}", " @5", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("ALPM_out ", "{ 2D coordinates}", " @6", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("start ", "{ numeric }", " @7", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("end ", "{ 3D coordinates }", " @8", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("mech.len ", "{ numeric }", " @9", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("num.sects ", "{ numeric }", " @10", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("num.lands ", "{ numeric }", " @11", "\n", sep = ""), file = file, append = TRUE, sep = "")
  
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("### Data ###", "\n", sep = ""), file = file, append = TRUE, sep = "")
  coo3Dout<-morphomapArray2matrix(mpShapeObject$`3D_out`)
  coo3Dinn<-morphomapArray2matrix(mpShapeObject$`3D_inn`)
  coo2Dout<-morphomapArray2matrix(mpShapeObject$`2D_out`)
  coo2Dinn<-morphomapArray2matrix(mpShapeObject$`2D_inn`)
  cooALPM_inn<-morphomapArray2matrix(mpShapeObject$ALPM_inn)
  cooALPM_out<-morphomapArray2matrix(mpShapeObject$ALPM_out)
  coostart<-mpShapeObject$start
  cooend<-mpShapeObject$end
  coomlen<-mpShapeObject$mech.len
  num.sects<-dim(mpShapeObject$`3D_out`)[3]
  num.lands<-dim(mpShapeObject$`3D_out`)[1]
  
  cat(paste("@1", "\n", sep = ""), file = file, append = TRUE,sep = "")
  write.table(format(coo3Dout, scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@2", "\n", sep = ""), file = file, append = TRUE, sep = "")
  write.table(format(coo3Dinn, scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@3", "\n", sep = ""), file = file, append = TRUE, sep = "")
  write.table(format(coo2Dout, scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@4", "\n", sep = ""), file = file, append = TRUE, sep = "")
  write.table(format(coo2Dinn, scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@5", "\n", sep = ""), file = file, append = TRUE, sep = "")
  write.table(format(cooALPM_inn, scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@6", "\n", sep = ""), file = file, append = TRUE, sep = "")
  write.table(format(cooALPM_out, scientific = F, trim = T), 
              file = file, sep = " ", append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, na = "", 
              eol = "\n")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@7", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste(coostart, "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@8", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste(cooend, "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@9", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste(coomlen, "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@10", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste(num.sects, "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("@11", "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste(num.lands, "\n", sep = ""), file = file, append = TRUE, sep = "")
  cat(paste("\r", sep = ""), file = file, append = TRUE, sep = "")
  
}