#' morphomapReadMorphologika
#' 
#' Import an array stored in a morphologika file
#' @param file path of the file to be read
#' @return out list containing an array, labels, groups and variables
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @export

morphomapReadMorphologika<-function(file){
  out<-list("array"=NULL,"labels"=NULL,"groups"=NULL,"variables"=NULL)
  allines<-as.vector(read.delim(file,header = FALSE)[,1])
  nind<-as.numeric(allines[grep("individuals",tolower(allines))+1])
  nlan<-as.numeric(allines[grep("landmarks",tolower(allines))+1])
  ndim<-as.numeric(allines[grep("dimensions",tolower(allines))+1])
  names_s<-grep("names",allines)
  coor<-grep("rawpoints",allines)+1
  groups_pos<-grep("group",allines)+1
  if(length(groups_pos)!=0){
    group_all<-allines[groups_pos:(names_s-1)]
    group_all<-strsplit(group_all," ")[[1]]
    group_names<-group_all[seq(1,length(group_all),2)]
    group_sizes<-group_all[seq(2,length(group_all),2)]
    groups<-c(rep(group_names,group_sizes))
    out$groups<-groups  
  }
  
  names<-allines[(names_s+1):(names_s+nind)]
  out$labels<-names  
  
  labvar<-grep("labels",allines)+1
  labval<-grep("labelvalues",allines)+1
  if(length(labvar)!=0){
    labelnames<-allines[c(labvar:(labval-2))]
    labelvalues<-allines[c(labval:(coor-2))]
    label_levels<-strsplit(labelnames," ")[[1]]
    label_values<-unlist(strsplit(labelvalues," "))
    
    list_vars<-list()
    for(i in 1:length(label_levels)){
      list_vars[[i]]<- label_values[seq(i,length(label_values),by=length(label_levels))]
    }
    names(list_vars)<-label_levels
    out$variables<-list_vars  
  }
  
  tempdata <-allines[coor:length(allines)]
  spec_pos<-which(substr(allines,1,1)=="'")
  alldataland<-allines[coor:(spec_pos[nind]+nlan)]
  alldataland <- alldataland[-which(substr(alldataland,1,1)=="'")]
  
  start<-spec_pos
  array<-array(NA,dim=c(nlan,ndim,nind))
  for(i in 1:length(start)){
    array[,,i]<-as.matrix(read.table(file,skip=start[i],nrows=nlan,dec="."))
  }
  out$array<-array  
  
  return(out)
}