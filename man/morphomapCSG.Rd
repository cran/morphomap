% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morphomapCSG.R
\name{morphomapCSG}
\alias{morphomapCSG}
\title{morphomapCSG}
\usage{
morphomapCSG(
  cp,
  mp,
  translate = FALSE,
  center = c("I", "E", "CCA"),
  delta = 0.1,
  Cx = NULL,
  Cy = NULL,
  I_xy = TRUE,
  I_minmax = TRUE,
  Zxy = TRUE
)
}
\arguments{
\item{cp}{matrix: coordinates of the external outline}

\item{mp}{matrix: coordinates of the internal outline}

\item{translate}{logical: if TRUE the section will be centered}

\item{center}{how to define the center of each section. The method allowed are "CCA" (center of cortical area), "E" (barycenter of the external outline) and "I" (barycenter of the internal outline)}

\item{delta}{numeric: picture elements of adjustable side length}

\item{Cx}{numeric: new x center coordinate}

\item{Cy}{numeric: new y center coordinate}

\item{I_xy}{logical: if TRUE the product of inertia around the x and y axis is calculated}

\item{I_minmax}{logical: if TRUE the Imin and Imax will be calculated}

\item{Zxy}{logical: if TRUE the polar moment of inertia will be calculated}
}
\value{
Cx x coordinate of the centered section

Cy y coordinate of the centered section

T_area total area

M_area medullar area

CA cortical area

Ext_perim external perimeter

Med_perim medullar perimiter

Mean_thick mean thickness of the section

Sd_thick thickness standard deviation

Min_thick minimum thickness

Max_thick maximum thickness

Ix numeric: moment of inertia around the x axis

Iy numeric: moment of inertia around the y axis

Zx numeric: moment of inertia around the x axis

Zy numeric: moment of inertia around the y axis

Zpol numeric: polar moment of inertia

dx new centered coordinates of the internal outline

dy new centered coordinates of the internal outline

Imin numeric: minimum moment of inertia

Imax numeric: maximum moment of inertia

J numeric: polar moment of inertia

Zmax numeric: the maximum polar section

Zmin numeric: the minimum polar section

theta numeric: theta angle
}
\description{
Tool for Cross-sectional geometry
}
\examples{
\donttest{
#calculation of csg parameter on a human femur cross section
data(HomFem38023)
meshes<-morphomapSegm(HomFem38023, param1=4)
perMesh<-meshes$external
endMesh<-meshes$internal
mech_length<-380.23
rawSections<-morphomapCore(out.sur=perMesh,
                           inn.sur=endMesh,num.sect=61,mech.len = mech_length, 
                           start = 0.2,end=0.8)
shapeSections<-morphomapShape(rawSections,250,sects_vector=NULL,cent.out="CCA",
delta=0.1, side="left")
csgSect31<-morphomapCSG(cp = shapeSections$`2D_out`[,,31],
                        mp=shapeSections$`2D_inn`[,,31],
                        translate = FALSE,center="CCA")

#Cross sectional geometry along the entire femur bone
results<-matrix(NA,ncol=24,nrow=61)
rownames(results)<-paste("section",c(1:61))
colnames(results)<-c("Cx","Cy","T_area","M_area","CA",
                     "Ext_perim","Med_perim","Mean_thick","Sd_thick" ,
                     "Min_thick","Max_thick","Ix","Iy","Zx" ,"Zy","Zpol" ,
                     "dx","dy","Imin","Imax","J","Zmax","Zmin","theta")

for(i in 1:61){
  results[i,]<-unlist(morphomapCSG(cp = shapeSections$`2D_out`[,,i],
                                     mp=shapeSections$`2D_inn`[,,i],
                                   translate = FALSE,center="CCA",delta = 0.5))
}

plot(c(1:61),results[,24],type="b",main="Theta",cex=1,
     xlab="section",ylab="radians")


#calculation of csg parameter on a chimpanzee femur cross section
data(PanFem27713)
meshes<-morphomapSegm(PanFem27713, param1=3)
perMesh<-meshes$external
endMesh<-meshes$internal
mech_length<-277.13
rawSections<-morphomapCore(out.sur=perMesh,
                           inn.sur=endMesh,num.sect=61,mech.len = mech_length, 
                           start = 0.2,end=0.8)
shapeSections<-morphomapShape(rawSections,250,sects_vector=NULL,cent.out="CCA",
delta=0.1, side="left")
csgSect31<-morphomapCSG(cp = shapeSections$`2D_out`[,,31],
                        mp=shapeSections$`2D_inn`[,,31],
                        translate = FALSE,center="CCA")

#Cross sectional geometry along the entire femur bone
results<-matrix(NA,ncol=24,nrow=61)
rownames(results)<-paste("section",c(1:61))
colnames(results)<-c("Cx","Cy","T_area","M_area","CA",
                     "Ext_perim","Med_perim","Mean_thick","Sd_thick" ,
                     "Min_thick","Max_thick","Ix","Iy","Zx" ,"Zy","Zpol" ,
                     "dx","dy","Imin","Imax","J","Zmax","Zmin","theta")

for(i in 1:61){
  results[i,]<-unlist(morphomapCSG(cp = shapeSections$`2D_out`[,,i],
                                   mp=shapeSections$`2D_inn`[,,i],
                                   translate = FALSE,center="CCA",delta = 0.5))
  }

plot(c(1:61),results[,24],type="b",main="Theta",cex=1,
     xlab="section",ylab="radians")
}
}
\author{
Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
}
