#' morphomapSegm
#'
#' Separate a mesh from its visible and not visible components by using CA-LSE method
#' @param mesh object of class mesh3d
#' @param views numeric: number of points of view
#' @param param1 numeric: first parameter for spherical flipping
#' @param num.cores numeric: number of cores
#' @return external mesh3d of the visible facets from the points of view
#' @return internal mesh3d of the not visible facets from the points of view
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @references Profico A., Schlager S., Valoriani V., Buzi C., Melchionna M., Veneziano A.,
#' Raia P., Moggi-Cecchi J. and Manzi G., 2018. Reproducing the internal and external anatomy of fossil bones: Two new automatic digital tools.
#' American Journal of Physical Anthropology 166(4): 979-986.
#' @examples
#' \donttest{
#' #automatic separation of external and medullar femur components
#' library(rgl)
#' data(HomFem38023)
#' meshes<-morphomapSegm(HomFem38023)
#' perMesh<-meshes$external
#' endMesh<-meshes$internal
#' open3d()
#' wire3d(perMesh,col="grey")
#' wire3d(endMesh,col="red")
#' }
#' @export

morphomapSegm<-function(mesh,views=30,param1=3,num.cores=NULL){

  ca_lse<-ext.int.mesh(mesh, views=views, param1=param1,
                       default=TRUE,import_pov = NULL,
                       num.cores = num.cores)
  out_inn<-out.inn.mesh(ca_lse,mesh,plot=FALSE)

  inner_sur<-vcgIsolated(out_inn$invisible)
  outer_sur<-out_inn$visible

  out<-list("external"=outer_sur,"internal"=inner_sur)

}
