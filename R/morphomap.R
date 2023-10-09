#' @name morphomap-package
#' @docType package
#' @aliases morphomap
#' @title 2D and 3D cortical thickness maps and cross sectional geometry
#' @author Antonio Profico, Luca Bondioli, Pasquale Raia, Paul O'Higgins, Damiano Marchi
#' @description Tool to process long bone meshes (shape data, morphometric maps and cross-sectional geometry)
#' @import Arothron
#' @import lattice
#' @import mgcv
#' @import Rvcg
#' @import Morpho
#' @import geometry
#' @import colorRamps
#' @import DescTools
#' @import grDevices
#' @import graphics
#' @import utils
#' @importFrom colorRamps blue2green2red
#' @importFrom grDevices dev.off tiff
#' @importFrom rgl open3d layout3d triangles3d next3d rotate3d rotationMatrix wire3d axis3d title3d bbox3d par3d view3d planes3d spheres3d lines3d
#' @importFrom oce matrixSmooth
#' @importFrom sp point.in.polygon CRS
#' @importFrom graphics plot points polygon text
#' @importFrom stats IQR quantile sd weighted.mean dist prcomp
#' @importFrom utils setTxtProgressBar txtProgressBar
NULL

