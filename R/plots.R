#' @title Find outlying occurrence data
#'
#' @description Spatial or environmental outliers
#'
#' @details
#' See Examples.
#'
#' @param pres a `SpatialPoints` or `SpatialPointsDataFrame` object describing the locations of species records. A `SpatialPointsDataFrame` containing the values of environmental variables to be used must be supplied if `envOutliers=TRUE`
#' @param env an optional raster object to be plotted
#' @param outlierNames character vector indicating the name of the columns in `pres` to plot. e.g., c('spOutlier','envOutlier')
#' @param shpToPlot an optional spatialPolygons object to be plotted
#' @param legLoc legend location, of the form 'bottomleft', 'topright', etc.
# @keywords
#' @export
#'
#' @examples
#' myPres=read.csv(system.file('extdata/SpeciesCSVs/Camissonia_tanacetifolia.csv',
#'                             package='occOutliers'))
#' myPres=myPres[complete.cases(myPres),]
#' sp::coordinates(myPres)=c(1,2)
#' myEnv=raster::stack(system.file('extdata/AllEnv.tif',package='occOutliers'))
#' names(myEnv)=read.table(system.file('extdata/layerNames.csv',package='occOutliers'))[,1]
#' myPresDF=sp::SpatialPointsDataFrame(myPres,data.frame(raster::extract(myEnv,myPres)))
#' presOut=findOutlyingPoints(pres=myPresDF,
#'                            spOutliers=TRUE,
#'                            envOutliers=TRUE,
#'                            pval=1e-5)
#' world.shp=readRDS(system.file('extdata/worldShpMollwide.rds',package='occOutliers'))
#' plotOutliers(presOut,shpToPlot = world.shp)
#' 
#' @return nothing
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

# will need to specify which cols to plot
plotOutliers=function(pres,env=NULL,outlierNames=c('spOutlier','envOutlier'),
                      shpToPlot=NULL,legLoc='bottomleft'){
  #  for testing
  #  pres=presOut; env=NULL; shpToPlot=NULL;plotFile=NULL; outlierNames=c('spOutlier','envOutlier')
  
  #if(!is.null(plotFile)) png(plotFile,h=1000,w=1000)
  pe=raster::extent(pres)
  ydif=pe[4]-pe[3] ; xdif=pe[2]-pe[1] # for plotting
  pe[1]=pe[1]-xdif; pe[2]=pe[2]+xdif; pe[3]=pe[3]-ydif; pe[4]=pe[4]+ydif
  
  if(length(outlierNames)==1){ toss=which(pres@data[,outlierNames])
  } else { toss=which(apply(pres@data[,outlierNames],1,any))}
  if(length(toss)>1) { pres0=pres[-toss,]
  } else {pres0=pres}
  if(!is.null(env)){
    env.bg=raster::crop(env[[1]],pe)
    raster::plot(env.bg,col='grey50',legend=FALSE)
    graphics::points(pres0,pch=3,cex=1.5,col='black')
  } else {
    sp::plot(pres0,pch=3,cex=1.5,col='black',xlim=c(pe[1],pe[2]),ylim=c(pe[3],pe[4]))
  }
  if(!is.null(pres$spOutlier)){
    graphics::points(pres[pres$spOutlier,],pch=16,cex=1.3,col='red3')
  }
  if(!is.null(pres$envOutlier)){
    graphics::points(pres[pres$envOutlier,],pch=16,cex=.9,col='steelblue3')
  }
  if(!is.null(shpToPlot)) sp::plot(shpToPlot,add=TRUE,lwd=.7,border='grey40')
  
  graphics::legend(legLoc,legend=c('good presence','spatial outlier','environmental outlier'),
         pch=c(3,16,16), col=c('black','red3','steelblue3'),bty='n',cex=1,pt.cex=1)
  #if(!is.null(plotFile)) dev.off()

}
