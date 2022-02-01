#' @title Find outlying occurrence data
#'
#' @description Spatial or environmental outliers
#'
#' @details
#' See Examples.
#'
#' @param pres a `SpatialPoints` or `SpatialPointsDataFrame` object describing the locations of species records. A `SpatialPointsDataFrame` containing the values of environmental variables to be used must be supplied if `envOutliers=TRUE`
#' @param spOutliers logical; perform spatial outlier analysis
#' @param envOutliers logical; perform environmental outlier analysis
#' @param method character; options are 'iqr', 'grubbs', 'dixon', 'rosner'
#' @param pval user-specified p-value for assessing the significance of Grubbs test statistic.
#' @param checkPairs logical; check for a single pair of outliers. This can only be performed for sample sizes <30. Only a single test is used because repeating it tends to throw out more points than seem reasonable, by eye.
#' @param verbose logical; print messages
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
#' @return Returns the SpatialPointsDataFrame provided with additional logical columns indicating spatial outliers (`spOutliers`) and environmental outliers ('`envOutliers`).
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.


findOutlyingPoints=function(pres,
                            spOutliers=TRUE,
                            envOutliers=TRUE,
                            method='grubbs',
                            pval=1e-5,
                            checkPairs=TRUE,
                            verbose=TRUE){

  #  for testing
  #  pres=myPresDF; env=myEnv; verbose=T
  #  envOutliers=T; spOutliers=T; doPlot=T; plotFile=NULL; pval=1e-5; envVar=names(env)
  # #plotFile=paste0(dirs$sp.diag.path, '/', species,'__OccurrenceOutliers.png')

  if(!any(class(pres)==c('SpatialPoints','SpatialPointsDataFrame'))) stop('Please make your presence data a SpatialPoints or SpatialPointsDataFrame object and try again')
  
  if(spOutliers) { 
  	sp.toss.id=findSpatialOutliers(pres=pres,pvalSet=pval,checkPairs=checkPairs,method=method)
  	if(verbose) print(paste0(length(sp.toss.id),' geographic outlier(s) found'))
  } else {sp.toss.id=NULL}
  
  if(envOutliers) { 
  	env.toss.id=findEnvOutliers(pres=pres,pvalSet=pval,checkPairs=checkPairs,method=method)
  	pres$envOutlier=FALSE
  	pres$envOutlier[env.toss.id]=TRUE
  	if(verbose) print(paste0(length(env.toss.id),' environmental outlier(s) found'))
  	if(nrow(pres)-length(env.toss.id) < 2){
  		warning(paste0('pretty much all the presences came up as environmental outliers. The only case Ive seen this in was when there were two obvious outliers and all the other presences had the same exact environment. In any case, I cant find any outliers.'))
  	}
  }

  if(spOutliers) {
    pres$spOutlier=FALSE
    pres$spOutlier[sp.toss.id]=TRUE
  }

  return(pres)
}
