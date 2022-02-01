#======================================================================
#======================================================================
#' @title Find outlying occurrence data in geographic space
#'
#' @description Spatial outliers
#'
#' @details
#' See Examples.
#' @param pres a `SpatialPoints` or `SpatialPointsDataFrame` object describing the locations of species records. A `SpatialPointsDataFrame` containing the values of environmental variables to be used must be supplied if `envOutliers=TRUE`
#' @param method character; options are 'iqr', 'grubbs', 'dixon', 'rosner'
#' @param pvalSet user-specified p-value for assessing the significance of Grubbs test statistic.
#' @param checkPairs logical; check for a single pair of outliers using the Grubbs test. This can only be performed for sample sizes <30. Only a single test is used because repeating it tends to throw out more points than seem reasonable, by eye.
# @keywords
#' @export
#'
#' @examples
#' myPres=read.csv(system.file('extdata/SpeciesCSVs/Camissonia_tanacetifolia.csv',
#'                             package='occOutliers'))
#' myPres=myPres[complete.cases(myPres),]
#' sp::coordinates(myPres)=c(1,2)
#' presOut=findSpatialOutliers(pres=myPres, pvalSet=1e-5)
#' 
#' @return Returns the indices of spatial outliers
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

findSpatialOutliers=function(pres,
                             pvalSet=1e-5,
                             method='grubbs',
                             checkPairs=TRUE){
  #  for testing
  #  pvalSet=1e-5
  #  pres=mcp.pres; pvalSet=1e-5; checkPairs=T
  
  pres.inliers=pres
  sp.toss.coord=NULL
  
  if(any(method=='grubbs')){
    pval=0
    #toss singles
    while(pval<pvalSet){
    	dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
      gt=outliers::grubbs.test(dists)
      #gt=dixon.test(dists)
    	#cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
      (pval=gt$p.value)
      # conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
      if(gt$p.value<pvalSet){
        toss=which.max(dists)
        # IDs in the original data frame
        sp.toss.coord=rbind(sp.toss.coord,sp::coordinates(pres.inliers)[toss,])
        pres.inliers=pres.inliers[-toss,]
        dists=dists[-toss]
      }
    }
    # toss pairs
    if(checkPairs){
    	if(length(pres.inliers)<31){
  			pval=0
  			# By turning off this loop, I'm ensuring that you can only toss 1 pair of outliers. with the loop, it tends to find lots of supposed outliers very confidently, but by eye, it tends to omit clusters
  			#while(pval<pvalSet){
  				dists=.presPercentile(pres.inliers, percent=NULL)[[1]]$dist.from.centroid
  				gt=outliers::grubbs.test(dists,type=20)
  
  				#gt=dixon.test(dists)
  				#cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
  				(pval=gt$p.value)
  				# conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
  				if(gt$p.value<pvalSet){
  					toss=utils::tail(order(dists),2)
  					# IDs in the original data frame
  					sp.toss.coord=rbind(sp.toss.coord, sp::coordinates(pres.inliers)[toss,])
  					pres.inliers=pres.inliers[-toss,]
  				}
  			#}
  		}	
    }
    if(!is.null(sp.toss.coord)){
      coor=sp::coordinates(pres)
      sp.toss.id= apply(sp.toss.coord,1,function(x) which(x[1]==coor[,1] & x[2]==coor[,2]))
    } else {sp.toss.id=NULL}
  } # end grubbs
  
  if(any(method=='iqr')){
    dists=.presPercentile(pres.inliers, percent=NULL)[[1]]$dist.from.centroid
    sp.toss.id=.iqrOutlier(dists)
  }
  
  if(any(method=='dixon')){
    
  }
  
  if(any(method=='rosner')){
    
  }
  

  
  
  sp.toss.id
}