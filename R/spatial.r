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
#' @param checkPairs logical; check for a single pair of outliers using the Grubbs test. This can only be performed for sample sizes <30. Only a single test is used because repeating it tends to throw out more points than seem reasonable, by eye. The value has no effect unless `method='grubbs'`.
#' @param kRosner integer between 1 and 10. Determines the number of outliers suspected with a Rosner test. The value has no effect unless `method='rosner'`.
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
                             checkPairs=TRUE,
                             kRosner=NULL){
  #  for testing
  #   pvalSet=1e-5; checkPairs=T
  
  pres.inliers=pres
  sp.toss.coord=NULL
  sp.toss.id=NULL
  
  if(any(method=='grubbs')){
    pval=0
    #tmp.dists=dists
    #toss singles
    while(pval<pvalSet & length(pres.inliers)>3){
    	# i want to recompute the distances once an outlier is removed because the outlier biased the centroid of the group, which influences distances
      dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
      gt=outliers::grubbs.test(dists)
      #gt=dixon.test(dists)
    	#cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
      (pval=gt$p.value)
      # conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
      if(is.na(pval)) {
      	warning('p-value for grubbs test was NA. sample size is too small to implement the test')
      	break
      }
      if(gt$p.value<pvalSet){
        toss=which.max(dists)
        # IDs in the original data frame
        sp.toss.coord=rbind(sp.toss.coord,sp::coordinates(pres.inliers)[toss,])
        pres.inliers=pres.inliers[-toss,]
        #tmp.dists=tmp.dists[-toss]
      }
    }
    if(length(dists)<3) warning('All but two records were deemed outliers. The Grubbs test may not be appropriate for these data.')
    # toss pairs
    if(checkPairs){
    	if(length(pres.inliers)<31 & length(pres.inliers)>3){
  			pval=0
  			#tmp.dists=dists
      	dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid

  			# By turning off this loop, I'm ensuring that you can only toss 1 pair of outliers. with the loop, it tends to find lots of supposed outliers very confidently, but by eye, it tends to omit clusters
  			#while(pval<pvalSet){
  				gt=outliers::grubbs.test(dists,type=20)
  
  				#gt=dixon.test(dists)
  				#cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
  				(pval=gt$p.value)
  				# conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
  				if(is.na(pval)) {
      			warning('p-value for grubbs test checking for pairs of outliers was NA. sample size is too small to implement the test')
      			break
      	}
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
    } 
  } # end grubbs
  
  if(any(method=='iqr')) {
  	dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
  	sp.toss.id=.iqrOutlier(dists)
  }
  
  if(any(method=='dixon')){
  	if(length(pres.inliers)<3 | length(pres.inliers)>30 ){
  		warning('dixon test only applies to sample sizes between [3,30]. skipping this taxon')
    	return(sp.toss.id)
  	}
    dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
    if(length(unique(dists))==1) {
    	warning('all records are the same distance from the spatial centroid so outliers cannot be detected. maybe your records come from gridded data. skipping this taxon')
    	return(sp.toss.id)
    }
    dt=outliers::dixon.test(dists,type=0,two.sided = FALSE)
    #cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
    if(dt$p.value<pvalSet) sp.toss.id=which.max(dists)
  }
  
  if(any(method=='rosner')){
    dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
    if(kRosner <= length(dists)) {
    	warning('kRosner must be an integer less than the number of presence records, skipping this taxon')
    	return(sp.toss.id)
    }
    rt=EnvStats::rosnerTest(dists,kRosner,alpha=pvalSet)
    if(any(rt$all.stats$Outlier)) sp.toss.id=utils::tail(order(dists),kRosner)
  }
  
  sp.toss.id
}