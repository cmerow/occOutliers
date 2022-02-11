
#======================================================================
#======================================================================
#' @title Find outlying occurrence data in environmental space
#'
#' @description Environmental outliers
#'
#' @details
#' See Examples.
#'
#' @param pres a `SpatialPoints` or `SpatialPointsDataFrame` object describing the locations of species records. A `SpatialPointsDataFrame` containing the values of environmental variables to be used must be supplied if `envOutliers=TRUE`
#' @param method character; options are 'iqr', 'grubbs', 'dixon', 'rosner'
#' @param pvalSet user-specified p-value for assessing the significance of Grubbs test statistic.
#' @param checkPairs logical; check for a single pair of outliers using the Grubbs test. This can only be performed for sample sizes <30. Only a single test is used because repeating it tends to throw out more points than seem reasonable, by eye. The value has no effect unless `method='grubbs'`.
#' @param kRosner integer between 1 and 10. Determines the number of outliers suspected with a Rosner test. The value has no effect unless `method='rosner'`.# @keywords
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
#' presOut=findEnvOutliers(pres=myPresDF,pvalSet=1e-5)
#' @return Returns the indices of environmental outliers
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

findEnvOutliers=function(pres,
                         #myEnv=NULL,
                         pvalSet=1e-5,
                         method='grubbs',
                         checkPairs=TRUE,
                         kRosner=NULL){
  #  for testing
  #  pres=presDF; pvalSet=1e-5; checkPairs=F; myEnv=NULL
  #  myEnv=env
  #  env=myEnv; pvalSet=1e-5
  
  #if(!is.null(myEnv)){ p.env=raster::extract(myEnv,pres)
  #} else {p.env=pres}
  p.env=pres@data
  p.env=base::scale(p.env)
  # remove variables that are the same for all observations
  f=which(apply(p.env,2,function(x) !all(is.nan(x))))
  p.env=p.env[,f]
  pres.inliers=p.env
  row.id=apply( p.env, 1 , paste , collapse = "-" )
  env.toss.id=NULL
  #dists=apply(p.env,1,function(x) sqrt(sum((x)^2)) )
  
  if(any(method=='grubbs')){
    pval=0
    #tmp.dists=dists
    while(pval<pvalSet & length(pres.inliers)>3){
    	# do this within the loop to recompute the centroid ofter each outlier is removed
      dists=apply(p.env,1,function(x) sqrt(sum((x)^2)) )
    	# this is close enough to fixed, in the outer function it'll report that everything was tossed
    	if(nrow(p.env)<3) break 
      gt=outliers::grubbs.test(dists)
      pval=gt$p.value
      # conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
      
      if(gt$p.value<pvalSet){
        toss=which.max(dists)
        # IDs in the original data frame
        thisID=paste(p.env[toss,],collapse='-')
        env.toss.id=c(env.toss.id,which(row.id == thisID))
        #tmp.dists=tmp.dists[-toss]
        p.env=p.env[-toss,]
      }  
    }
    
    #if(checkPairs) print('checkPairs not yet implemented')
    # toss pairs
    if(checkPairs){ # taken from spatial - must test
      if(length(pres.inliers)<31 & length(pres.inliers)>3){
        pval=0
			  dists=apply(p.env,1,function(x) sqrt(sum((x)^2)) )
        # By turning off this loop, I'm ensuring that you can only toss 1 pair of outliers. with the loop, it tends to find lots of supposed outliers very confidently, but by eye, it tends to omit clusters
        #while(pval<pvalSet){
        # duplicated
        #dists=.presPercentile(pres.inliers, percent=NULL)[[1]]$dist.from.centroid
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
  }# end grubbs
  
  if(any(method=='iqr')){
    dists=apply(p.env,1,function(x) sqrt(sum((x)^2)) )
    env.toss.id=.iqrOutlier(dists)
  }
  
  if(any(method=='dixon')){
    if(length(pres.inliers)<3 | length(pres.inliers)>30 ){
  		warning('dixon test only applies to sample sizes between [3,30]. skipping this taxon')
    	return(env.toss.id)
  	}
    dists=apply(p.env,1,function(x) sqrt(sum((x)^2)) )
    if(length(unique(dists))==1) {
    	warning('all records are the same distance from the environmental centroid so outliers cannot be detected. maybe your records come from gridded data. skipping this taxon')
    	return(env.toss.id)
    }
    dt=outliers::dixon.test(dists,type=0,two.sided = FALSE)
    #cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
    if(dt$p.value<pvalSet) env.toss.id=which.max(dists)
  }
  
  if(any(method=='rosner')){
    dists=apply(p.env,1,function(x) sqrt(sum((x)^2)) )
    if(kRosner <= length(dists)) {
    	warning('kRosner must be an integer less than the number of presence records, skipping this taxon')
    	return(env.toss.id)
    }
    rt=EnvStats::rosnerTest(dists,kRosner,alpha=pvalSet)
    if(any(rt$all.stats$Outlier)) env.toss.id=utils::tail(order(dists),kRosner)
  }
  
  # if(!is.null(sp.toss.coord)){ is this needed from spatial?
  #   coor=coordinates(pres)
  #   sp.toss.id= apply(sp.toss.coord,1,function(x) which(x[1]==coor[,1] & x[2]==coor[,2]))
  # } else {sp.toss.id=NULL}
  
  unique(env.toss.id)
}
