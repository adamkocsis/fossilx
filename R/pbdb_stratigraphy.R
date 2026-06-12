# stg:1.0

#' Assign stage-level stratgiraphy values to PBDB occurrnces \code{stg}
#'
#' This protocol uses the \code{stg} assignments of \code{divDyn}
#'
#' The framework is described in detail
#' @param x A Paleobiology Database occurrence data.frame (full output)
#' @param version A version number of the binning scheme.
#' @return An occurrence data.frame with the column \code{stg} added to it. 
#' @export
pbdb_strat_stg <- function(x, version="1.0"){
	# match the arguments
	version <- match.arg(version, "1.0")

	if(version=="1.0"){
		# grab the keys from divDyn
		e <- environment()
		data(keys, package="divDyn", envir=e)

		# execute the main categorization, resolve interval names
		x$stg <- pbdb_bin_interval(x, categs=e$keys$stgInt)

		# the Cambrian
		data(strat.camb_20180831,package="fossilx", envir=e)
		x <- replace_entries(x, where="stg", with=e$strat.camb_20180831, by="collection_no")

		# the Ordovician
		data(strat.ord_20180831,package="fossilx", envir=e)
		x <- pbdb_strat_ordovician_stg(
			dat=x,
			format=e$strat.ord_20180831$format,
			max.int=e$strat.ord_20180831$max.int,
			zones=e$strat.ord_20180831$zones
		)
	}else{
		stop("No other versions exist yet.")
	}

	return(x$stg)
}

#' Assign ten million year stratgiraphic bin values to PBDB occurrnces \code{stg}
#'
#' This protocol uses the \code{ten} assignments of \code{divDyn}
#'
#' The framework is described in detail 
#' @param x A Paleobiology Database occurrence data.frame (full output)
#' @param version A version number of the binning scheme.
#' @return A vector of bin assignments.
#' @export
pbdb_strat_ten <- function(x, version="1.0"){

	# grab the keys from divDyn
	e <- environment()
	data(keys, package="divDyn", envir=e)

	# execute the main categorization, resolve interval names
	bin <- pbdb_bin_interval(x, categs=e$keys$tenInt)

	return(bin)
}



# General utily function for creating a binning as descrbied in the divDyn SOM (ddPhanero)
#
# @param x AN occurrence data frame.
# @param categs Category list, as the elements of divDyn keys$stgInt 
# @param out column name of the resulting stratigraphiy assignment
# @param early_interval column name of the early interval
# @param late_interval column name of the late interval
pbdb_bin_interval <- function(x, categs, early_interval="early_interval", late_interval="late_interval"){

	# do the categorization
	binMin <- divDyn::categorize(x[,early_interval], categs)
	binMax <- divDyn::categorize(x[,late_interval], categs)

	# make them numeric
	binMin <- as.numeric(binMin)
	binMax <- as.numeric(binMax)

	# empty container
	res <- rep(NA, nrow(x))

	# select entries, where
	binCondition <- c(
	# the early and late interval fields indicate the same stg
		which(binMax==binMin),
	# or the late_intervarl field is empty
		which(binMax==-1))

	res[binCondition] <- binMin[binCondition]

	return(res)
}

#' Function to generate stage-lookup table from the Paleobiology Database occurrence data.frame
#'
#' Processess time-bin data downloaded from the database
#'
#' @param pbdb data.frame from database.
#' @param stagecolumn Charater string, of the time-bin information column
#' @param min_ma Character string, column name of the minimum ages 
#' @param max_ma Character string, column name of the maximum ages 
#' @param bin Character string, this will be the column name of the bin information
#' @return A data.frame with the stages' name bottom, mid and top age and bin number.
#' @export
pbdb_timebins <- function(pbdb, stagecolumn="time_contain", min_ma="min_ma", max_ma="max_ma", bin="stb"){
	# get the appropriate colulmns out of the PBDB download
	timeinfo <- unique(pbdb[,c(stagecolumn, max_ma, min_ma)])

	# omit those entires that are not meaningful
	timeinfo <- timeinfo[timeinfo[,stagecolumn]!="-",]
	
	# search for stage minima and maxima
	ma <- tapply(INDEX=timeinfo[,stagecolumn], X=timeinfo[, max_ma], max)
	min <- tapply(INDEX=timeinfo[,stagecolumn], X=timeinfo[, min_ma], min)

	# create stages object
	stagesDF <- data.frame(stage=names(ma), bottom=ma, mid=(ma+min)/2, top=min)
	stagesDF <-stagesDF[order(stagesDF$bottom, decreasing=TRUE),]
	rownames(stagesDF) <- NULL

	# check consistency - Cambrian is not consistent
	# stages$top[1:(nrow(stages)-1)]==stages$bottom[2:nrow(stages)]


	# create integer bin numbers
	stagesDF<- cbind(stagesDF, 1:nrow(stagesDF))
	colnames(stagesDF)[ncol(stagesDF)] <- bin
	return(stagesDF)
}


#' Function to bin PBDB data to stb stages 
#'
#' @param pbdb data.frame from database.
#' @param ts data.frame produced by StagesFromPBDB.
#' @param bin Character string, this will be the column name of the bin information
#' @param stagecolumn Charater string, of the time-bin information column in the PBDB
#' @return A vector with integer bin numbers
#' @examples
#' pbdb <- chronosphere::fetch("pbdb")
#' # new stages object
#' stages <- StagesFromPBDB(pbdb)
#' # bins from the pbdb output
#' pbdb$stb <- BinPBDB(pbdb, stages)
#' @export
pbdb_strat_stb <- function(pbdb,ts,  bin="stb",stagecolumn="time_contain"){
	
	# new vector
	y <- rep(NA, nrow(pbdb))

	# simple for loop to execute the binning
	for(i in 1:nrow(ts)){
		y[which(pbdb[, stagecolumn]==ts$stage[i])] <-ts[i,bin]
	}

	return(y)
}

# Script to assign Orodician collections to 'stg' stages based on tables format, max.int, and zones
# last checked with data of 2018-08-31 - Adam Kocsis
pbdb_strat_ordovician_stg <- function(dat, format, max.int,zones){
	# Normalize to collections
	new <- unique(dat[is.na(dat$stg), c(
		"collection_no", 
		"early_interval", 
		"late_interval", 
		"zone", 
		"formation",
		"max_ma",
		"min_ma",
		"stg")])
	
	# check always
	
	# Looping formation (added condition that period is same)
    for (i in 1:nrow(format))  {
		ix <- which((as.character(new$formation) == as.character(format$formation[i])))
		new$stg[ix] <- format$stg[i]
    }
	
#	x <- table(new$stg) # control
#   sum(x)

    # Looping early_intervals 
    for (i in 1:nrow(max.int))  {
        ix <- which(as.character(new$early_interval) == as.character(max.int$Max.int[i]))
        new$stg[ix] <- max.int$stg.1[i]
    }

    #  Looping late intervals (to check if different)
	stg2 <- rep(NA, nrow(new))
    for (i in 1:nrow(max.int))  {
		ix <- which(as.character(new$late_interval) == as.character(max.int$Max.int[i]))
        stg2[ix] <- max.int$stg.1[i]
    }

	ix <- which(new$stg<stg2) # should ignore NAs in second column
    new$stg[ix] <- NA


 #   x <- table(new$stg) # control
 #   sum(x)

	
	# Looping zones
	for (i in 1:nrow(zones))  {
		ix <- which(as.character(new$zone) == as.character(zones$zone[i]))
		new$stg[ix] <- zones$stg[i]
	}

	
#	x <- table(new$stg) # control
#	sum(x)
#	View(new[new$min_ma>400,])
	
	# only that part, which has stg assignments now
	new2 <- new[!is.na(new$stg), ]

	# vector: stg numbers, names:collection numbers
	ord <- new2$stg
	names(ord) <- new2$collection_no

	# the collection identifiers of occurrences in the total dataset
	colls <- as.character(dat$collection_no)

	# which are present in the newly gathered data?
	bool <- colls%in%names(ord)

	# collection identifiers of the occurrences of only these collections
	subColls <- colls[bool]
	
	# order/assign the stg accordingly
	subStg<-ord[subColls]
	
	# copy original
	newStg <- dat$stg
	
	# replace the missing entries
	newStg[bool]  <- subStg
	
	# make sure things are OK
#	origTab <- table(dat$stg)
#	newTab <- table(newStg)
#	newStg-origTab # should be all positive!!!

	# add to the full table
	dat$stg <- newStg

	return(dat)
}
