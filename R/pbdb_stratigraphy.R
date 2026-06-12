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
		x <- pbdb_assign_bin_interval(x, categs=e$keys$stgInt, out="stg")

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

	return(x)
}

#' Assign ten million year stratgiraphic bin values to PBDB occurrnces \code{stg}
#'
#' This protocol uses the \code{ten} assignments of \code{divDyn}
#'
#' The framework is described in detail 
#' @param x A Paleobiology Database occurrence data.frame (full output)
#' @param version A version number of the binning scheme.
#' @return An occurrence data.frame with the column \code{ten} added to it. 
#' @export
pbdb_strat_ten <- function(x, version="1.0"){

	# grab the keys from divDyn
	e <- environment()
	data(keys, package="divDyn", envir=e)

	# execute the main categorization, resolve interval names
	x <- pbdb_assign_bin_interval(x, categs=e$keys$tenInt, out="ten")

	return(x)
}



# General utily function for creating a binning as descrbied in the divDyn SOM (ddPhanero)
#
# @param x AN occurrence data frame.
# @param categs Category list, as the elements of divDyn keys$stgInt 
# @param out column name of the resulting stratigraphiy assignment
# @param early_interval column name of the early interval
# @param late_interval column name of the late interval
pbdb_assign_bin_interval <- function(x, categs, out, early_interval="early_interval", late_interval="late_interval"){

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

	# add to the rest and rename
	x<- cbind(x,res)
	colnames(x)[ncol(x)] <- out

	return(x)
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
