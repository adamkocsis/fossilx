#' Extend and post-process Paleobiology Database records
#'
#' Declarative definition of post-processing steps of Paleobiology Database data
#' 
#' @param include Multiple element mean successive filtering, i.e. the intersection of sets. Character entries represent pre-defined datasets.
#' @param tax.level Taxonomic resolution. Currently only the genus level is defined. Arguments passed to pbdb_taxon_quality_filter.
#' @param tax.combine Combine columsn for homonomy evasion?
#' @param env.categories Should environment categories be defined?
#' @param strat Stratigraphic assignment definitions. Available are 'stg', 'ten' and 'stb', colons separate differnet versions.
#' @param omit What records should be omitted from the set? 
#' @export
pbdb_extend <- function(dat,
	include=list("tax.marine_1.0", "env.marine_1.0"),
	tax.level=c("genus"),
	tax.combine="clgen",
	env.categories="divDyn",
	strat=c("stg:1.0", "ten:1.0", "stb"),
	omit=list("lithification1:unlithified"), verbose=TRUE){

	# quality filters
	if(!is.null(tax.level)){
		cat("Taxon quality filtering...                                           \r")
		flush.console()
		# omit pooor quality occs
		dat <- pbdb_taxon_quality_filter(dat, level=tax.level)
		cat("Taxon quality filtering...OK                                           \n")
		flush.console()
	}

	# tax.combine
	if(is.null(tax.combine)){
		cat("Combining taxon columns...                                           \r")
		flush.console()
		if("clgen"%in% tax.combine) dat$clgen <- paste0(dat$class, dat$genus, sep="_")
		cat("Combining taxon columns...OK                                           \n")
		flush.console()
	}


	# Core subset definition. i.e. include
	if(!is.null(include)){
		cat("Defining included subsets...                                           \r")
		flush.console()

		# what to include
		globalInclude <-rep(FALSE, nrow(dat))
		
		for(i in 1:length(include)){ 

			# whatever this is...
			cur <- include[[i]]

			# use a pre-defined dataset
			if(inherits(cur, "character")){
				# define environment to load these..
				e <- environment()

				# load data 
				data(list=cur, package="fossilx", envir=e)

				# execute filter
				incIndex <- fossilx:::column_filter(dat, e[[cur]])

				# create a union
				globalInclude <- globalInclude | incIndex
			}
		}

		dat<- dat[globalInclude,]
		cat("Defining included subsets...OK                                           \n")
		flush.console()
	 }

	# omission
	if(!is.null(omit)){
		cat("Defining omitted subsets...                                           \r")
		flush.console()
		for(i in 1:length(omit)){
			# whatever this is...
			cur <- omit[[i]]
			# use a pre-defined dataset
			if(inherits(cur, "character")){
				# colon-syntax
				if(grepl(":", cur)){
					# the key-value pair
					omission <-  colonfilter_to_dataframe(cur)
					# and do the omissions
					incInd <- fossilx:::column_filter(dat, omission, include=FALSE)
					dat<- dat[incInd,]
				}else{# built in dataset
					# define environment to load these..
					e <- environment()

					# load data 
					data(list=cur, package="fossilx", envir=e)

					# execute filter
					incIndex <- fossilx:::column_filter(dat, e[[cur]], include=FALSE)
					dat<- dat[incIndex,]

				}
			}
		}
		cat("Defining omitted subsets...OK                                           \n")
		flush.console()
	}


	# environmental categories
	if(!is.null(env.categories)){
		cat("Categorizing environmental data...                                           \r")
		flush.console()
		if(env.categories=="divDyn"){

			# create a new env for storing the results 
			e <- environment()
			data(keys, package="divDyn", envir=e)

			dat$lith<-divDyn::categorize(dat$lithology1,e$keys$lith)

			# batyhmetry
			dat$bath <- divDyn::categorize(dat$environment,e$keys$bath) 

			# grain size
			dat$grain <- divDyn::categorize(dat$lithology1,e$keys$grain) 

			# reef or not?
			dat$reef <- divDyn::categorize(dat$environment, e$keys$reef) 
			dat$reef[dat$lith=="clastic" & dat$environment=="marine indet."] <- "non-reef" # reef or not?/2

			# onshore - offshore
			dat$depenv <- divDyn::categorize(dat$environment,e$keys$depenv) 
		}
		cat("Categorizing environmental data...OK                                           \n")
		flush.console()
	}

	# stratigraphic assignments
	if(!is.null(strat)){
		# stg - binnning
		stg <- grep("stg",strat)
		if(length(stg)>0){
			cat("Assigning stratigraphic bins - 'stg'...                                           \r")
			flush.console()
			# grab version
			ver <- unlist(lapply(strsplit(strat[stg],"_"), function(x) x[2]))

			# do the assignments
			stgBin <- pbdb_strat_stg(dat, version=ver)

			# add to the data.frame 
			dat <- cbind(dat, stg=stgBin)

		}
		# ten - binnning
		ten <- grep("ten",strat)
		if(length(ten)>0){
			cat("Assigning stratigraphic bins - 'ten'...                                           \r")
			flush.console()
			# grab version
			ver <- unlist(lapply(strsplit(strat[ten],"_"), function(x) x[2]))

			# do the assignments
			tenBin <- pbdb_strat_ten(dat, version=ver)

			# add to the data.frame 
			dat <- cbind(dat, ten=tenBin)
			
		}

		# stb - binning
		stb <- grep("stb",strat)
		if(length(stb)>0){
			cat("Assigning stratigraphic bins - 'stb'...                                           \r")
			flush.console()
			# generate a timescale object
			stb_stages <- pbdb_timebins(dat, stagecolumn="time_contain", min_ma="min_ma", max_ma="max_ma", bin="stb")
			# resolve based on timescale
			stbVec <- pbdb_strat_stb(dat, ts=stb_stages, bin="stb", stagecolumn="time_contain")
			# add to the data.frame 
			dat <- cbind(dat, stb=stbVec)
		}
		cat("Assigning stratigraphic bins...OK                                           \n")
		flush.console()
		
	}

	return(dat)

}
