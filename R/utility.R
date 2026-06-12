# Utility function to filter occurrences
# @param x An occurrence data.frame
# @param y A data frame with two columns: column and entries. The first column specifies which columns will be searched and the second for what entries.
# @param include If TRUE, then the entries specied in y are retained only. If FALSE, these are omitted.
# @return A logical vector with row indices to keep.
column_filter <- function(x, y, include=TRUE){
	# the ranks to filter with
	theColumns <- levels(factor(y$column))

	# which are to be included
	index <- rep(FALSE, nrow(x))

	for(i in 1:length(theColumns)){
		currCol <- theColumns[i]
		# the current subset to compare with
		current <- y[y$column==currCol,"entries"]
		if(include){
			index[x[,currCol]%in%current] <- TRUE
		}else{
			index[!x[,currCol]%in%current] <- TRUE
		}
	}

	return(index)

}

# replace entries in a column with a named vector
replace_entries <- function(x, where="stg", with=strat.camb_20180831,by="collection_no"){
# the collections entered
	colls <- as.character(dat[,by])
	
	# which are cambrian?
	bool <- colls%in%names(with)

	# the cambrian collections - as they are in the occurrence dataset
	subColls <- colls[bool]
	
	# order/assign the stg accordingly
	subStg <- with[subColls]
	
	# copy original
	newStg <- x[,where] 
	
	# replace the missing entries
	newStg[bool]  <- subStg
	
	#
	x[,where] <- newStg

	return(x)
}

