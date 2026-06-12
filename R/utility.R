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

# General utily function for creating a binning as descrbied in the divDyn SOM (ddPhanero)
#
# @param x AN occurrence data frame.
# @param categs Category list, as the elements of divDyn keys$stgInt 
# @param out column name of the resulting stratigraphiy assignment
# @param early_interval column name of the early interval
# @param late_interval column name of the late interval
pbdb_assign_bin <- function(x, categs, out, early_interval="early_interval", late_interval="late_interval"){

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
