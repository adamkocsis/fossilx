# Utility function to filter occurrences
 
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
