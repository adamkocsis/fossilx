#' Filter PBDB to a specific level of taxonomy
#' 
#' The function ensures that 
#' 
#' @param x A PBDB data.frame.
#' @param level A taxonomic range
#' @param omit Should the filter omit occurrences, or create a new logical column that indicate which should be kept?
#' @export
pbdb_taxon_quality_filter <- function(x, level="genus", omit=TRUE){
	# match the arguments
	level <- match.arg(level, "genus")

	if(level=="genus"){
		# filter records not identified at least to genus
		include <- x$accepted_rank %in% c("genus", "species", "subgenus", "subspecies")

		# omit non-informative genus entries
		include <- include & x$genus!=""
	}

	if(omit){
		x <- x[include,]
	}else{
		x <- cbind(x, include)
		colnames(x)[ncol(x)] <- paste0("filter_", level)
	}
	return(x)
}


