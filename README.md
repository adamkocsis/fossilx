# fossilx

Extension and Post-processing of Fossil Occurrence Records


## Goal 

The purpose of this package is to define standard (non-analytical) post-processing steps of data downloaded from online fossil occurrence occurrence databases. Currently the project focuses on the Paleobiolgy Database but its is aimed at having a broader scope. The post-processing is documented declaratively and optimized for maximum portability.

## Installation

The project is in development, you can install it with:

``` R
library(devtools)
install_github("adamkocsis/fossilx")
```

## Example use

The following snippet recreates the data procesing steps of the [divDyn SOM](https://github.com/divDyn/ddPhanero/). 

``` R
library(fossilx)
library(chronosphere)

# download whole pbdb from Zenodo
dat <- chronosphere::fetch("pbdb", ver="20260412")

# example use
refined <- pbdb_extend(dat,
	tax.level="genus",
	tax.combine="clgen",
	include=list("tax.marine_1.0"),
	env.categories="divDyn",
	strat=c("stg_1.0", "ten_1.0", "stb"),
	omit=list("env.nonmarine_1.0", "lithification1:unlithified")
)
```

Syntax guidelines:
- Colons define key-value pairs: `lithification1:unlithified`, in this case records where the lithification1 column is unlithified will be omitted.
- Underscores define versions, in stratigraphic assignments these are version of procedures. In filtering steps these refer to specific items, cf. `data(env.nonmarine_1.0)`.

The framework `stb` represents a solution for the discrete binning based on the `time_contain` output of the pbdb. You can construct a relevant timescale with:

``` R
timebins <- pbdb_timebins(dat, stagecolumn="time_contain", min_ma="min_ma", max_ma="max_ma", bin="stb")
```

## To do...

- [x] Taxonomic and environmental filtering 
- [x] Environmental categorization
- [x] Stratigraphic assignments 
- [ ] GPM reconstructions with `rgplates`
- [ ] Adding body size data from Heim et al. 2015.
- [ ] Select data with abundance records 
- [ ] Species-level taxonomy
- [ ] Cleaning the data 
- [ ] Include climate proxy data
- [ ] Categorizing functional traits
- [ ] Recent list of genera - Taxon list from the PBDB 

