library(CAMERA)

##
## Load processed data
##
load("/home/sneumann/tex/papers/2009tandemms/R/acquisitionCYP2011.Rdata")

##
## Load manual annotation, remove CAMERA-only annotation
##
#esm3.orig <- read.csv("11306_2012_401_MOESM3_ESM.csv", skip=5, stringsAsFactors=FALSE)
esm3.orig <- read.csv("Supplemental_S6_compound_annotation.csv", skip=5, stringsAsFactors=FALSE)
esm3 <- esm3.orig[sapply(esm3.orig[,"ion.type.1"], function(x) nchar(x)>0, USE.NAMES=FALSE), ]

##
## load ISA assay files
##
a.isa <- read.delim("a_mtbl2_metabolite profiling_mass spectrometry.txt", stringsAsFactors=FALSE)
a.samples <- a.isa$Sample.Name

##
## load existing maf files
##
maf.orig <- read.delim("a_mtbl2_metabolite profiling_mass spectrometry_maf.csv")

##
## These columns are defined by mzTab
##

maf.std.colnames <- c("identifier", "chemical_formula", "description",
"mass_to_charge", "fragmentation", "charge", "retention_time",
"taxid", "species", "database", "database_version", "reliability",
"uri", "search_engine", "search_engine_score", "modifications",
"smallmolecule_abundance_sub", "smallmolecule_abundance_stdev_sub",
"smallmolecule_abundance_std_error_sub")

all.colnames <- c(maf.std.colnames, a.samples)

##
## Get preprocessed data
##

pl <- getPeaklist(an)
charge <- sapply(an@isotopes, function(x) {ifelse( length(x) > 0, x$charge, NA) })
l <- nrow(pl)
abundance <- groupval(an@xcmsSet, value="into")

##
## Now assemble new maf
##

maf <- data.frame(identifier = character(l), 
                  chemical_formula = character(l), 
                  description = character(l), 
                  mass_to_charge = pl$mz, 
                  fragmentation = character(l), 
                  charge = charge, 
                  retention_time = pl$rt, 
                  taxid = character(l), 
                  species = character(l), 
                  database = character(l), 
                  database_version = character(l), 
                  reliability = character(l), 
                  uri = character(l), 
                  search_engine = character(l), 
                  search_engine_score = numeric(l), 
                  modifications = character(l), 
                  smallmolecule_abundance_sub = numeric(l), 
                  smallmolecule_abundance_stdev_sub = numeric(l), 
                  smallmolecule_abundance_std_error_sub = numeric(l),
                  abundance, stringsAsFactors=FALSE)

##
## Now overlay information from supplemental
## Map manual annotation into XCMS peaklist
##

## pre-calculate fixed-digit character representation
## for comparison with CSV

esm3.char <- cbind(m.z = formatC(esm3[,"m.z"], format="f", digits=4, flag = "#"),
                  ret..time..s. = formatC(esm3[,"ret..time..s."], format="f", digits=1, flag = "#"))

pl.char <- cbind(mz = formatC(pl[,"mz"], format="f", digits=4, flag = "#"),
                 rt = formatC(pl[,"rt"], format="f", digits=1, flag = "#"))

w <- apply(esm3.char, 1, function(x) {
  which (x["m.z"] == pl.char[,"mz"]
         & x["ret..time..s."] == pl.char[,"rt"])
})

## Christoph's manual annotation:
maf[w, "modifications"] <- esm3[,"ion.type.1"]
maf[w, "chemical_formula"] <- esm3[,"elemental.composition"]
maf[w, "description"] <- esm3[,"compound.name"]

## We have PubChem for a few of them
maf[w, "database"] <- "PubChem"
maf[w, "identifier"] <- esm3[,"PubChem.CID"]

maf_character <- apply(maf, 2, as.character)

write.table(maf_character, file="a_mtbl2_metabolite profiling_mass spectrometry_maf2.csv",
            row.names=FALSE, col.names=all.colnames, quote=TRUE, sep="\t", na="\"\"")
