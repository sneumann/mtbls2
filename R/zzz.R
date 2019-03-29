# save compressed to avoid warnings:
# save(mtbls2Set, file="data/mtbls2.rda", compress='xz') 
load(file.path(find.package("mtbls2"), "data", "mtbls2.rda"))

# The following path fixup requires that the mtbls2.rda in BioC has the relative path
#

.onAttach <- function(libname, pkgname) {
    attr(mtbls2Set, "filepaths") <- file.path(find.package("mtbls2"), gsub("/", .Platform$file.sep, attr(mtbls2Set, "filepaths")))    
}
