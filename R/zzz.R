load(file.path(find.package("mtbls2"), "data", "mtbls2.rda"))

# This requires that the mtbls2.rda in BioC has the relative path
attr(mtbls2Set, "filepaths") <- file.path(find.package("mtbls2"), gsub("/", .Platform$file.sep, attr(mtbls2Set, "filepaths")))
