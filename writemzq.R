library(faahKO)

require(XML) || stop("We need library(XML) to write mzData")

write.mzq <- function(object, filename) {
    mzq <- buildMzq(object)

    ## the sink() workaround seems to be needed for proper indenting.
    sink(file=filename)
    cat(saveXML(mzq)) ##, prefix = '<?xml version="1.0"?>\n'))
    sink(NULL)
}

verify.mzq <- function(xmlfilename, xsdfilename) {
    xsd = xmlTreeParse(xsdfilename, isSchema =TRUE, useInternal = TRUE)
    doc = xmlInternalTreeParse(xmlfilename)
    xmlSchemaValidate(xsd, doc)    
}

buildCVlist <- function() {

    CvList <- newXMLNode("CvList")

    CVs <- list(
        newXMLNode("Cv", attrs=c(
                             id="PSI-MS",
                             fullName="Proteomics Standards Initiative Mass Spectrometry Ontology",
                             version="2.29.0",
                             uri="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo")),
        
        
        newXMLNode("Cv", attrs=c(
                             id="UO",
                             fullName="Unit Ontology",
                             version="1.20",
                             uri="http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo")))
    
    addChildren(CvList, kids=CVs)   
}

buildAuditCollection <- function() {
    newXMLNode("AuditCollection",
               .children=list(
                   newXMLNode("Person",
                              attrs=c(
                                  firstName="Steffen",
                                  lastName="Neumann",
                                  id="PERS_STN"),
                              .children=list(newXMLNode("Affiliation",
                                  attrs=c(organization_ref="ORG_IPB")))),
                   newXMLNode("Organization",
                              attrs=c(
                                  name="Leibniz Institute of Plant Biochemistry",
                                  id="ORG_IPB"))
                   ))
}

buildCvParams <- function(cvparams) {
    lapply(cvparams, function(cv) {
        newXMLNode("cvParam", attrs=cv)
    })               
}

buildAnalysisSummary <- function() {
    newXMLNode("AnalysisSummary",
               .children=buildCvParams(
                   list(c(accession="MS:1001834", cvRef="PSI-MS", name="LC-MS label-free quantitation analysis"),
                        c(accession="MS:1002019", cvRef="PSI-MS", name="label-free raw feature quantitation",  value="false"))
                   )) 
}

buildInputFiles <- function(xs) {
    f <- filepaths(xs)
    rfgs <- lapply(1:length(f), function(i) {
        x <- f[i]
        name <- basename(x)
        newXMLNode("RawFilesGroup",
                   attrs=c(id=paste("rfg_", i, sep="")),
                   .children=newXMLNode("RawFile",
                       attrs=c(location=x, name=name, id=paste("rf_", i, sep=""))))
    })
    newXMLNode("InputFiles", .children=rfgs)    
}

buildSoftwareList <- function() {
    newXMLNode("SoftwareList",
               .children=list(
                   newXMLNode("Software",
                              attrs=c(
                                  id="xcms",
                                  version=packageVersion("xcms")),
                              .children=buildCvParams(
                                  list(c(accession="MS:1001830", cvRef="PSI-MS", name="XCMS")))
                              )))
}

buildDataProcessingList <- function() {
    newXMLNode("DataProcessingList",
               .children=newXMLNode("DataProcessing",
                   attrs=c(order="1",
                       software_ref="xcms", id="DP1"),
                   .children=newXMLNode("ProcessingMethod",
                   attrs=c(order="1"))))
}

buildAssayList <- function(xs) {
    f <- filepaths(xs)
    assays <- lapply(1:length(f), function(i) {
        x <- f[i]
        name <- basename(x)
        newXMLNode("Assay",
                   attrs=c(
                       rawFilesGroup_ref=paste("rfg_", i, sep=""),
                       name=x,
                       id=paste("assay_", i, sep="")),
                   .children=list(newXMLNode("Label",
                       .children=newXMLNode("Modification",
                           buildCvParams(list(c(accession="MS:1002038", cvRef="PSI-MS", name="unlabeled sample")))
                           ))))
    })
    newXMLNode("AssayList",
               attrs=c(id="AssayList_1"),
               .children=assays)        
}



buildStudyVariableList <- function(xs) {
    p <- phenoData(xs)
    sc <- sampclass(xs)

    svlevels <- levels(sampclass(faahko))
    StudyVariables <- lapply (seq(1,length(svlevels)), function(i) {
        sv <- svlevels[i]
        ##        assays <- paste (sampnames(xs)[sampclass(xs) == sv], collapse="-")
        assay_ref <- paste("assay", which(sampclass(xs) == sv), sep="_", collapse=" ")
        
        newXMLNode("StudyVariable",
                   attrs=c(
                       name=svlevels[i],
                       id=paste("SV_COLLAPSED", i, sep="_")),
                   .children=list(newXMLNode("Assay_refs", assay_ref)))
    })
    
    newXMLNode("StudyVariableList",
               .children=StudyVariables)
}


buildMzq <- function(xs) {
    mzqVersion="1.0.0"
    schemaLocation="http://psidev.info/psi/pi/mzQuantML/1.0.0 ../../../schema/mzQuantML_1_0_0.xsd"
    
    mzq = xmlTree(tag="MzQuantML",
        attrs=c(
            version=mzqVersion,
            id="mzq-generated-from-xcms",
            creationDate="2012-11-26T15:13:08.330Z",
            "xsi:schemaLocation"=schemaLocation),
        namespaces = c(
            "http://psidev.info/psi/pi/mzQuantML/1.0.0",
            xsi="http://www.w3.org/2001/XMLSchema-instance")    
        )    
    
    mzq$addNode(buildCVlist())
    mzq$addNode(buildAuditCollection())
    mzq$addNode(buildAnalysisSummary())
    mzq$addNode(buildInputFiles(xs))
    mzq$addNode(buildSoftwareList())
    mzq$addNode(buildDataProcessingList())
    mzq$addNode(buildAssayList(xs))
    mzq$addNode(buildStudyVariableList(xs))
##    mzq$addNode(buildSmallMoleculeList(xs))
##    mzq$addNode(buildFeatureList(xs))    
    mzq$closeTag()                                        
}

write.mzq(faahko, "writemzq.mzq.xml")

verify.mzq(xmlfilename="writemzq.mzq.xml",
           xsdfilename="mzQuantML_1_0_0.xsd")# $errors[[1]]$msg



