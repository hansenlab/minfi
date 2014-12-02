\name{readTCGA}
\alias{readTCGA}
\title{
Read in tab deliminited file in the TCGA format
}
\description{
Read in tab deliminited file in the TCGA format
}
\usage{
readTCGA(filename, sep = "\t", keyName = "Composite Element REF", Betaname = "Beta_value", pData = NULL, array = "IlluminaHumanMethylation450k", annotation = .default.450k.annotation, mergeManifest = FALSE, showProgress = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
The name of the file to be read from.
}
  \item{sep}{
 The field separator character.  Values on each line of the
          file are separated by this character.  
}
  \item{keyName}{
The column name of the field containing the feature IDs.
}
  \item{Betaname}{
The character string contained all column names of the beta value fields.
}
  \item{pData}{
A \code{DataFrame} or \code{data.frame} describing the samples represented by the columns of
\code{mat}. If the rownames of the \code{pData} don't match the colnames
of \code{mat} these colnames will be changed. If \code{pData} is not
supplied, a minimal \code{DataFrame} is created.
}
  \item{array}{
Array name.
}
  \item{annotation}{
The feature annotation to be used. This includes the location of
features thus depends on genome build.
}
  \item{mergeManifest}{
Should the Manifest be merged to the final object.
}
  \item{showProgress}{
    TRUE displays progress on the console. It is
    produced in fread's C code.
  }
}
\details{
  This function is a wrapper for
  \code{\link{makeGenomicRatioSetFromMatrix}}. It assumes a very specific
  format, used by TCGA, and then uses the \code{\link[data.table]{fread}}
  function in the \pkg{data.table} package to read the data. To install
  \pkg{data.table} type \code{install.packages("data.table")}. We use this
  package because the files too large for \code{read.table}.
    
  Currently, an example of a file that this function reads is here:
  \url{http://gdac.broadinstitute.org/runs/stddata__2014_10_17/data/UCEC/20141017/gdac.broadinstitute.org_UCEC.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2014101700.0.0.tar.gz}. Note it is a 8.1 GB archive.
}
\value{
  A \linkS4class{GenomicRatioSet} object.
}
\author{
  Rafael A. Irizarry\email{rafa@jimmy.harvard.edu}.
}
\seealso{
  \code{\link{makeGenomicRatioSetFromMatrix}}
}
\examples{
\dontrun{
    filename <- "example.txt" ##file must be in the specicif TCGA format
    readTCGA(filename)
}
}
