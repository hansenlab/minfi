\name{readGEORawFile}
\alias{readGEORawFile}
\title{
  Read in Unmethylated and Methylated signals from a GEO raw file.
}
\description{
  Read in Unmethylated and Methylated signals from a GEO raw file.
}
\usage{
readGEORawFile(filename, sep = ",", Uname = "Unmethylated signal", Mname = "Methylated signal", row.names = 1, pData = NULL, array = "IlluminaHumanMethylation450k", annotation = .default.450k.annotation, mergeManifest = FALSE, showProgress = TRUE)
}

\arguments{
  \item{filename}{
    The name of the file to be read from.
  }
  \item{sep}{
    The field separator character.  Values on each line of the
    file are separated by this character.  
  }
  \item{Uname}{
    A string that uniquely identifies the columns containing the
    unmethylated signals. 
  }
  \item{Mname}{
    A string that uniquely identifies the columns containing the
    methylated signals.}
  \item{row.names}{
    The column containing the feature (CpG) IDs. 
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
 450K experiments uploaded to GEO typically include a raw data file as
 part of the supplementary materials. Unfortunately there does not
 appear to be a standard format. This function provides enough
 flexibility to read these files. Note that you will likely need to
 change the \code{sep}, \code{Uname}, and \code{Mname} arguments and
 make sure the first column includes the feature (CpG) IDs. You can use the
 \code{\link{readLines}} function to decipher how to set these
 arguments.
 
 Note that the function uses the \code{\link[data.table]{fread}}
function in the \pkg{data.table} package to read the data. To install
\pkg{data.table} type \code{install.packages("data.table")}. We use this
  package because the files too large for \code{read.table}.
}
\value{
  A \linkS4class{GenomicMethylSet} object.
  
}
\author{
  Rafael A. Irizarry\email{rafa@jimmy.harvard.edu}.
}
\seealso{
  \code{\link{getGenomicRatioSetFromGEO}}
}
\examples{
\dontrun{
    ##Here is an example for GSE29290
    gmset=read.GEORawFile("GSE29290_Matrix_Signal.txt",Uname="Signal_A",Mname="Signal_B",sep="\t")
    ##Note: I changed 0,00 to 0.00 in the file first using sed:
    ##sed s/,/\\./g GSE29290_Matrix_Signal.txt >  GSE29290_Matrix_Signal.txt
}
}
