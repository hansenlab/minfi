\name{minfinews}
\title{Minfi News}
\encoding{UTF-8}

\section{Version 1.3}{
  \itemize{
    \item Updated preprocessSwan to fix a bug when mSet was not set to
    the default value of NULL.  Specifically, now the "counts" tables is
    used to construct "subset".
    \item Changed the function manifestNew() to
    IlluminaMethylationManifest().
    \item Added IlluminaMethylationAnnotation().
    \item Added placeholders for unit testing based on RUnit.
    \item Introduced a new show method for MethylSet and RGChannelSet,
    derived from the eSet method in Biobase.
    \item The annotation slot of a MethylSet/RGChannelSet is now
    intended to _not_ be a scalar, but instead have length 2 with
    components 'array' and 'annotation'.  This foreshadows introdution of
    annotation packages for use with minfi.
    \item Reorganization of R files; rewriting of the man pages for
    MethylSet, RGChannelSet.
    \item getMeth, getUnmeth, getBeta, getM are now methods.
  }
}

\section{Version 1.1}{
  \itemize{
    \item Changed NAMESPACE file
    \item Defined constructors for MethylSet, RGChannelSet,
    RGChannelSetExtended.
    \item Included a version number in the class definition for
    MethylSet and RGChannelSet.  Old objects can be updated by calls of
    the form updateObject(Mset).
    \item read.manifest (not exported) updated to include nCpGs.
    \item preprocessSwan was added.  Still work in progress.
    \item Changed background calculation in preprocessSwan.
    \item Added a section to the vignette describing preprocessSwan.
    \item Bug fix: ilogit2 is now in base (it used to be base e).  Thanks
    to Time Triche, Jr <tim.triche@gmail.com>.
    \item Added and dcoumented the IlluminaMethylationAnnotation class;
    still work in progess.
    \item Moved package vignette from inst/doc to vignettes.
  }
}

\section{Version 0.99}{
  \itemize{
    \item Initial release to Bioconductor.
    \item Added NEWS file.
    \item Bugfix to vignette.
    \item readIDAT is now exported by crlmm, implying that we can import
  this function through NAMESPACE.
  }
}
