# Record the active pinned reference and BLAS identity before fitting.
args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(susieR))
stopifnot(as.character(packageVersion('susieR')) == '0.16.6')
dyn.load(args[[1L]])
probe <- .Call('prusie_blas_info')
stopifnot(!is.na(probe$threads), probe$threads == 1L)
write_json(list(R=R.version.string,version=as.character(packageVersion('susieR')),
  path=normalizePath(find.package('susieR')),
  blas=normalizePath(trimws(extSoftVersion()[['BLAS']])),blas_probe=probe,
  session=I(capture.output(sessionInfo()))),args[[2L]],auto_unbox=TRUE,pretty=TRUE)
