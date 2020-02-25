## First load the function
source("../downsampleFileList.R")
## Specify required parameters
sample.n <- c(100, 80, 60, 40, 30, 20)
file.list <- "HG00733.files.txt"
## Specify 'outputdirectory' if not the same as 'file.list' directory
## Run the function
downsampleFileList(file.list = file.list, file.ID.field = 1, sample.n = sample.n)


#' This function takes file names listed as a single column in a .txt, .csv or .tsv file
#' and reports a user defined random subset(s) of the original number of file names.
#'
#' @param file.list A path to a file that contains a list of files to be downsampled.
#' @param outputdirectory Uset defined directory to export downsampled file lists. (default: 'file.list' directory)
#' @param file.ID.field A field number to extract filenames from the file. (default: 1)
#' @param sample.n User defined number of subsets to create from files in 'file.list'.
#' @param sep A field separator used in 'file.list' file.
#' @return A value \code{NULL}.
#' @author David Porubsky
#'
downsampleFileList <- function(file.list, outputdirectory=NULL, file.ID.field=1, sample.n=NULL, sep=NULL) {
  ## Load file list
  if (!is.null(sep) & is.character(sep)) {
    file.names <- utils::read.table(file.list, sep = sep, stringsAsFactors = FALSE)
  } else {
    file.names <- utils::read.table(file.list, stringsAsFactors = FALSE)
  }
  file.names <- file.names[,file.ID.field]
  
  ## Create user defined outputdirectory or report into file.list directory
  if (!is.null(outputdirectory) & is.character(outputdirectory)) {
    if (!dir.exists(outputdirectory)) {
      dir.create(outputdirectory)
    }
  } else {
    outputdirectory <- dirname(file.list)
  }
  
  ## Check if user defined sample size are not larger then total number of files in file.list
  total.n <- length(file.names)
  if (any(sample.n > total.n)) {
    to.remove <- sample.n[sample.n > total.n]
    warning("Removing sample size ", to.remove, ", larger than the total ", total.n, " files in ", file.list)
    sample.n <- sample.n[sample.n <= total.n]
  }
  
  ## Get filename prefix
  file.prefix <- basename(file.list)
  file.prefix <- gsub(file.prefix, pattern = ".txt|.csv|.tsv", replacement = "")
  
  ## Report various sample size in a file
  for (n in sample.n) {
    sub.sample <- sample(file.names)[1:n]
    destination <- file.path(outputdirectory, paste0(file.prefix, '.sample', n, '.txt'))
    utils::write.table(sub.sample, file = destination, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  return(NULL)
}
