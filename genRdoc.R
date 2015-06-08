## To compile the package use RStudio. Note that roxygen2 is NOT used for 
## generating documentation insted R.oo's internal Rdoc class is used. 
## 1) Build the package (no files in the man folder yet) 
## 2) Run the code below to generate the Rd files
## 3) Build the package (now with files in the man folder) 
## 4) Run the code below to generate the Rd files

# Generate Rd files. MUST be run with working dir equal the R source folder!
library(discretizeGaussian)
setwd("./discretizeGaussian/R")
doc<-Rdoc()
Rdoc$compile(filename=".*[.]R$", verbose=F, source=F, check=TRUE)
setwd("../..")


# f <- list.files(".",pattern="\\.[R]",full.names=TRUE)
# a <- lapply(f,source)
# doc$compile(check=F, source=FALSE) #check=TRUE, debug=TRUE
# setwd("../..")

