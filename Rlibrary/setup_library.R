############################################################
####
####    SET-UP THE FILES AND FOLDERS TO CREATE R LIBRARY
####
############################################################

library(Rcpp)

# Input the name of the package here
mypackage <- "stigene"

# path to the model code:
code.model.folder <- "../code_model/"

# Retrieve all relevant C++ files:
cppfiles <- system(paste0("ls ",code.model.folder,"*.cpp"),intern = TRUE)
cppfiles <- cppfiles[!grepl(pattern = "main",x = cppfiles)]   # <-- Remove "main" files
hfiles <- system(paste0("ls ",code.model.folder,"*.h"),intern = TRUE)
c.path <- c(cppfiles, hfiles)

# This generates all the necessary files 
# when creating an R package from scrath 
# with the view of interfacing with C++ (using Rcpp)
#?Rcpp::Rcpp.package.skeleton
Rcpp.package.skeleton(name =  mypackage,
                      example_code = FALSE,
                      author = "David Champredon",
                      cpp_files = c.path, 
                      force = TRUE
)
