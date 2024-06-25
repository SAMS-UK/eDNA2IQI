#the only purpose of this project is
#to compile the package
#drag editing/updated files into the ...push r folder
#from ...debugging project.
#devtools will compile from there
#assuming header information given
#use roxygen2 for header templates
#once the package has been re-compiled, check the version by running in project 'edna2iqi_apply'.
#may be necessary to delete previous versions under the 'man' folder here
#devtools then rebuilds them.

#remove previously installed version.
unlink(paste(.libPaths(),"/eDNA2IQI",sep=""),recursive=TRUE)
#update.packages(ask = FALSE)

#install.packages("pacman")         # Install pacman package
#library("pacman")                  # Load pacman package
#p_update(update=FALSE)
#p_update()


library(devtools)
devtools::check_man()
devtools::check()
devtools::document()
#devtools::test()
#?devtools::test
install()#press 3 for ignore all library updates.
library(eDNA2IQI)
