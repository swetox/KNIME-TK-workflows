########################################
######### Package installation #########
#######################################

# Run script in local, fully updated R
# to ensure correction installations 
# for running KNIME workflows

# Installation of ggplot:
install.packages("ggplot2") 	
packageVersion("ggplot2")		# version 2.2.1 used during study

# Installation of reshape2:
install.packages("reshape2")	
packageVersion("reshape2")		# version 1.4.2 used during study

# Installation of Rserve
install.packages("Rserve")
packageVersion("Rserve")      # version 1.7.3 used during study

# Installation of deSolve
install.packages("deSolve")
packageVersion("deSolve")     # version 1.20 used during study

# Installation of MESS
install.packages("MESS")
packageVersion("MESS")        # version 0.4.15 used during study

# Installation of dplyr
install.packages("dplyr")
packageVersion("dplyr")       # version 0.7.2 used during study

# Installation of grid
install.packages("grid")
packageVersion("grid")        # version 3.4.1 used during study

# Installation of fitdistrplus
install.packges("fitdistrplus")
packageVersion("fitdistrplus") # version 1.0.9 used during study

# Installation of tmvtnorm
install.packages("tmvtnorm")
packageVersion("tmvtnorm")    # version 1.4.10 used during study

########## Location of R ############
R.home()			
R.version.string 				# version 3.4.1 used for during study
