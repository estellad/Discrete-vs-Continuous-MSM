infjack.glm <- function(glm.obj,groups)
{
##############################################
# Original author:  T. Lumley                #
#                   Biostatistics            #
#                   University of Washington #
##############################################
#
## 
# Run estfun.glm to get GLM score function
##       
	umat <- estfun.glm(glm.obj)
##
# Sum scores within each cluster (individual)
##        
	usum <- rowsum(umat,groups,reorder=F)
##
# Calculate the empirical variance of the summed scores
# and then compute sandwich variance-covariance matrix
##        
	modelv <- summary(glm.obj)$cov.unscaled        
	output <- modelv%*%(t(usum)%*%usum)%*%modelv
##
# Output empirical variance-covariance matrix 
##
	output
}
