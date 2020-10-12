
estfun.glm <- function(glm.obj)
{
##############################################
# Original author:  T. Lumley                #
#                   Biostatistics            #
#                   University of Washington #
##############################################
#
##
# Create X matrix from glm object 
##
	if(is.matrix(glm.obj$x))
		xmat <- glm.obj$x
	else {
		mf <- model.frame(glm.obj)
		xmat <- model.matrix(terms(glm.obj), mf)
	} 
##
# Calculate variance weight and residual (Y-fitted)
##
	output <- residuals(glm.obj, "working") * glm.obj$weights * xmat
##
# Output this matrix
##
	output
}
