##will need: https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/bdiag.html
	###OK so this is a just a block matrix where I can already know the values.
	###The values are based on # of variables.
#projmat.maker2 <- function(t.dat.center,colDesign,tol=.Machine$double.eps*100){

projmat.maker2 <- function(t.dat.center,tol=sqrt(.Machine$double.eps)){
	

	tJ <- t(t.dat.center) %*% ginv(tcrossprod(t.dat.center)) %*% t.dat.center	
	tJ[abs(tJ) < tol] <- 0
	rownames(tJ) <- colnames(tJ) <- colnames(t.dat.center)
	return(tJ)

}