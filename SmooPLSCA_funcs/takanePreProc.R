###Is this the same as the projmaker?
takanePreProc <- function(t.dat.nom,lambda,projmat=NULL){
	####I can update this later by expanding what L can be.
		###Which, technically, can be a diagonal matrix. But requires more effort for LW and RW
	
	t.dat.center <- scale(t.dat.nom,scale=F)
	if( is.null(projmat) ){
		projmat <- t(t.dat.center) %*% svd.pinv_sq( tcrossprod(t.dat.center) ) %*% t.dat.center
	} ##I can just make this a diagonal... probably.
	LW <- diag(nrow(t.dat.nom)) + (lambda * svd.pinv_sq( tcrossprod(t.dat.center) ) )	
	rownames(LW) <- colnames(LW) <- rownames(t.dat.nom)
	RW <- diag(colSums(t.dat.nom)) + (lambda * projmat)
	rownames(RW) <- colnames(RW) <- colnames(t.dat.nom)
	
	return(list(dat=t.dat.center,LW=LW,RW=RW,projmat=projmat))
}