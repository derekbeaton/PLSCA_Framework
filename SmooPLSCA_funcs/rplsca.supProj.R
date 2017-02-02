rplsca.supProj <- function(X,Y,lambda.x=NULL,lambda.y=NULL,projmat.X=NULL,projmat.Y=NULL,fi.scores,fj.scores,Dv,p.scores=NULL,q.scores=NULL){

	if( is.null(lambda.x) & is.null(lambda.y) ){
		lambda.x <- lambda.y <- 0
	}else if( !is.null(lambda.x) & is.null(lambda.y) ){
		lambda.y <- lambda.x
	}else if( is.null(lambda.x) & !is.null(lambda.y) ){
		lambda.x <- lambda.y
	}

	X.preproc <- takanePreProc(X,lambda.x,projmat.X)
	Y.preproc <- takanePreProc(Y,lambda.y,projmat.Y)
	X.dat <- sqrt.mat(X.preproc$LW) %*% X.preproc$dat
	Y.dat <- sqrt.mat(Y.preproc$LW) %*% Y.preproc$dat

	R.dat <- t( X.dat ) %*% ( Y.dat )
	sup.fj <- svd.pinv_sq(Y.preproc$RW) %*% t(R.dat) %*% (fi.scores %*% diag(Dv^-1))
	sup.fi <- svd.pinv_sq(X.preproc$RW) %*% R.dat %*% (fj.scores %*% diag(Dv^-1))
	if( !is.null(p.scores) & !is.null(p.scores)){
		sup.LX <- X.dat %*% p.scores
		sup.LY <- Y.dat %*% q.scores
		return(list(sup.fi=sup.fi,sup.fj=sup.fj,sup.LX=sup.LX,sup.LY=sup.LY))
	}else{
		return(list(sup.fi=sup.fi,sup.fj=sup.fj))
	}
}
