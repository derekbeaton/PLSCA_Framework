rplsca.supProj <- function(X,Y,lambda.x=NULL,lambda.y=NULL,projmat.X=NULL,projmat.Y=NULL,fi.scores,fj.scores,Dv,X.RW,Y.RW,p.scores=NULL,q.scores=NULL){

	if( is.null(lambda.x) & is.null(lambda.y) ){
		lambda.x <- lambda.y <- 0
	}else if( !is.null(lambda.x) & is.null(lambda.y) ){
		lambda.y <- lambda.x
	}else if( is.null(lambda.x) & !is.null(lambda.y) ){
		lambda.x <- lambda.y
	}

	X.preproc <- takanePreProc(X,lambda.x,projmat.X)
	Y.preproc <- takanePreProc(Y,lambda.y,projmat.Y)
	X.dat <- X.preproc$dat
	Y.dat <- Y.preproc$dat
		rm(X.preproc)
		rm(Y.preproc)

	R.dat <- t( X.dat ) %*% ( Y.dat )
	
	##actually, the weights should be the original weights...
	sup.fj <- as.matrix( matrix(1/Y.RW,ncol(R.dat),ncol(fi.scores),byrow=F) * (t(R.dat) %*% (fi.scores %*% diag(Dv^-1))) )
	sup.fi <- as.matrix( matrix(1/X.RW,nrow(R.dat),ncol(fj.scores),byrow=F) * (R.dat %*% (fj.scores %*% diag(Dv^-1))) )


	if( !is.null(p.scores) & !is.null(p.scores)){
		sup.LX <- X.dat %*% p.scores
		sup.LY <- Y.dat %*% q.scores
		return(list(sup.fi=sup.fi,sup.fj=sup.fj,sup.LX=sup.LX,sup.LY=sup.LY,R= R.dat))
	}else{
		return(list(sup.fi=sup.fi,sup.fj=sup.fj,R= R.dat))
	}
}
