

library("Rglpk")
Objvec <- c(2, 3)
varnames <- c("Lane_Paint","Sign_Paint")

RHSvec <- c(6, 8, -1, 2)
constraintnames <- c("Paint A limit","Paint B limit","L can't exceed S by >1","L less than 2 tons")

Amat <- rbind( c(1,2),
               c(2,1),
               c(1,-1),
               c(1,0))


## solve using glpk via the API
library(glpkAPI)
## initialize model
lp<- initProbGLPK()
## model dimensions
nrows <- dim(Amat)[1]
ncols <- dim(Amat)[2]

## use sparse format for model data
nnonzero <- 7
rindexnonzero <- c(1, 1, 2, 2, 3, 3, 4)
cindexnonzero <- c(1, 2, 1, 2, 1, 2, 1)
valuesnonzero <- c(1, 2, 2, 1, 1, -1, 1)

# row upper and lower bounds
rlower <- rep(0,nrows)
rlower[3] <- -1000
rupper <- RHSvec

# column upper and lower bounds
clower <- rep(0,ncols)
cupper <- rep(1000,ncols)

# maximize objective GLP_Max (minimize with GLP_MIN)
setObjDirGLPK(lp,GLP_MAX)

# tell model how many rows and columns
addRowsGLPK(lp, nrows)
addColsGLPK(lp, ncols)

# add column limits
setColsBndsGLPK(lp,c(1:ncols), clower, cupper)
setRowsBndsGLPK(lp,c(1:nrows),rlower,rupper)
setObjCoefsGLPK(lp,c(1:ncols),Objvec)

# load constraint matrix coefficients
loadMatrixGLPK(lp,nnonzero, rindexnonzero, cindexnonzero,valuesnonzero)

# solve LP problem using Simplex Method
solveSimplexGLPK(lp)

# get results of solution
# solve status 5 = optimal solution found
getSolStatGLPK(lp)
status_codeGLPK(getSolStatGLPK(lp))

# objective function value
getObjValGLPK(lp)
# value of variables in optimal solution
getColsPrimGLPK(lp)
# status of each variable in optimal solution 1 = basic variable
getColsStatGLPK(lp)

# get dual values for each row/constraint
getRowDualGLPK(lp,1)
getRowDualGLPK(lp,2)
getRowDualGLPK(lp,3)
# this is supposed to get all at once, but doesn't seem to work
#getRowsDualGLPK(lp)

printRangesGLPK
##printRangesGLPK(lp, numrc = 0, rowcol = NULL, fname = "test.txt")