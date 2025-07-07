#Functional regression - Canadian weather example
library(fda)

# Load data and extract components


# Response basis (Fourier)

phi <- function(t) eval.basis(t, basis_Y)
Phi <- phi(1:12)

# Coefficient matrix for functional response
Y <- temp
C <- t(solve(t(Phi) %*% Phi) %*% t(Phi) %*% Y)



# Basis for beta functions (B-splines)
K_beta <- 8
basis_beta <- create.bspline.basis(rangeval = range, nbasis = K_B)

# Penalty matrices and inner products
J_YY <- inprod(basis_Y, basis_Y)
J_BB <- inprod(basis_B, basis_B)
J_YB <- inprod(basis_Y, basis_B)
R <- eval.penalty(basis_B, int2Lfd(2))

# Regularization
lambda <- 0.1
Lambda <- diag(rep(lambda, ncol(Z)))

# Estimate B
vec_B <- solve(J_BB %x% (t(Z) %*% Z) + R %x% Lambda) %*% as.vector(t(Z) %*% C %*% J_YB)
B <- matrix(vec_B, nrow = ncol(Z), ncol = basis_B$nbasis, byrow = TRUE)
rownames(B) <- colnames(Z)
colnames(B) <- basis_B$names

# Functional representation of beta functions (optional, no plotting)
beta <- fd(coef = t(B), basisobj = basis_B,
           fdnames = list("t", colnames(Z), "beta"))

#We want to obtain predictions, we look at the chain of mappings

#Y2cMap, takes Y and returns the smoothing coefficients c--------
K_y <- 7 #nr of basis to smooth y (fourier)
basis_Y <- create.fourier.basis(rangeval = c(0,12), nbasis = K_y)

phi <- function(t) eval.basis(t, basis_Y) #this takes any t and gives the basis coordinates
Phi <- phi(1:12) #This contains the basis coordinates for all the points in T

Y <- t(temp) #Matrix N times n with the original observations for Y
#Each y(t) is a random vector of size N (observations for time t)

#This gives the smoothing coordinates of y (least squares)
S <- t(solve(t(Phi) %*% Phi) %*% t(Phi))

Y2CMap <- S %x% diag(rep(1,7)) #This mapping matrix returns vec(C) from vec(Y)


#C2BMap - This map takes smoothing coefficients C and returns smoothing coef B----
#This is where the regression model actually occurs
#The design matrix Z is involved in this step

Theta <- eval.basis(1:12, basis_B)
Z %*% B %*% t(Theta)


