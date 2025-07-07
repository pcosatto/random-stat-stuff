library(fda)


# Extract monthly temperature curves (response)
temp <- t(as.matrix(data_split[[7]][,7:18]))  # 12 x 35 matrix

# Predictors
latitude <- data_split[[7]]$Latitude
region <- data_split[[7]]$Region
altitude <-  data_split[[7]]$Altitude

range_months <- c(0, 12)
nbasis_Y <- 7  # can be tuned

basis_Y <- create.fourier.basis(rangeval = range_months, nbasis = nbasis_Y)

# Time points (midpoint of each month)
month_points <- 1:12

# Smooth each curve using basis_Y
temp_fd <- Data2fd(argvals = month_points, y = temp, basisobj = basis_Y)

X <- model.matrix(~ latitude + region + altitude)
xfdlist <- lapply(1:ncol(X), function(j) X[, j])


nbasis_beta <- 5  # can be different from nbasis_Y
basis_beta <- create.fourier.basis(rangeval = range_months, nbasis = nbasis_beta)

lambda <- 1e-2  # smoothing parameter
beta_fdPar_list <- vector("list", ncol(X))
for (j in seq_along(beta_fdPar_list)) {
  beta_fdPar_list[[j]] <- fdPar(basis_beta, Lfdobj = int2Lfd(2), lambda = lambda)
}

model_fit <- fRegress(y = temp_fd, xfdlist = xfdlist, betalist = beta_fdPar_list)

beta_fd <- model_fit$betaestlist
# Example: plot latitude coefficient function
plot(beta_fd[[2]]$fd, main = "Functional coefficient for Latitude")

Yhat_fd <- model_fit$yhatfd
plot(temp_fd, col = "gray", main = "Observed vs Fitted")
lines(Yhat_fd, col = "blue")


