# Concordance-Ordinal-Classification
Concordance-Ordinal-Classification (COC) is a R framework that implements the method from the paper "Sparse Concordance-based Ordinal Classification".

cv.cocnet
fit a COC model with lasso regularization

Description
Does k-fold cross-validation for coc

Fit a penalized smoothed concordance ordinal classification model. The regularization path is computed for the lasso penalty at a grid of values for the regularization parameter.

Usage
cv.cocnet(trainx,trainy,K,weights,lambda,nfolds,foldid,ini,no.it=10)

Arguments
trainx input matrix, of dimension nobs x nvars; each row is an observation vector.
trainy response variable
K number of categories
weights Separate penalty factors can be applied to each coefficient
lambda A user supplied lambda sequence
nfolds number of folds
foldid an optional vector of values between 1 and nfold identifying what fold each observation is in.

