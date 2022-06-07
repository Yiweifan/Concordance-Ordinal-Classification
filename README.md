# Concordance-Ordinal-Classification
Concordance-Ordinal-Classification (COC) implements the method from the paper "Sparse Concordance-based Ordinal Classification".

## Usage

- `cv.cocnet(trainx, trainy, K, weights, lambda, nfolds, foldid, ini, no.it=10)`
  - **Description:** Does k-fold cross-validation for a COC model with lasso regularization at a grid of values for the regularization parameter.
  - **Arguments:**
    - <font face="逐浪圆体">trainx</font>: Input matrix, of dimension nobs x nvars; Each row is an observation vector.
    - trainy: Response variable. Should be a vector of integers.
    - K: Number of categories.
    - weights: Separate penalty factors that can be applied to each coefficient.
    - lambda: A user supplied lambda sequence.
    - nfolds: Number of folds.
    - foldid: A list of fold position integers corresponding to the training data.
    - ini: 	Initial values for coefficients.
    - no.it: The maximum number of iterations. Defaults to 10.
  - **Value:**
    - lambda: Sequence of lambda values
    - cvm: The mean cross-validated error - a vector of length length(lambda).

- `coc(trainx, trainy, K, weights, lambda, loss=c("01","abs","rps","all"), rescale=1, ini, no.it=10)`
  - **Description:** Fit a penalized smoothed COC model. The regularization path is computed for the lasso penalty at a given value for the regularization parameter.
  - **Arguments:**
    - trainx: Input matrix, of dimension nobs x nvars; Each row is an observation vector.
    - trainy: Response variable. Should be a vector of integers.
    - K: Number of categories.
    - weights: Separate penalty factors that can be applied to each coefficient.
    - lambda: A user supplied lambda.
    - loss: Loss function to mimize. The default "all" gives results of all three loss functions.
    - rescale: The constraint of the norm of coefficients. Defaults to 1.
    - ini: 	Initial values for coefficients.
    - no.it: The maximum number of iterations. Defaults to 10.
  - **Value:**
    - coefficients: A vector of coefficients.
    - thresholds: A list of thresholds.

