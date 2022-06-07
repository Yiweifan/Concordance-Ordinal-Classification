# Concordance-Ordinal-Classification
Concordance-Ordinal-Classification (COC) implements the method from the paper "Sparse Concordance-based Ordinal Classification".

## Usage

- `coc(trainx, trainy, K, weights, lambda, loss=c("01","abs","rps","all"), rescale=1, ini, no.it=10)`
  - **Description:** Fits a COC model with lasso penalty at a given value for the regularization parameter.
  - **Arguments:**
    - trainx: Input matrix, of dimension nobs x nvars; Each row is an observation vector.
    - trainy: Response variable. Should be a vector of integers.
    - K: Number of categories.
    - weights: Separate penalty factors that can be applied to each coefficient.
    - lambda: A user supplied lambda.
    - loss: Loss function to mimize. The default "all" minimizes all three loss functions.
    - rescale: The constraint of the norm of coefficients. Defaults to 1.
    - ini: 	Initial values for coefficients.
    - no.it: The maximum number of iterations. Defaults to 10.
  - **Value:**
    - coefficients: A matrix of coefficients.
    - thresholds: A list of thresholds.

- `cv.cocnet(trainx, trainy, K, weights, lambda, nfolds, foldid, ini, no.it=10)`
  - **Description:** Does k-fold cross-validation for a COC model with lasso penalty at a grid of values for the regularization parameter.
  - **Arguments:**
    - trainx: Input matrix, of dimension nobs x nvars; Each row is an observation vector.
    - trainy: Response variable. Should be a vector of integers.
    - K: Number of categories.
    - weights: Separate penalty factors that can be applied to each coefficient.
    - lambda: A user supplied lambda sequence.
    - nfolds: Number of folds.
    - foldid: A list of fold position integers corresponding to the training data.
    - ini: 	Initial values for coefficients.
    - no.it: The maximum number of iterations. Defaults to 10.
  - **Value:**
    - lambda: Sequence of lambda values.
    - cvm: The mean cross-validated error - a vector of length length(lambda).

- `evaluation(predy, testy, cumprob, K)`
  - **Description:** Evaluates ordinal classification performance.
  - **Arguments:**
    - predy: A vector of predicted labels.
    - testy: A vector of true labels.
    - cumprob: A matrix of dimension K x nobs of predicted cumulative probabilities.
    - K: Number of categories.
  - **Value:** A vector of misclassification error rate, mean absolute error, Kendall's coefficient, and mean ranked probability score.

- `mycut(beta, test_X, cutoff)`
  - **Description:** Predict method for a coc model.
  - **Arguments:**
    - beta: A matrix of coefficients.
    - test_X: Covariate matrix.
    - cutoff: A vector of thresholds.
  - **Value:**
    - predy: A vector of predicted labels.
    - prop: Proportion of non-monotone cases for out-of-sample in test_X.

- `solu(initial, train_X, train_Y, k, Lambda, w, no.it=10, mu=0)`
  - **Description:** Maximize a penalized smoothed concordance function for category k with lasso or ridge penalty at a given value for the regularization parameter.
  - **Arguments:**
    - initial: 	Initial values for coefficients.
    - train_X: Input matrix, of dimension nobs x nvars; Each row is an observation vector.
    - train_Y: Response variable. Should be a vector of integers.
    - k: A user supplied category.
    - Lambda: A user supplied lambda for the lasso penalty.
    - w: Separate penalty factors that can be applied to each coefficient.    
    - no.it: The maximum number of iterations. Defaults to 10.
    - mu: A user supplied mu for the ridge penalty.
  - **Value:**
    - beta: A vector of coefficients.
 
## Examples
 
The Balance Scale Data is available from the UCI repository. It was generated from psychological experiments, where each sample was classified as having a balance scale tip to the right, tip to the left, or being balanced. The attributes included the left weight, the left distance, the right weight, and the right distance.


