function [Z0, alpha0, beta0] = init_SVT_X(A, k, X)
% Initialization by Singular Value Thresholding. Algorithm 3 in Ma and Ma, 2017
%
% Function Prototype:
%   function [Z0, alpha0, beta0] = init_SVT_X(A, k, X)

% INPUT
%   A               binary adjacency matrix, shape = n * n
%   k               dimension of the latent vector, integer
%   X               observed symmetric feature matrix X, shape = n * n

% OUTPUT
%   Z0              initial estimate of the latent vector matrix, shape = n * k
%   alpha0          initial estimate of the degree heteroge
%   beta0           initial estimate of the homophily effect, real number

   % estimate logit matrix
   n = size(A, 1);
   logit = A2logit(A);
   logit = (logit + logit') / 2;

   %solving least squares
   logit_sum = sum(logit, 2);
   X_sum = sum(X, 2);
   u = (logit_sum - 1 / 2  * ones(n, 1) * mean(logit_sum))/n;
   v = (X_sum - 1 / 2  * ones(n, 1) * mean(X_sum)) / n;
   c = 1 / (norm(X, 'fro')^2 / 2 - X_sum' * v);
   logit_X = sum(sum(X .* logit));
   alpha0 = u + c * v * (v' * logit_sum) - c / 2 * logit_X * v;
   beta0 = - c * v' * logit_sum + c / 2 * logit_X;

   %finding Z0
   R = logit - alpha0 * ones(1, n) - ones(n, 1) * alpha0' - beta0 * X; 
   G0 = center(R);
   [U, D] = eig((G0 + G0') / 2);
   [D, order] = sort(diag(D), 'descend');  
   U = U(:, order);
   Z0 = U(:, 1:k) * sqrt(diag(D(1:k)));
end