function [Z0, alpha0] = init_SVT(A, k)
% Initialization by Singular Value Thresholding. Algorithm 3 in Ma and Ma, 2017
% without including covariate matrix X. 
%
% Function Prototype:
%   function [Z0, alpha0] = init_SVT(A, k)

% INPUT
%   A               binary adjacency matrix, shape = n * n
%   k               dimension of the latent vector, integer

% OUTPUT
%   Z0              initial estimate of the latent vector matrix, shape = n * k
%   alpha0          initial estimate of the degree heteroge

   % estimate logit matrix
   n = size(A, 1);
   logit = A2logit(A);
   logit = (logit + logit') / 2;

   % solving least squares
   logit_sum = sum(logit, 2);
   alpha0 = (logit_sum - 1 / 2  * ones(n, 1) * mean(logit_sum)) / n;

   % finding Z0
   R = logit - alpha0 * ones(1, n) - ones(n, 1) * alpha0'; 
   G0 = center(R);
   [U, D] = eig((G0 + G0') / 2);
   [D, order] = sort(diag(D), 'descend');  
   U = U(:, order);
   Z0 = U(:, 1:k) * sqrt(diag(D(1:k)));
end