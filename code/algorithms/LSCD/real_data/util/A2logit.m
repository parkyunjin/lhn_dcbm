function [logit] = A2logit(A, varargin)
% Estimate the logit matrix by singular value thresholding.
%
% Function Prototype:
%   function [logit] = A2logit(A, varargin)

% INPUT
%   A               binary adjacency matrix, shape = n * n

% OPTIONAL INPUT
%   tau             cutoff of the singular value thresholding
%                   default = sqrt(n * p_ave) where p_ave = sum(sum(A))/(n*(n-1))
%   M               -M specifies the lower bound of the logits
%                   default = 4

% OUTPUT
%   logit           estimated logit matrix, shape = n * n

n = size(A, 1);
p_ave = sum(sum(A))/(n*(n-1));
tau0 = sqrt(n * p_ave);
numvarargs = find(~cellfun('isempty',varargin));
optargs = {tau0, 4};
optargs(numvarargs) = varargin(numvarargs);
[tau, M] = optargs{:};

    [U,D,V] = svd(A);
    d = diag(D);
    d = d.*(d > tau);
    D = diag(d);
    temp = U * D * V';
    temp = max(temp, 1/(1 + exp(M)));
    temp = min(temp, 1/(1 + exp(-M)));
    logit = -log(1 ./ temp - 1);
end