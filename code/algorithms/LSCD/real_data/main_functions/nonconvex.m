function [Z1, alpha1] = nonconvex(A, Z0, alpha0, varargin)
% Non-convex estimation of the latent space model with adjacency matrix A. 
% Algorithm 1 in Ma and Ma, 2017 without covariate matrix X. 
%
% Function Prototype:
%   function [Z1, alpha1] = nonconvex_method_covariate(A, Z0, alpha0, 
%                                                      maxiter, report_every)
%
% INPUT
%   A               binary adjacency matrix, shape = n * n
%   Z0              initial estimate of the latent vector matrix, shape = n * k
%   alpha0          initial estimate of the degree heterogeneity parameter, length = n
%
% OPTIONAL INPUT
%   maxiter         maximum number of iterations, integer, default = 2000
%   report_every    report change of estimates every report_every steps, integer
%                   default = 500
%   miniter         minimum number of iterations, integer, default = 100
%   epsilon         target precision of this optimization algorithm, real number
%                   default = 10^(-5)

%     
% OUTPUT
%   Z1          final estimate of the latent vector matrix, shape = n * k
%   alpha1      final estimate of the degree heterogeneity parameter, length = n

numvarargs = find(~cellfun('isempty',varargin));
optargs = {2000, 500, 100, 10^(-5)};
optargs(numvarargs) = varargin(numvarargs);
[maxiter, report_every, miniter, epsilon] = optargs{:};

n = length(alpha0);
%flip A
A_inv = 1.-A;  
A_inv(1:length(A_inv)+1:numel(A_inv)) = 0;

% step size
eta_alpha = 2.0 / n; 
eta_Z = 1.0 / (norm(Z0, 2))^2;

iter = 0;
G0 = Z0 * Z0';
change_alpha = 1.0;
change_G = 1.0;


while (((max(abs(change_alpha), abs(change_G)) > epsilon) && (iter<maxiter))||(iter<miniter))
    iter = iter + 1;
    logit0 = G0 + alpha0 * ones(1, n) + ones(n, 1) * alpha0';
    E = sigmoid(logit0);
    Z1 = Z0 - eta_Z * (E - A) * Z0;
    Z1 = Z1 - ones(n,1) * (ones(1,n) * Z1) / n;  
    alpha1 = alpha0 - eta_alpha * ((E - A) * ones(n,1)); 
    G1 = Z1 * Z1';
    if mod(iter, report_every) == 0
        change_alpha = norm(alpha1 - alpha0, 2) / norm(alpha0, 2);
        change_G = norm(G1 - G0,'fro')/norm(G0,'fro');
        fprintf('iteration = %d\n', iter);        
        fprintf('change-G = %f\n', change_G) 
        fprintf('change-alpha = %f\n', change_alpha)
    end  
    Z0 = Z1;
    G0 = G1;    
    alpha0 = alpha1;
end 
end