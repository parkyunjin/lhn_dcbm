function [error] = clustering_error_Z(Z, label)
% Computing the clustering error. 
%
% Function Prototype:
%   function [error] = clustering_error_perm(id, label, k)

% INPUT
%   Z               latent vector matrix, shape = n * p
%	label    		true label, length = n
%   k               number of clusters, integer 

% OUTPUT
%   error           clustering error

n = length(label);
labelset = unique(label);
k = length(labelset);

[U, D, V] = svd(Z, 0);
Z_k = U(:,1:k) * D(1:k,1:k);
[id, C] = kmeans(Z_k, k, 'Replicates', 200);
idset = unique(id);

P = perms(labelset);
N = size(P, 1);
confusion = zeros(1, N);
for i = 1 : N
    for node = 1 : n
        if (label(node) ~= P(i, id(node)))
               confusion(i) = confusion(i) + 1;
        end
    end
end
error = min(confusion) / n;