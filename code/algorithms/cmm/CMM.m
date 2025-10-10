function label = CMM(A, r, lambda)
% Convexified Modularity Maximization
%
% [idx, result] = cmm(A, r, lambda)
% Input:
%   A      - Adjacency matrix of the graph
%   r      - Number of clusters to extract
%   lambda - (Optional) Tuning parameter
% Output: 
%   label  - Cluster indices of each node
%
% Reference:
%   "Convexified Modularity Maximization for Degree-corrected Stochastic Block Models"
%   http://arxiv.org/abs/1512.08425
%   Yudong Chen (Cornell University; yudong.chen@cornell.edu)
%   Xiaodong Li (University of California, Davis; xdgli@ucdavis.edu) 
%   Jiaming Xu  (University of Pennsylvania; jiamingx@wharton.upenn.edu)


%% Algorithm parameters
iter  = 200;		% Max number of iterations for ADMM
replicates = 100;	% Number of random initializations for k-median

%%
n = size(A,1);	% Number of nodes
d = sum(A,2);	% Column vector of node degrees
if nargin < 3   % If no tuning parameter specified
	mu = 1/sum(d); % This corresponds to the modularity maximization
else
	mu = lambda;    %This corresponds to generalized mm  
end


%% Solve convex relaxation by ADMM

E = - A + mu*d*d';
Z=zeros(n);
Lambda=zeros(n);
k=0;
rho=1;
while k<iter
    k=k+1;
    ZZ=(Z-Lambda-E/rho)/2+(Z-Lambda-E/rho)'/2;
    [U, D]=eig(ZZ);
    Y = U*max(D, 0)*U';
    Z=min(max(Y+Lambda, 0), 1);   %Project all entries into [0,1].
    Z(1:n+1:end)=1;               %Change all diagonal entries into 1.
    Lambda = Lambda+Y-Z;
end


%% Weighted L1-norm k-medoid (TBD)

feature_weights = max(d,2);
weights = max(d,2);

% Pairwise weighed distances
D = squareform( pdist(Y*diag(feature_weights), 'cityblock') );

% Run the greedy algorithm from multiple random initialiozations,
% and choose the best result
label = ones(n,1); cost = 1/eps;

for t = 1:max(replicates,1)
    % Random initial r medoids (cluster centers)
    % Assign nodes to its closest cluster center
    % label_t[i] is the cluster index of node i
    [~, label_t] = min(D(randsample(n,r),:),[],1); 
    
    last = 0;
    while any(label_t ~= last) % While the clustering has been changed
        % Assign as new cluster center node which has the minimum sum 
        % of the weighted L_1 distances to points in cluster k 
        % index_t[k] is the index of the new center for cluster 1<=k<=r
        [~, index_t] = min(D*diag(weights)*sparse(1:n,label_t,1,n,r,n),[],1);
        last = label_t;
        % recluster nodes according to the new cluster centers
        [val, label_t] = min(D(index_t,:),[],1);
    end
    cost_t = sum(val);
    
    % Choose the best result
    if cost_t < cost
        label  = label_t;
        cost = cost_t;
    end
    
end
