function label = SCOREplus(A, k, c, r)
% This function is a community detection method called SCORE+
%
% Input: 
% A: adjacency matrix, which should be symmetric and connected.
% k: number of communities, which should be > 1.
%
% Optional input: 
% c: Tuning parameter for Graph Laplacian. Default is 0.1.
% r: Latent dimension (>1), if not given, chosen between k and k+1 determined by eigen gap
%
% Output:
% label: Predicted community labels. 

 if ~exist('c')
     % tuning parameter c does not exist, so default it to 0.1
      c = 0.1;
 end

 if ~exist('r')
     % if r not give, set to be k+1
     r = k+1;
     fix.latent.dim = false;
 else
     fix.latent.dim = true;
 end
 
n = size(A,1);	             % Number of nodes
G = graph(A);

if max(conncomp(G))>1        % Check connectivity
    fprintf('The adjacency matrix is not connected! \n');
    return                      
end
if ~issymmetric(A)           % Check symmetry
    fprintf('The adjacency matrix is not symmetric! \n');
    return                      
end
    
if k < 2                     % Check whether k>1
    fprintf('The number of communities k should be >1!\n');  
    return
end       
if c < 0                     % Check whether c>0
    fprintf('The tuning parameter c should be positive!\n');  
    return
end     
delta = c*max(sum(A));     % tuning parameter for Graph Laplacian


d = 1./sqrt(delta+sum(A));

A2 = A.*repmat(d', [1 n]);
A1 = (A2)'.*repmat(d', [1 n]);    % A1 is the Graph Laplacian with tuning c.


[vb,db] = eigs(A1,r);        % PCA

vb2 = vb*db;                 % Re-weighting
vb3 = vb2(:,2:r);
rb = vb3./repmat(vb2(:,1),[1,r-1]); % eigen wise ratios


% decide latent dimension by eigen-gap.
if ~fix.latent.dim
    temp = (db(k,k) - db(k+1,k+1))/db(k,k); 
    if temp > 0.1
            rb = rb(:,1:(k-1));
    end
end

% use k-means to predict labels
label = kmeans(rb,k,'replicates',1000);
end

