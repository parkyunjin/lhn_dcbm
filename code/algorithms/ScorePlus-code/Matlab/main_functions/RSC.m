function label = RSC(A, k, c)
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


[vb,db] = eigs(A1, k);        % PCA


%rb = diag(sqrt(sum(vb.*vb)))\vb;
norms = sqrt(sum(vb.*vb));
rb = vb ./ norms;

% use k-means to predict labels
label = kmeans(rb,k,'replicates',1000);
end

