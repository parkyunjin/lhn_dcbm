function [value] = center(X)
% center both the columns and rows of a matrix
     n = size(X, 1);
     X = X - ones(n, 1) * (ones(1, n) * X / n); 
     X = X - X / n * ones(n, 1) * ones(1, n);
     value = X;
end

     