function [value] = sigmoid(x)
% compute the sigmoid function
	value = 1 ./ (1 + exp(-x));
end
