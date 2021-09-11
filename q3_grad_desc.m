function [beta,grad_desc] = q3_grad_desc(X,y,beta0,alpha,n)
% The following function does the gradient descent needed to do the
% regression 
% alpha is the learning rate parameter
% n is the number of iterations needed
% Initialize some useful values
m = length(y); % number of training examples
grad_desc = zeros(n, 1);

for iter = 1:n
    delta = ((X*beta0)-y)' * X;
    temp = beta0 - (alpha/m)*delta';
    beta0 = temp;
    grad_desc(iter) = q3_cost(X, y, beta0);
end
beta = beta0;
end

