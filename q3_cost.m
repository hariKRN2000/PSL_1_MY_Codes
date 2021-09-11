function q3_cost = q3_cost(X,y,beta0) ; 
% This function computes cost function
% Initialize some useful values
m = length(y); % number of training examples
predictions = X*beta0;
sqrErrors = (predictions - y).^2 ; 
q3_cost = (1/(2*m))*sum(sqrErrors) ; 
end