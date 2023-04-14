function [r,f] = eval_cirt_reference(obj, x, z)
% Evaluate the inverse of the conditional squared Rosenblatt transport 
% Y | X = R^{-1}(Z, X), where X is given, (X, Y) jointly follow the target 
% distribution represented by SIRT, Z is uniform. 
%   [Y,f] = EVAL_CIRT(irt, X, Z)
%
%   Z - uniform random variables, m x n
%   X - given variables following the marginal of the target distribution
%   Y - random variable drawn from Y | X
%   f - pdf of Y | X
%
% The conditioning depends on irt.order.
%   * With >0 (marginalised from right to left), X = (X_1, ..., X_k) is the
%     the first k coordinates in the joint and Y = (X_{k+1}, ... X_d)
%   * With <0 (marginalised from left to right), X = (X_m, ..., X_d) is the
%     last (d-m+1) coordinates in the joint and Y = (X_1, ... X_{m-1})

r = x;
f = z;

end
