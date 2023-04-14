function J = eval_rt_jac_reference(obj, r, z)
% Evaluate the jacobian of the squared Rosenblatt transport Z = R(X), where
% Z is the uniform random variable and X is the target random variable.
%   J = EVAL_RT_JAC_REFERENCE(irt, X, Z)
%
%   X - random variable drawn form the pdf defined by SIRT
%   Z - uniform random variables, d x n
%   J - Jacobian, d x (d x n), each d x d block is the Jabocian for X(:,j)

d = ndims(obj.approx);
n = size(r,2);
if size(r,1) ~= d
    error('Input variable must be d x n')
end
J = zeros(d,n*d);

nbatch = ceil(n/obj.batch_size);
start_i = 1:obj.batch_size:n;
end_i = start_i+(obj.batch_size-1);
end_i(end) = n;

for batch = 1:nbatch
    ind = start_i(batch):1:end_i(batch);
    ni = numel(ind);
    jnd = (1:(d*ni)) + (batch-1)*(d*obj.batch_size);
    J(:,jnd) = eval_rt_jac_reference_inner(obj, r(:,ind), z(:,ind));
end

end