function [z, mlogf] = eval_rt(obj, x, k)
% Evaluate the deep Rosenblatt transport Z = T(X), where Z is the reference
% random variable and X is the target random variable.
%   [z, f] = EVAL_RT(dirt, x, k)
%
%   x    - random variable drawn form the pdf defined by DIRT
%   k    - number of layers in the evaluations
%   z    - reference random variables, d x n
%   mlogf- negative log of the DIRT density

if nargin <= 2
    k = num_layers(obj);
else
    k = min(k, num_layers(obj));
end
z = x;
mlogf = zeros(1,size(z,2));
for l = 1:k
    if obj.basis_r(l) > 0 && obj.basis_r(l) < obj.d
        zr = obj.basis{l}'*z;
        zn = z - obj.basis{l}*zr;
        %
        % evaluate rt
        ur = eval_rt(obj.irts{l}, zr);
        mlogt = eval_potential(obj.irts{l}, zr);
        % evaluate the diagonal transform
        zr = invert_cdf(obj.ref, ur);
        logd = log_joint_pdf(obj.ref, zr);
        %
        z = obj.basis{l}*zr + zn;
    elseif obj.basis_r(l) == obj.d
        % reparametrize, double check
        z = obj.basis{l}'*z;
        % evaluate rt
        u = eval_rt(obj.irts{l}, z);
        % eval pdf (det Jacobian)
        mlogt = eval_potential(obj.irts{l}, z);
        % evaluate the diagonal transform
        z = invert_cdf(obj.ref, u);
        logd = log_joint_pdf(obj.ref, z);
    else
        % evaluate rt
        u = eval_rt(obj.irts{l}, z);
        % eval pdf (det Jacobian)
        mlogt = eval_potential(obj.irts{l}, z);
        % evaluate the diagonal transform
        z = invert_cdf(obj.ref, u);
        logd = log_joint_pdf(obj.ref, z);
    end
    % update density
    mlogf = mlogf + mlogt + logd;
end
logf = log_joint_pdf(obj.ref, z);
mlogf = mlogf - logf;

end
