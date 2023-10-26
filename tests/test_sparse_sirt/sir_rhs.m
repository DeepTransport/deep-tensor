function [f]=sir_rhs(x,W,beta,gamma)
% Compartmental SIR:
% dS_k/dt = -beta_k S_k I_k + W_{kj} S_j
% dI_k/dt = beta_k S_k I_k - gamma_k I_k + W_{kj} I_j
% dR_k/dt = gamma_k I_k + W_{kj} R_j
%
% W*1 = 0, e.g. W = circ(-1,1) or W = tridiag(1,-2,1)
% 
% x is of size 3d x I
% beta,gamma are of size 1 x d x I
d = size(W,1);
x = reshape(x, 3, d, []);
f = permute(x, [2,1,3]);
f = reshape(f, d, []);
f = W*f; % all movements are done here
f = reshape(f, d, 3, []);
f = permute(f, [2,1,3]);
f(1,:,:) = f(1,:,:) - beta.*x(1,:,:).*x(2,:,:);
f(2,:,:) = f(2,:,:) + beta.*x(1,:,:).*x(2,:,:) - gamma.*x(2,:,:);
f(3,:,:) = f(3,:,:) + gamma.*x(2,:,:);
f = f(:);
end

