d = 5; a = 0.5;
A = diag(-sqrt(1-a^2)*ones(d-1,1), -1) + eye(d);
D = diag([1, a*ones(1,d-1)]);
B = D\A;
Q = B'*B;   % precision matrix
C = inv(Q); % covariance matrix
z = sqrt((2*pi)^d/det(Q)); % normalising constant
% The joint distribution, unnormalised


nsteps = 1E4;
tic;
out = NUTS(@(x) log_target_nuts(B,Q,x), randn(d,1), nsteps);
toc

tic;
out = pCN(@(x) log_target_pcn(B,x), randn(d,1), nsteps, 1);
toc

norm(cov(out.samples') - C, 'fro')/norm(C, 'fro')
mean(out.samples')

%%%%%

function [f,g] = log_target_nuts(B,Q,x)

f = 0.5*sum((B*x).^2,1);
g = Q*x;

end

function [f1,f2] = log_target_pcn(B,x)

f1 = 0.5*sum((B*x).^2,1);
f2 = 0.5*sum(x.^2,1);
end