function [z,w] = cc_rule(log_order)
% Generating Clenshaw-Curtis rule using the method of Jorg Waldvogel

N = 2^log_order;
z = cos((0:N)'*pi/N);

N2 = mod(N,2);
u0 = 1/(N^2-1+N2); % Boundary weights of CC
%
% Clenshaw-Curtis nodes: k = 0,1,...,N;
% vector of weights w_0 = w_n = u0
%
% auxiliary vectors
L = (0:N-1)';
m = min(L,N-L);
r = 2./(1-4*m.^2);
%
w = [ifft(r-u0); u0]; % Clenshaw-Curtis weights
end