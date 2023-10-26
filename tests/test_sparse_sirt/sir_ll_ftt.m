function [lF, lP] = sir_ll_ftt(betagamma, data, W, sigma_n, x0, tobs, masks)
if (nargin>6) && (~isempty(masks)) && (all(masks>0))
    data = data(masks);
    data = data(:);
end
betagamma = betagamma'; % FTT
d = size(W,1);
beta = betagamma(:,1:2:2*d-1)';
beta = reshape(beta, 1, d, []);
gamma = betagamma(:,2:2:2*d)';
gamma = reshape(gamma, 1, d, []);
x0 = repmat(x0, 1, size(beta,3));
[t,x] = ode45(@(t,x)sir_rhs(x,W,beta,gamma), tobs, x0, odeset('AbsTol', 1e-6, 'RelTol', 1e-6));
x = x(2:end, :);
x = reshape(x, [], 3*d, size(beta,3));
x = x(:, 2+(0:d-1)*3, :);
x = reshape(x, [], size(beta,3));

lF = ((x - data).^2)/(2*sigma_n^2); % - - for FFT
if (nargin<7) || (isempty(masks)) || (all(masks>0))
    lF = sum(lF, 1);
% else
    % DON'T sum over data! Return a matrix of pointwise likelihoods
end

lP = zeros(1, size(betagamma,1)); % FTT
end
