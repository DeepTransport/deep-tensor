function [x,err,qA,rA] = ls_solve(A, y, w)

% A num_samples x cardinality
% y num_samples x num_outputs
%

% m: number of samples, n: number of coefficients
[m,n] = size(A);
nout = size(y,2);

if ~isempty(w)
    A = A.*sqrt(w);
    y = y.*sqrt(w);
end

[qA,rA] = qr(A,0);
x = rA\(qA'*y);

second_moment = sum(y.^2,1)/size(y,1);

% Compute the relative leave-one-out cross-validation error the fast
% leave-one-out cross-validation procedure [Cawlet & Talbot 2004] based on
% the Sherman-Morrison-Woodbury formula
%

err = zeros(1,nout);

if m-1 < n
    warning('Not enough samples for cross-validation');
    err(:) = Inf;
    return
end

% Compute the predicted residuals using the Bartlett matrix inversion formula
% (special case of the Sherman-Morrison-Woodbury formula)
T = sum(qA.^2,2);
del = (y-A*x)./(1-T);
% Compute the absolute cross-validation error
err = sum(del.^2,1)/m;

% Compute the relative cross-validation error
err = err./second_moment;


if rcond(rA) > 1E-16
    %invrA = inv(rA);
    %C = invrA*invrA';
    % (1-n/m)^(-1)*(1+trace(C))
    % Direct Eigenvalue Estimator [Chapelle, Vapnik & Bengio, 2002], [Blatman & Sudret, 2011]
    % -> accurate even when m is not >> n (corr ~ 1+2*n/m when m->Inf, as trace(C) ~ n/m when m->Inf)
    %corr2 = (m/(m-n))*(1+trace(C));
    %norm(sum(diag(rA).^(-1))^2 - trace(C))
    %norm(sum(diag(rA).^(-1)) - trace(invrA))
    %norm(corr-corr2)
    corr = (m/(m-n))*(1+sum(diag(rA).^(-1))^2);
    if corr > 0
        err = err*corr;
    end
end

err = sqrt(err);

end