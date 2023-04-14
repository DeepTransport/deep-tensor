function [x,err,qA,rA] = ls_solve_update(qA, rA, Aadd, y, w)

% A num_samples x cardinality
% y num_samples x num_outputs
%

% m: number of samples, n: number of coefficients
[m,n] = size(Aadd);
nout = size(y,2);

if ~isempty(w)
    Aadd = Aadd.*sqrt(w);
    y = y.*sqrt(w);
end

R1 = qA'*Aadd;
[Q2, R2] = qr(Aadd - qA*R1, 0);

rA = [rA, R1; zeros(size(Aadd,2),size(qA,2)), R2];
qA = [qA, Q2];

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
del = (y-qA*(rA*x))./(1-T);
% Compute the absolute cross-validation error
err = sum(del.^2,1)/m;

% Compute the relative cross-validation error
err = err./second_moment;


if rcond(rA) > 1E-16
    corr = (m/(m-n))*(1+sum(diag(rA).^(-1))^2);
    if corr > 0
        err = err*corr;
    end
end

err = sqrt(err);

end