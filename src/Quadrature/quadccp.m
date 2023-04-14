function [y,fcnt,z,w,mid,dxdz] = quadccp(fun,theta,a,b,tol)
% Adaptive numerical integration using Clenshaw-Curtis,
% integrates the function from a to b.
%   [y,fcnt,refx,refw,] = QUADCCP(fun,theta,a,b,tol)
%
%   fun     - given as either a string or an inline function
%   theta   - parameters that defines the function
%   a       - left boundary, 1 x m vector
%   b       - right boundary, 1 x m vector
%   tol     - tolerance, default is 1E-10
%   y       - the result of the integration
%   fcnt    - the number of function evaluations
%   z       - quadrature points used, in the reference interval [-1,1]
%   w       - reference quadrature weights
%
% Example:
%   fun = @(c,x) 1./(x.^3-2*x-c); % define function
%   [y, n] = quadcc(@(theta,x) fun(theta,x), 5, 0, 2) % integrate
%
%   % We can also return the points where the function is evaluated at
%   % and the corresponsing weights
%   [y, n, x, w, jac] = quadcc(@(theta,x) fun(theta,x), 5, 0, 2);
%
%   % Vector valued boundaries
%   [y, n, x, w] = quadcc(@(x) fun(x,5), [0,-20], [2,2]);
%   % User specified tolerance (note: the following does not converge)
%   [y, n] = quadcc(@(x) fun(x,5), [0,-20], [2,2], 1E-16)
%
% Reference: 
%   Jörg Waldvogel. Fast Construction of the Fejer and Clenshaw-Curtis 
%   Quadrature Rules. BIT Numerical Mathematics 46 (2006), 195-202.

max_log_order = 16;

if nargin < 3
    error('Need to specify boundary points')
elseif nargin == 3
    tol = 1E-10;
end
m = length(a);
if m ~= length(b)
    error('Boundary points have mismatch dimensions')
end
% change of variable: x in [a, b] and z in [-1, 1]
%   z = 2 (x-a)/(b-a) - 1
%   x = z (b-a)/2 + (a+b)/2
dxdz = reshape((b-a)/2, 1, []); % dx/dz
%
mid = reshape((b+a)/2, 1, []); % shift by the mid point

y = 0;
conv_flag = false;
for lo = 1:max_log_order
    yp = y;
    %
    [z,w] = cc_rule(lo);
    if lo == 1
        x = z.*dxdz + mid;
        f_cache = reshape(feval(fun,theta,x),[],m);
    else
        n = 2^lo + 1;
        ind1 = 1:2:n;
        ind2 = 2:2:n;
        x = z(ind2).*dxdz + mid;
        f_cache([ind1,ind2],:) = cat(1, f_cache, reshape(feval(fun,theta,x),[],m));
    end
    y = sum(f_cache.*w,1).*dxdz;
    %
    if lo > 1 && norm(y - yp, Inf) < tol
        conv_flag = true;
        break;
    end
end

fcnt = (2^lo+1)*m;

if ~conv_flag
    warning(['CC quad does not converge, final error = ' num2str(norm(y-yp, Inf))])
end

end
