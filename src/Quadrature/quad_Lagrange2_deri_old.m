function [f,g] = quad_Lagrange2_deri_old(poly, theta, x)
% compute the integration of the derivative of piecewise Lagrange 2
% polynomials from poly.domain(1) to x
%
%   [f,g] = QUAD_LAGRANGE1_DERI(poly, theta, x)
%
%   theta - interpolation function values, d x n, where n is the number of
%           functions and d is the number of interpolation
%   x     - 1 x n, number of terminal points
%   f     - integration result
%   g     - gradient of the integration with respect to theta


% the function to be integrated log(2.^b + 1)/log(2); with piecewise b
% derivative of the function b./log(2.^b + 1)
%
% size of each element is given by poly.elem_size, except ghost cells
%
% the derivative of ghost cells are non-zero if 'Dirichlet' bc is used
% the coefficient is given by the following
%
% For lagrange 2: determine the function a+bx+cx^2, then the derivative is
% piecewise linear with value b+2cx
%

% L = [1, 0, 0; 1, 1, 0; 1, 2, 1];
% U = [1, 0, 0; 0, dhh, dhh^2; 0, 0, dhh^2*2];
% LU = Vandermore
% Vandermore = [1, 0, 0; 1, dhh, dhh^2; 1, dhh*2, 4*dhh^2];
% inverse of the Vandermore matrix
dhh = poly.elem_size/2;
iV = [1, 0, 0; -3/(2*dhh), 2/dhh, -1/(2*dhh); 1/(2*dhh^2), -1/dhh^2, 1/(2*dhh^2)];

if poly.order ~= 2
    error('polynomial order should be 2')
end

% index in which the boundary x falls into
% ei = 0: left ghost cell, ei = poly.num_elems+1: right ghost cell
ei = ceil( (x-poly.grid(1))./poly.elem_size );
%
f = zeros(size(x));
g = zeros(size(theta));
for k = 0:(poly.num_elems+1)
    ind1 = (ei>k);  % integrate the whole element
    ind2 = (ei==k); % integrate to x
    n1 = sum(ind1);
    n2 = sum(ind2);
    %
    if k == 0 && poly.gs > 0 && strcmp(poly.bc, 'Dirichlet')
        %b = theta(1,:)./poly.gs;
        %tmp = 2.^b;
        tmp = 2.^(theta(1,:)./poly.gs);
        if n1 > 0
            f(ind1) = f(ind1) + log2(tmp(ind1)+1)*poly.gs;
            g(1,ind1) = g(1,ind1) + tmp(ind1)./(tmp(ind1)+1);
        end
        if n2 > 0
            dx = x(ind2) - poly.domain(1);
            f(ind2) = f(ind2) + log2(tmp(ind2)+1).*dx;
            g(1,ind2) = g(1,ind2) + tmp(ind2)./(tmp(ind2)+1).*(dx/poly.gs);
        end
    elseif k == (poly.num_elems+1) && poly.gs > 0 && strcmp(poly.bc, 'Dirichlet')
        %b = -theta(end,:)./poly.gs;
        %tmp = 2.^b;
        tmp = 2.^(-theta(end,:)./poly.gs);
        if n1 > 0
            f(ind1) = f(ind1) + log2(tmp(ind1)+1)*poly.gs;
            g(end,ind1) = g(k,ind1) - tmp(ind1)./(tmp(ind1)+1);
        end
        if n2 > 0
            dx = x(ind2) - poly.grid(end);
            f(ind2) = f(ind2) + log2(tmp(ind2)+1).*dx;
            g(end,ind2) = g(end,ind2) - tmp(ind2)./(tmp(ind2)+1).*(dx/poly.gs);
        end
    elseif k >=1 && k <= poly.num_elems
        ind = (1:3) + 2*(k-1);
        %coef = iV*theta(ind,:); % a b c for each elements
        b = iV(2,:)*theta(ind,:);
        c = iV(3,:)*theta(ind,:);
        % int(g) = -dilog(2^(b + 2*c*x) + 1)/(2*c*log(2)^2)
        % dint(g)/db = log(2^(b + 2*c*x) + 1)/(2*c)
        % dint(g)/dc = dilog(2^(b + 2*c*x) + 1)/(2*c^2*log(2)) + (x*log(2^(b + 2*c*x) + 1))/c
        if n1 > 0
            tmp = dilog2(-2.^b(ind1)) - dilog2(-2.^(b(ind1)+2*c(ind1)*poly.elem_size));
            %tmp3 = polylog(2,-2.^b(ind1)) - polylog(2,-2.^(b(ind1)+2*c(ind1)*poly.elem_size));
            %norm(tmp - tmp3)
            %tmp2 = (dilog(2.^b(ind1)+1) - dilog(2.^(b(ind1)+2*c(ind1)*poly.elem_size)+1))./(2*c(ind1)*log(2)^2);
            %norm(tmp - tmp2)
            f(ind1) = f(ind1) + tmp./(2*c(ind1)*log(2)^2);
            %g(end,ind1) = g(k,ind1) - tmp(ind1)./(tmp(ind1)+1);
            %
            dintdb = (log(2.^(b(ind1)+2*c(ind1)*poly.elem_size)+1)-log(2.^b(ind1)+1))./(2*c(ind1)*log(2));
            dintdc = -tmp./(2*c(ind1).^2*log(2)^2) + poly.elem_size*log(2.^(b(ind1)+2*c(ind1)*poly.elem_size)+1)./(c(ind1)*log(2));
            g(ind,ind1) = g(ind,ind1) + dintdb.*iV(2,:)' + dintdc.*iV(3,:)';
        end
        if n2 > 0
            dx = x(ind2) - poly.grid(k);
            tmp = dilog2(-2.^b(ind2)) - dilog2(-2.^(b(ind2)+2*c(ind2).*dx));
            %tmp3 = polylog(2,-2.^b(ind2)) - polylog(2,-2.^(b(ind2)+2*c(ind2).*dx));
            %norm(tmp - tmp3)
            %tmp2 = (dilog(2.^b(ind2)+1) - dilog(2.^(b(ind2)+2*c(ind2).*dx)+1))./(2*c(ind2)*log(2)^2);
            %norm(tmp - tmp2)
            f(ind2) = f(ind2) + tmp./(2*c(ind2)*log(2)^2);
            %g(end,ind2) = g(end,ind2) - tmp(ind2)./(tmp(ind2)+1).*(dx/poly.gs);
            %
            dintdb = (log(2.^(b(ind2)+2*c(ind2).*dx)+1)-log(2.^b(ind2)+1))./(2*c(ind2)*log(2));
            dintdc = -tmp./(2*c(ind2).^2*log(2)^2) + dx.*log(2.^(b(ind2)+2*c(ind2).*dx)+1)./(c(ind2)*log(2));
            g(ind,ind2) = g(ind,ind2) + dintdb.*iV(2,:)' + dintdc.*iV(3,:)';
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = dilog2(z)
% DILOG = di-Logarithm.
%
% d = dilog(z) = Li_2(z) 
%   = -Int From t=0 To t=z    log(1-t) dt/t         for all z.
%   =  Sum From n=1 To n=Inf  z^n/n^2               for |z|<=1.
%
% INPUT  z: real or complex, scalar, vector or matrix.
% OUTPUT d: component-wise dilogarithm of z.
% References:
% [1] Lewin, L. 1958. Dilogarithms and associated functions. Macdonald.
% [2] Wood, D. C. 1992. Technical Report 15-92. University of Kent computing laboratory.
% [3] http://en.wikipedia.org/wiki/Polylog
% Didier Clamond, February 28th, 2006.
% Initialization.
d  = zeros(size(z));     
s  =  ones(size(z));     
      
% For large moduli: Mapping onto the unit circle |z|<=1.
j = find(abs(z)>1);
d(j) = -1.64493406684822643 - 0.5*log(-z(j)).^2; 
s(j) = -s(j);
z(j) = 1./z(j); 
% For large positive real parts: Mapping onto the unit circle with Re(z)<=1/2.
j = find(real(z)>0.5);
d(j) = d(j) + s(j).*( 1.64493406684822643 - log((1-z(j)).^log(z(j))) );
s(j) = -s(j);
z(j) = 1 - z(j);
% Transformation to Debye function and rational approximation.
z = -log(1-z);                                                                                
s = s.*z;                                                                                      
d = d - 0.25*s.*z;                                                                            
z = z.*z;                                                                                     
s = s.*(1+z.*(6.3710458848408100e-2+z.*(1.04089578261587314e-3+z*4.0481119635180974e-6)));    
s = s./(1+z.*(3.5932681070630322e-2+z.*(3.20543530653919745e-4+z*4.0131343133751755e-7)));    
d = d + s;                                                                                    

end