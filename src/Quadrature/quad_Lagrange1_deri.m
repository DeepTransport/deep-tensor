function [f,g] = quad_Lagrange1_deri(poly, theta, x)
% compute the integration of the derivative of piecewise Lagrange 1
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
% For lagrange 1: determine the function a + bx, then the derivative is
% piecewise constant with value b
%

if poly.order ~= 1
    error('polynomial order should be 1')
end

% index in which the boundary x falls into
% ei = 0: left ghost cell, ei = poly.num_elems+1: right ghost cell
ei = ceil( (x-poly.grid(1))./poly.elem_size );
%
f = zeros(size(x));
g = zeros(size(theta));
for k = 1:poly.num_elems
    ind1 = (ei>k);  % integrate the whole element
    ind2 = (ei==k); % integrate to x
    n1 = sum(ind1);
    n2 = sum(ind2);
    %
    %b = (theta(k+1,:)-theta(k,:))./poly.elem_size;
    %tmp = 2.^b + 1;
    tmp = 2.^((theta(k+1,:)-theta(k,:))./poly.elem_size);
    if n1 > 0
        f(ind1) = f(ind1) + log2(tmp(ind1)+1)*poly.elem_size;
        g(k,ind1) = g(k,ind1) - tmp(ind1)./(tmp(ind1)+1);
        g(k+1,ind1) = g(k+1,ind1) + tmp(ind1)./(tmp(ind1)+1);
    end
    if n2 > 0
        dx = x(ind2) - poly.grid(k);
        f(ind2) = f(ind2) + log2(tmp(ind2)+1).*dx;
        g(k,ind2) = g(k,ind2) - tmp(ind2)./(tmp(ind2)+1).*(dx/poly.elem_size);
        g(k+1,ind2) = g(k+1,ind2) + tmp(ind2)./(tmp(ind2)+1).*(dx/poly.elem_size);
    end
end

end

%{
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
        %b = (theta(k+1,:)-theta(k,:))./poly.elem_size;
        %tmp = 2.^b + 1;
        tmp = 2.^((theta(k+1,:)-theta(k,:))./poly.elem_size);
        if n1 > 0
            f(ind1) = f(ind1) + log2(tmp(ind1)+1)*poly.elem_size;
            g(k,ind1) = g(k,ind1) - tmp(ind1)./(tmp(ind1)+1);
            g(k+1,ind1) = g(k+1,ind1) + tmp(ind1)./(tmp(ind1)+1);
        end
        if n2 > 0
            dx = x(ind2) - poly.grid(k);
            f(ind2) = f(ind2) + log2(tmp(ind2)+1).*dx;
            g(k,ind2) = g(k,ind2) - tmp(ind2)./(tmp(ind2)+1).*(dx/poly.elem_size);
            g(k+1,ind2) = g(k+1,ind2) + tmp(ind2)./(tmp(ind2)+1).*(dx/poly.elem_size);
        end
    end
end
%}
