function [f,grad] = quad_piecewise_deri(poly, theta, x)
% compute the integration of the derivative of piecewise Lagrange 1 and 
% Lagrange 2 polynomials from poly.domain(1) to x
%
%   theta - interpolation function values, d x n, where n is the number of 
%           functions and d is the number of interpolation 
%   x     - 1 x n, number of terminal points


% the function g(y) = log(2.^y + 1)/log(2);
%
% size of each element is given by poly.elem_size, except ghost cells
%
% the derivative of ghost cells are non-zero if 'Dirichlet' bc is used
% the coefficient is given by the following
gh_left = 0;
gh_right = 0;
if poly.gs > 0 && strcmp(poly.bc, 'Dirichlet')
    gh_left = theta(1,:)./poly.gs;
    gh_right = -theta(end,:)./poly.gs;
end

% For lagrange 1: determine the function a + bx, then the derivative is
% piecewise constant with value b
%
% For lagrange 2: determine the function a+bx+cx^2, then the derivative is
% piecewise linear with value b+2cx
%
if poly.order == 1
    a = (theta(2:end,:) - theta(1:end-1,:))./poly.elem_size;
elseif poly.order == 2
    %
    dhh = poly.elem_size/2;
    %
    ii = zeros(3,poly.num_elems);
    jj = zeros(3,poly.num_elems);
    %
    % L = [1, 0, 0; 1, 1, 0; 1, 2, 1];
    % U = [1, 0, 0; 0, dhh, dhh^2; 0, 0, dhh^2*2];
    % LU = Vandermore
    % Vandermore = [1, 0, 0; 1, dhh, dhh^2; 1, dhh*2, 4*dhh^2];
    % inverse of the Vandermore matrix
    iV = [1, 0, 0; -3/(2*dhh), 2/dhh, -1/(2*dhh); 1/(2*dhh^2), -1/dhh^2, 1/(2*dhh^2)];
    for i = 1:obj.num_elems
        ind = (1:3)+(i-1)*3;
        ii(:,i) = ind;
        jj(:,i) = (i-1)*2 + (1:3);
    end
    T = sparse(ii(:), jj(:), ones(poly.num_elems*3,1), poly.num_elems*3, poly.num_nodes);
    coef = iV*reshape(T*theta, 3, []); % a b c for each elements
    % we need to integrate each local element from 0
else
    error('polynomial order more than 2, not supported')
end

% index in which the boundary x falls into
% ei = 0: left ghost cell, ei = poly.num_elems+1: right ghost cell
ei = ceil( (x-poly.grid(1))./poly.elem_size );

for k = 0:(poly.num_elems+1)
    ind1 = (ei>k);  % integrate the whole element
    ind2 = (ei==k); % integrate to x
    n1 = sum(ind1);
    n2 = sum(ind2);
    %
    if k == 0
        
    elseif k == (poly.num_elems+1)
    else
    end
end



end

