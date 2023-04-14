

% define some polynomials in the test
% it does not be a good idea to integrate Langrage polynomials using CC
%
% coefficients of spectral polynomials natually have a L2 innner product
%

polys{1} = Legendre(10, [-4, 4]);
polys{2} = Fourier(5, [-4, 4]);

% piecewise polynomials 
% coefficients of piecewise Lagrange polynomials (nodal values) natually 
% have an innner product weighted by poly.mass_R' * poly.mass_R, where 
% poly.mass_R is the upper triangular Choleskt factor of the mass matrix
%
polys{3} = Lagrange1(10, [-4, 4]);
polys{4} = Lagrangep(2, 5, [-4, 4]);

% integrate g(f'(x, theta)) from domain(1) to x
% generate x as a row vector
x = rand(1,100)*4;

theta = randn(polys{1}.num_nodes, 1);
[y1, n1] = quadcc(@(z) g(polys{1}, theta, z) , polys{1}.domain(1)*ones(size(x)), x, 1E-5);

theta = randn(polys{2}.num_nodes, 1);
[y2, n2] = quadcc(@(z) g(polys{2}, theta, z) , polys{2}.domain(1)*ones(size(x)), x, 1E-5);

theta = randn(polys{3}.num_nodes, 1);
%directly applying CC does not work
%[y3, n3] = quadcc(@(z) g(polys{3}, theta, z) , polys{3}.domain(1)*ones(size(x)), x);
%use the piecewise implementation
[y3, n3] = quad_piecewise(polys{3}, theta, x);
%
%Simpson rule by matlab
y3d = y3;
n3d = 0;
for i = 1:length(x)
    [y3d(i),n] = quad(@(z) g(polys{3}, theta, z) , polys{3}.domain(1), x(i), 1E-5 );
    n3d = n3d + n;
end
norm(y3d - y3)/norm(y3d)

theta = randn(polys{4}.num_nodes, 1);
[y4, n4] = quad_piecewise(polys{4}, theta, x);
y4d = y4;
n4d = 0;
for i = 1:length(x)
    [y4d(i),n] = quad(@(z) g(polys{4}, theta, z) , polys{4}.domain(1), x(i), 1E-5 );
    n4d = n4d + n;
end
norm(y4d - y4)/norm(y4d)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y, fcnt] = quad_piecewise(poly, theta, x)

% only use Neumann boundary condition for the ghost cell here

% index in which the boundary x falls into
% ei = 0: left ghost cell, ei = poly.num_elems+1: right ghost cell
ei = ceil( (x-poly.grid(1))./poly.elem_size );

y = zeros(size(x));
fcnt = 0;
for k = 0:(poly.num_elems+1)
    ind1 = (ei>k);  % integrate the whole element
    ind2 = (ei==k); % integrate to x
    n1 = sum(ind1);
    n2 = sum(ind2);
    c1 = 0;
    c2 = 0;
    if k == 0
        if n1 > 0
            [y1,c1] = quadcc(@(z) g(poly, theta, z), poly.domain(1), poly.grid(1), 1E-5);
            y(ind1) = y(ind1) + y1;
        end
        if n2 > 0
            [y2,c2] = quadcc(@(z) g(poly, theta, z), poly.domain(1)*ones(1,n2), x(ind2), 1E-5);
            y(ind2) = y(ind2) + y2;
        end
    else
        if n1 > 0
            [y1,c1] = quadcc(@(z) g(poly, theta, z), poly.grid(k), poly.grid(k+1), 1E-5);
            y(ind1) = y(ind1) + y1;
        end
        if n2 > 0
            [y2,c2] = quadcc(@(z) g(poly, theta, z), poly.grid(k)*ones(1,n2), x(ind2), 1E-5);
            y(ind2) = y(ind2) + y2;
        end
    end
    fcnt = fcnt + c1 + c2;
end

end

function y = g(poly, theta, z)

b = eval_basis_deri(poly, z(:));
y = reshape(b*theta, size(z));
y = log(2.^y + 1)/log(2);

end