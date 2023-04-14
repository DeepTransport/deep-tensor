

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

npts = 100;

% integrate g(f'(x, theta)) from domain(1) to x
% generate x as a row vector
x = rand(1,npts)*4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = randn(polys{1}.num_nodes, npts);
[y1, n1, z, w, dxdz, mid] = quadccp(@(theta,z) g_block(polys{1},theta,z), theta, polys{1}.domain(1)*ones(size(x)), x, 1E-5);

% fixed quad rule
nquad = 100;
log_order = ceil(log2(nquad-1));
[z,w] = cc_rule(log_order);
dxdz = (x - polys{1}.domain(1)*ones(size(x)))/2;
mid = (x + polys{1}.domain(1)*ones(size(x)))/2;
quadx = z.*dxdz + mid;
quadw = w.*dxdz;
gx = g_block(polys{1}, theta, quadx);
y1f = sum(gx.*quadw,1);

y1d = y1;
n1d = 0;
for i = 1:length(x)
    [y1d(i),n] = quad(@(z) g(polys{1}, theta(:,i), z) , polys{1}.domain(1), x(i), 1E-10 );
    n1d = n1d + n;
end
norm(y1d - y1)/norm(y1d)
norm(y1d - y1f)/norm(y1d)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = randn(polys{2}.num_nodes, npts);
[y2, n2] = quadccp(@(theta,z) g_block(polys{2},theta,z), theta, polys{2}.domain(1)*ones(size(x)), x, 1E-5);

% fixed quad rule
nquad = 100;
log_order = ceil(log2(nquad-1));
[z,w] = cc_rule(log_order);
dxdz = (x - polys{2}.domain(1)*ones(size(x)))/2;
mid = (x + polys{2}.domain(1)*ones(size(x)))/2;
quadx = z.*dxdz + mid;
quadw = w.*dxdz;
gx = g_block(polys{2}, theta, quadx);
y2f = sum(gx.*quadw,1);


y2d = y2;
n2d = 0;
for i = 1:length(x)
    [y2d(i),n] = quad(@(z) g(polys{2}, theta(:,i), z) , polys{2}.domain(1), x(i), 1E-10 );
    n2d = n2d + n;
end
norm(y2d - y2)/norm(y2d)
norm(y2d - y2f)/norm(y2d)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = randn(polys{3}.num_nodes, npts);
%directly applying CC does not work
%[y3, n3] = quadcc(@(z) g(polys{3}, theta, z) , polys{3}.domain(1)*ones(size(x)), x);
%use the piecewise implementation
%[y3, n3] = quad_piecewise(polys{3}, theta, x);

[quadx,quadw] = quadcc_piecewise(polys{3}, 9, x);
gx = g_block(polys{3}, theta, quadx);
y3 = sum(gx.*quadw,1);

%
%Simpson rule by matlab
y3d = y3;
n3d = 0;
for i = 1:length(x)
    [y3d(i),n] = quad(@(z) g(polys{3}, theta(:,i), z) , polys{3}.domain(1), x(i), 1E-10 );
    n3d = n3d + n;
end
norm(y3d - y3)/norm(y3d)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = randn(polys{4}.num_nodes, npts);
[quadx,quadw] = quadcc_piecewise(polys{4}, 15, x);
gx = g_block(polys{4}, theta, quadx);
y4 = sum(gx.*quadw,1);

y4d = y4;
n4d = 0;
for i = 1:length(x)
    [y4d(i),n] = quad(@(z) g(polys{4}, theta(:,i), z) , polys{4}.domain(1), x(i), 1E-10 );
    n4d = n4d + n;
end
norm(y4d - y4)/norm(y4d)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = g_block(poly, theta, z)
% z: m x n
% theta: dof x n

[m, n] = size(z);
dof = size(theta,1);

b = eval_basis_deri(poly, z(:)); % (mn) x dof
b = reshape(permute(reshape(full(b),m,n,dof), [1,3,2]), m, []);
B = reshape(permute(reshape(b.*theta(:)',m,dof,n),[1,3,2]),m*n,[]);
y = reshape(sum(B,2),m,n);
y = log(2.^y + 1)/log(2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = g(poly, theta, z)

b = eval_basis_deri(poly, z(:));
y = reshape(b*theta, size(z));
y = log(2.^y + 1)/log(2);

end
