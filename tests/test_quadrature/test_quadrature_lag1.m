

% piecewise polynomials
% coefficients of piecewise Lagrange polynomials (nodal values) natually
% have an innner product weighted by poly.mass_R' * poly.mass_R, where
% poly.mass_R is the upper triangular Choleskt factor of the mass matrix
%
poly = Lagrange1(5, [-4, 4]);

npts = 200;

% integrate g(f'(x, theta)) from domain(1) to x
% generate x as a row vector
x = rand(1,npts)*4;

theta = randn(poly.num_nodes, npts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
[y1,grad] = quad_Lagrange1_deri(poly, theta, x);
toc;

%
%Simpson rule by matlab
y2 = y1;
for i = 1:length(x)
    y2(i) = quad(@(z) g(poly, theta(:,i), z) , poly.domain(1), x(i), 1E-16 );
end
norm(y2 - y1)/norm(y2)

[quadx,quadw] = quadcc_piecewise(poly, 17, x);
gx = g_block(poly, theta, quadx);
y3 = sum(gx.*quadw,1);
norm(y3 - y1)/norm(y3)
norm(y3 - y2)/norm(y3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tol = 1E-5;
y1p = grad;
y1m = grad;
for j = 1:poly.num_nodes
    thetap = theta;
    thetap(j,:) = thetap(j,:) + tol;
    thetam = theta;
    thetam(j,:) = thetam(j,:) - tol;
    y1p(j,:) = quad_Lagrange1_deri(poly, thetap, x);
    y1m(j,:) = quad_Lagrange1_deri(poly, thetam, x);
end
gd = (y1p-y1m)/(2*tol);
plot(gd(:), grad(:), '.')
norm(gd-grad,'fro')
 
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
