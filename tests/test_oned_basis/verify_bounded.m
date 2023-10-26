% test the quadrature
max_order = 20;
polys = cell(4, max_order);
for order = 1:max_order
    polys{1, order} = Chebyshev1st(order);
    polys{2, order} = Chebyshev2nd(order);
    polys{3, order} = Legendre(order);
    polys{4, order} = Jacobi11(order);
end



%%%%%%%% case 1
fs = {@(y) sin(y*2*pi) + y.^2, @(y) sin(y*2*pi) + log(y+2), ...
    @(y) sin(y*1.5*pi) + log(y+2), @(y) sqrt(y+1), ...
    @(y) 1./sqrt(y+1), @(y) 1./(2+y.^2)};

for k = 1:length(fs)
    test_quad(polys, fs{k}, -1, 1);
end

% discrete orthogonality 
for i = 1:4
    disp(i)
    err = 0;
    for j = 1:max_order
        A = polys{i,j}.basis2node; 
        err = err + norm(A'*diag(polys{i,j}.weights)*A - eye(size(A)), 'fro');
    end
    disp(err)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% projection test
j = 10;
for i = 1:4
    test_proj_vs_matlab(polys{i,j}, @(y) log(y+1.1) , -1, 1, -1, 1 );
end
