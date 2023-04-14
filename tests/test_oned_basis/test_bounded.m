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
j = 5;
for i = 1:4
    disp(i)
    A = polys{i,j}.basis2node; 
    disp(A'*diag(polys{i,j}.weights)*A)
end

j = 20;
for i = 1:4
    disp(i)
    A = polys{i,j}.basis2node; 
    figure
    plot(A'*diag(polys{i,j}.weights)*A)
    disp(mean(abs(diag(A'*diag(polys{i,j}.weights)*A))))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% projection test
j = 10;
for i = 1:4
    test_proj_vs_matlab(polys{i,j}, @(y) log(y+1.1) , -1, 1, -1, 1 );
end
