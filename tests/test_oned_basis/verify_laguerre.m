
max_order = 20;
polys = cell(1, max_order);
for order = 1:max_order
    polys{order} = Laguerre(order);
end

f = @(y) exp(-6* y.^2 -y.^4);
test_quad(polys, f, 0, inf);

f = @(y) exp(-y.^2 - 0.5*y.^4 - 0.5*y.^2);
test_quad(polys, f, 0, inf);

f = @(y) log(y.^2+1) .* exp( - 0.5*y.^2);
test_quad(polys, f, 0, inf);

err = 0;
for j = 1:max_order
    A = polys{j}.basis2node;
    err = err + norm(A'*diag(polys{j}.weights)*A - eye(size(A)), 'fro');
end
disp(err)


test_proj_vs_matlab(polys{10}, f , 0, 10, 0, 6 );

test_proj_vs_matlab(polys{5}, @(x) 1./(1+x.^2) , 0, 10, 0, 6 );

