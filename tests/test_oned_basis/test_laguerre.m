
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

j = 5;
A = polys{j}.basis2node;
disp(A'*diag(polys{j}.weights)*A)

j = 20;
A = polys{j}.basis2node;
figure
plot(A'*diag(polys{j}.weights)*A)
disp(mean(abs(diag(A'*diag(polys{j}.weights)*A))))


test_proj_vs_matlab(polys{10}, f , 0, 10, 0, 6 );

test_proj_vs_matlab(polys{5}, @(x) 1./(1+x.^2) , 0, 10, 0, 6 );

