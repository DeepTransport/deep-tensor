% one dim gaussian

fun = @(x) 0.5*x.^2;
g_fun = @(x) x;

scales = 1:9;
for i = 1:numel(scales)
    map_A{i} = AlgebraicMapping(scales(i));
    map_L{i} = LogarithmicMapping(scales(i));
end
c1 = Chebyshev1st(20);
c2 = Chebyshev2nd(20);
le = Legendre(20);

p = 2;

oned_test(fun, g_fun, c1, map_A, p)

%oned_test(fun, g_fun, c1, map_A, p)
%oned_test(fun, g_fun, c2, map_A, p)
%oned_test(fun, g_fun, le, map_A, p)
