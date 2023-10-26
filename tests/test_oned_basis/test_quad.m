function test_quad(polys, f, a, b)

ana = integral(f, a, b);
    
max_quad = size(polys, 2);

err = zeros(4, max_quad);
for nquad = 1:max_quad
    for j = 1:size(polys,1)
        nodes = polys{j,nquad}.nodes;
        omega = eval_measure(polys{j,nquad},nodes);
        weights = polys{j,nquad}.weights;
        err(j, nquad) = abs( sum(f(nodes)./omega.*weights) - ana );
    end
end

names = cell(1, size(polys,1));
for i = 1:size(polys,1)
    names{i} = class(polys{i,1});
end

figure
subplot(1,2,1)
xs = linspace(max(a,-5),min(b,5),1E3);
plot(xs, f(xs))
subplot(1,2,2)
semilogy(err', '-o')
legend(names)
title(func2str(f))

end