tau = 1E-1; 

pdf = @(y) log(y.^2+1) .* exp( - 0.5*y.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poly = Laguerre(20);
poly_cdf = LaguerreCDF(poly);

xs = linspace(0, 10, 10000);

[A,dA] = eval_int_basis_newton(poly_cdf, xs);
disp(norm((A(3:end,:)-A(1:end-2,:))/(xs(3)-xs(1)) - dA(2:end-1,:), 'fro')/norm(dA(2:end-1,:), 'fro'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poly = Laguerre(8);
poly_cdf = LaguerreCDF(poly);

fp = poly.node2basis*(pdf(poly.nodes)./eval_measure(poly, poly.nodes));
pp = eval_radon(poly, fp, poly_cdf.nodes);

cdf_nodes = poly_cdf.cdf_basis2node*(poly_cdf.node2basis*pp);
cdf_nodes = cdf_nodes - cdf_nodes(1,:);
cdf_nodes = cdf_nodes./cdf_nodes(end,:);

xs = linspace(0, 10, 1E3);
fi = eval(poly, fp, xs);
Fi = eval_cdf(poly_cdf, pp+tau, xs);

figure
plot(xs, fi, xs, pdf(xs))
hold on
plot(xs, Fi, poly_cdf.sampling_nodes, cdf_nodes, 'o')
set(gca, 'xlim', [0, 10])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poly = Laguerre(20);
poly_cdf = LaguerreCDF(poly);

fp = poly.node2basis*( (pdf(poly.nodes)./eval_measure(poly, poly.nodes)).^(0.5) );
pp = eval_radon(poly, fp, poly_cdf.nodes).^2;
%
cdf_nodes = poly_cdf.cdf_basis2node*(poly_cdf.node2basis*pp);
cdf_nodes = cdf_nodes - cdf_nodes(1,:);
cdf_nodes = cdf_nodes./cdf_nodes(end,:);

xs = linspace(0, 10, 1E3);
fi = eval_radon(poly, fp, xs).^2.*eval_measure(poly, xs(:));
Fi = eval_cdf(poly_cdf, pp+tau, xs);

z = rand(5E5,1);
tic;r = invert_cdf(poly_cdf, pp+tau, z);toc
tic; norm(eval_cdf(poly_cdf, pp+tau, r)-z), toc

figure
plot(xs, fi, xs, pdf(xs))
hold on
plot(xs, Fi, poly_cdf.sampling_nodes, cdf_nodes, 'o')
ksdensity(r)
set(gca, 'xlim', [0, 10])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% multi cdf test
n = 5E5;

poly = Laguerre(20);
poly_cdf = LaguerreCDF(poly);

fp = poly.node2basis*( (pdf(poly.nodes)./eval_measure(poly, poly.nodes)).^(0.5) );
pp = eval_radon(poly, fp, poly_cdf.nodes).^2;
pp = repmat(pp, 1, n);

cdf_nodes = poly_cdf.cdf_basis2node*(poly_cdf.node2basis*pp);
cdf_nodes = cdf_nodes - cdf_nodes(1,:);
cdf_nodes = cdf_nodes./cdf_nodes(end,:);

xs = linspace(0, 10, n);
fi = eval_radon(poly, fp, xs).^2.*eval_measure(poly, xs(:));
Fi = eval_cdf(poly_cdf, pp+tau, xs);

z = rand(n,1);
taus = rand(1,n)*tau;
tic;r = invert_cdf(poly_cdf, pp+tau, z);toc
tic; norm(eval_cdf(poly_cdf, pp+tau, r)-z), toc


figure
plot(xs, fi, xs, pdf(xs))
hold on
plot(xs, Fi, poly_cdf.sampling_nodes, cdf_nodes(:,1), 'o')
ksdensity(r)
set(gca, 'xlim', [0, 10])