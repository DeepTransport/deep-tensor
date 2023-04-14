tau = 1E-2; 

pdf1 = @(x) exp(-8*(x-0.2).^2);
pdf2 = @(x) exp(-5*abs(x-0.2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

poly = Chebyshev2nd(20);
poly_cdf = Chebyshev2ndTrigoCDF(poly);

%poly = Chebyshev1st(20);
%poly_cdf = Chebyshev1stTrigoCDF(poly);

xs = linspace(-pi, 0, 10000);

[A,dA] = eval_int_basis_newton(poly_cdf, xs);
disp(norm((A(3:end,:)-A(1:end-2,:))/(xs(3)-xs(1)) - dA(2:end-1,:), 'fro')/norm(dA(2:end-1,:), 'fro'))

weight = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fp = poly.node2basis*( pdf2(poly.nodes)./eval_measure(poly, poly.nodes) );
pp = eval_radon(poly, fp, poly_cdf.nodes);

cdf_nodes = poly_cdf.cdf_basis2node*(poly_cdf.node2basis*pp);
cdf_nodes = cdf_nodes - cdf_nodes(1,:);
cdf_nodes = cdf_nodes./cdf_nodes(end,:);

xs = linspace(poly.domain(1), poly.domain(2), 10000);
fi = eval(poly, fp, xs);
Fi = eval_cdf(poly_cdf, pp+tau, xs);

figure
plot(xs, fi, xs, pdf2(xs))
hold on
plot(xs, Fi, poly_cdf.sampling_nodes, cdf_nodes, 'o')
set(gca, 'xlim', poly.domain)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fp = poly.node2basis*( (pdf2(poly.nodes)./eval_measure(poly, poly.nodes)).^0.5 );
pp = eval_radon(poly, fp, poly_cdf.nodes).^2;
%
cdf_nodes = poly_cdf.cdf_basis2node*(poly_cdf.node2basis*pp);
cdf_nodes = cdf_nodes - cdf_nodes(1,:);
cdf_nodes = cdf_nodes./cdf_nodes(end,:);

xs = linspace(poly.domain(1), poly.domain(2), 1000);
fi = eval_radon(poly, fp, xs).^2.*eval_measure(poly, xs(:));
Fi = eval_cdf(poly_cdf, pp+tau, xs);


z = rand(5E5,1);
tic;r = invert_cdf(poly_cdf, pp+tau, z);toc
tic; norm(eval_cdf(poly_cdf, pp+tau, r)-z), toc

figure
histogram(r, 'Normalization', 'pdf')
hold on
plot(xs, fi, xs, pdf2(xs))
plot(xs, Fi, poly_cdf.sampling_nodes, cdf_nodes, 'o')
set(gca, 'xlim', poly.domain)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% multi cdf test
n = 5E5;

%%%%%%

fp = poly.node2basis*( (pdf2(poly.nodes)./eval_measure(poly, poly.nodes)).^(0.5) );
pp = eval_radon(poly, fp, poly_cdf.nodes).^2;
pp = repmat(pp, 1, n);

cdf_nodes = poly_cdf.cdf_basis2node*(poly_cdf.node2basis*pp);
cdf_nodes = cdf_nodes - cdf_nodes(1,:);
cdf_nodes = cdf_nodes./cdf_nodes(end,:);

xs = linspace(poly.domain(1), poly.domain(2), n);
fi = eval_radon(poly, fp, xs).^2.*eval_measure(poly, xs(:));
Fi = eval_cdf(poly_cdf, pp+tau, xs);

z = rand(n,1);
taus = rand(1,n)*tau;
tic;r = invert_cdf(poly_cdf, pp+tau, z);toc
tic; norm(eval_cdf(poly_cdf, pp+tau, r)-z), toc


figure
histogram(r, 'Normalization', 'pdf')
hold on
plot(xs, fi, xs, pdf2(xs))
hold on
plot(xs, Fi, poly_cdf.sampling_nodes, cdf_nodes(:,1), 's')
set(gca, 'xlim', poly.domain)
