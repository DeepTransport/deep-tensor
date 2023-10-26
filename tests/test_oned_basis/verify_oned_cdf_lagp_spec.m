tau = 1E-2; 

pdf1 = @(x) exp(-8*x.^2);
pdf2 = @(x) exp(-5*abs(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

weight = 0.5;

lag = Lagrangep(6, 4);
lag_cdf = LagrangepCDF(lag);

fl = pdf2(lag.nodes)/weight;
pl = eval(lag, fl, lag_cdf.nodes);

xs = linspace(lag.domain(1), lag.domain(2), 1000);
fi = eval(lag, fl, xs);
Fi = eval_cdf(lag_cdf, pl+tau, xs);

z = rand(5E5,1);
tic;r = invert_cdf(lag_cdf, pl+tau, z);toc
tic; norm(eval_cdf(lag_cdf, pl+tau, r)-z), toc

dl = pdf2cdf(lag_cdf, pl);

figure
plot(xs, fi/dl.poly_norm, xs, Fi)
hold on
plot(lag_cdf.nodes, dl.cdf_poly_nodes/dl.poly_norm)
histogram(r, 'Normalization', 'pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lag = Lagrangep(6, 5);
lag_cdf = LagrangepCDF(lag);

fl = (pdf2(lag.nodes)/weight).^0.5;
pl = eval_radon(lag, fl, lag_cdf.nodes).^2;
dl = pdf2cdf(lag_cdf, pl);

xs = linspace(lag.domain(1), lag.domain(2), 10000);
fi = eval_radon(lag, fl, xs).^2*weight;
Fi = eval_cdf(lag_cdf, pl+tau, xs);

z = rand(5E5,1);
tic;r = invert_cdf(lag_cdf, pl+tau, z);toc
tic; norm(eval_cdf(lag_cdf, pl+tau, r)-z), toc


figure
plot(xs, fi/dl.poly_norm, xs, Fi)
hold on
plot(lag.nodes, pdf2(lag.nodes)/dl.poly_norm, 'o', lag_cdf.nodes, pl/dl.poly_norm, '-s')
histogram(r, 'Normalization', 'pdf')
set(gca, 'xlim', lag.domain)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


poly = Legendre(8);
poly_cdf = BoundedPolyCDF(poly);

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

poly = Legendre(20);
poly_cdf = BoundedPolyCDF(poly);

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

lag = Lagrangep(6, 4);
lag_cdf = LagrangepCDF(lag);

fl = (pdf2(lag.nodes)/weight).^(0.5);
pl = eval_radon(lag, fl, lag_cdf.nodes).^2;
%
pl = repmat(pl, 1, n);
dl = pdf2cdf(lag_cdf, pl);

xs = linspace(lag.domain(1), lag.domain(2), n);
fi = eval_radon(lag, fl, xs).^2*weight;
Fi = eval_cdf(lag_cdf, pl+tau, xs);

z = rand(n, 1);
tic;r = invert_cdf(lag_cdf, pl+tau, z);toc
tic; norm(eval_cdf(lag_cdf, pl+tau, r)-z), toc

figure
histogram(r, 'Normalization', 'pdf')
hold on
plot(xs, fi/dl.poly_norm(1), xs, Fi)
hold on
plot(lag.nodes, fl.^2/dl.poly_norm(1), 'o', lag_cdf.nodes, pl(:,1)/dl.poly_norm(1), 's')


%%%%%%

poly = Legendre(8);
poly_cdf = BoundedPolyCDF(poly);

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
