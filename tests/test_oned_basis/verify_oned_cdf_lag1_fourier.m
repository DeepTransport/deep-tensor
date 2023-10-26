tau = 1E-2; 

pdf1 = @(x) exp(-8*x.^2);
pdf2 = @(x) exp(-5*abs(x));

weight = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

lag = Lagrange1(10);
lag_cdf = Lagrange1CDF(lag);

fl = reshape(pdf2(lag.nodes),[],1)/weight;
pl = eval_radon(lag, fl, lag_cdf.nodes);

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
plot(lag_cdf.grid, dl.cdf_poly_grid/dl.poly_norm, 'o')
ksdensity(r)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lag = Lagrange1(40);
lag_cdf = Lagrange1CDF(lag);

fl = reshape((pdf2(lag.nodes)/weight).^0.5,[],1);
pl = eval_radon(lag, fl, lag_cdf.nodes).^2;

xs = linspace(lag.domain(1), lag.domain(2), 10000);
fi = eval_radon(lag, fl, xs).^2*weight;
Fi = eval_cdf(lag_cdf, pl+tau, xs);

z = rand(5E5,1);
tic;r = invert_cdf(lag_cdf, pl+tau, z);toc
tic; norm(eval_cdf(lag_cdf, pl+tau, r)-z), toc

dl = pdf2cdf(lag_cdf, pl);

figure
plot(xs, fi/dl.poly_norm, xs, Fi)
hold on
plot(lag.nodes, fl.^2/dl.poly_norm, 'o', lag_cdf.nodes, pl/dl.poly_norm, 's')
ksdensity(r)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% multi cdf test
n = 5E5;

lag = Lagrange1(40);
lag_cdf = Lagrange1CDF(lag);

fl = reshape((pdf2(lag.nodes)/weight).^(0.5),[],1);
pl = eval_radon(lag, fl, lag_cdf.nodes).^2;
%
pl = repmat(pl, 1, n);
dl = pdf2cdf(lag_cdf, pl);

xs = linspace(lag.domain(1), lag.domain(2), n);
fi = eval_radon(lag, fl, xs).^2*weight;
Fi = eval_cdf(lag_cdf, pl+tau, xs);

z = rand(n,1);
taus = rand(1,n)*tau;
tic;r = invert_cdf(lag_cdf, pl+tau, z);toc
tic; norm(eval_cdf(lag_cdf, pl+tau, r)-z), toc

figure
plot(xs, fi/dl.poly_norm(1), xs, Fi)
hold on
plot(lag.nodes, fl.^2/dl.poly_norm(1), 'o', lag_cdf.nodes, pl(:,1)/dl.poly_norm(1), 's')
ksdensity(r)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
poly = Fourier(8);
poly_cdf = FourierCDF(poly);

fp = poly.node2basis*pdf2(poly.nodes(:))/weight;
pp = eval_radon(poly, fp, poly_cdf.nodes(:));

cdf_nodes = poly_cdf.cdf_basis2node*(poly_cdf.node2basis*pp);
cdf_nodes = cdf_nodes - cdf_nodes(1,:);
cdf_nodes = cdf_nodes./cdf_nodes(end,:);

xs = linspace(poly.domain(1), poly.domain(2), 10000);
fi = eval(poly, fp, xs);
Fi = eval_cdf(poly_cdf, pp+tau, xs);

figure
plot(xs, fi, xs, Fi)
hold on
plot(poly_cdf.sampling_nodes, cdf_nodes, 'o')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poly = Fourier(8);
poly_cdf = FourierCDF(poly);

fp = poly.node2basis*( (pdf2(poly.nodes(:))/weight).^(0.5) );
pp = eval_radon(poly, fp, poly_cdf.nodes(:)).^2;
%
cdf_nodes = poly_cdf.cdf_basis2node*(poly_cdf.node2basis*pp);
cdf_nodes = cdf_nodes - cdf_nodes(1,:);
cdf_nodes = cdf_nodes./cdf_nodes(end,:);

xs = linspace(poly.domain(1), poly.domain(2), 1000);
fi = eval_radon(poly, fp, xs).^2*weight;
Fi = eval_cdf(poly_cdf, pp+tau, xs);


z = rand(5E5,1);
tic;r = invert_cdf(poly_cdf, pp+tau, z);toc
tic; norm(eval_cdf(poly_cdf, pp+tau, r)-z), toc

figure
plot(xs, fi, xs, Fi)
hold on
plot(poly_cdf.sampling_nodes, cdf_nodes, 'o')
ksdensity(r)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 5E5;

poly = Fourier(8);
poly_cdf = FourierCDF(poly);

fp = poly.node2basis*( (pdf2(poly.nodes(:))/weight).^(0.5) );
pp = eval_radon(poly, fp, poly_cdf.nodes(:)).^2;
pp = repmat(pp, 1, n);

cdf_nodes = poly_cdf.cdf_basis2node*(poly_cdf.node2basis*pp);
cdf_nodes = cdf_nodes - cdf_nodes(1,:);
cdf_nodes = cdf_nodes./cdf_nodes(end,:);

xs = linspace(poly.domain(1), poly.domain(2), n);
fi = eval(poly, fp, xs).^2*weight;
Fi = eval_cdf(poly_cdf, pp+tau, xs);

z = rand(n,1);
taus = rand(1,n)*tau;
tic;r = invert_cdf(poly_cdf, pp+tau, z);toc
tic; norm(eval_cdf(poly_cdf, pp+tau, r)-z), toc

figure
plot(xs, fi.^2, xs, Fi)
hold on
plot(poly_cdf.nodes, pp(:,1), 'o', poly_cdf.sampling_nodes, cdf_nodes(:,1), 's')
ksdensity(r)

