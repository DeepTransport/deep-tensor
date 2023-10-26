tau = 1E-2; 

pdf1 = @(x) exp(-8*x.^2);
pdf2 = @(x) exp(-5*abs(x));

weight = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

lag = Lagrange1(10);
lag_cdf = Lagrange1AbsCDF(lag);

fl = reshape(pdf2(lag.nodes),[],1)/weight;
pl = eval_radon(lag, fl, lag_cdf.nodes);
dl = pdf2cdf(lag_cdf, pl);

xs = linspace(lag.domain(1), lag.domain(2), 1000);
fi = eval_pdf_radon(lag_cdf, fl, xs)*weight;
Fi = eval_cdf(lag_cdf, pl+tau, xs);

z = rand(5E5,1);
tic;[r,p] = invert_cdf(lag_cdf, pl+tau, z);toc
%
tic; norm(eval_cdf(lag_cdf, pl+tau, r)-z), toc
%
norm(eval_pdf_radon(lag_cdf, fl+tau, r) - p)

figure
plot(xs, fi, xs, Fi)
hold on
plot(lag_cdf.grid, [0;dl.cdf_grid_unnormalised_right]/dl.norm, 'o')
ksdensity(r)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lag = Lagrange1(40);
lag_cdf = Lagrange1AbsCDF(lag);

fl = reshape(pdf2(lag.nodes),[],1)/weight;
pl = eval_radon(lag, fl, lag_cdf.nodes);
dl = pdf2cdf(lag_cdf, pl);

xs = linspace(lag.domain(1), lag.domain(2), 1000);
fi = eval_pdf_radon(lag_cdf, pl, xs)*weight;
Fi = eval_cdf(lag_cdf, pl+tau, xs);

z = rand(5E5,1);
tic;[r,p] = invert_cdf(lag_cdf, pl+tau, z);toc
%
tic; norm(eval_cdf(lag_cdf, pl+tau, r)-z), toc
%
norm(eval_pdf_radon(lag_cdf, fl+tau, r) - p)


figure
plot(xs, fi, xs, Fi)
hold on
plot(lag_cdf.grid, [0;dl.cdf_grid_unnormalised_right]/dl.norm, 'o')
ksdensity(r)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
n = 5E5;

lag = Lagrange1(40);
lag_cdf = Lagrange1AbsCDF(lag);

fl = reshape(pdf2(lag.nodes),[],1)/weight;
pl = eval_radon(lag, fl, lag_cdf.nodes);
pl = repmat(pl, 1, n);
dl = pdf2cdf(lag_cdf, pl);

xs = linspace(lag.domain(1), lag.domain(2), n);
fi = eval_pdf_radon(lag_cdf, pl+tau, xs)*weight;
Fi = eval_cdf(lag_cdf, pl+tau, xs);

z = rand(n,1);
tic;r = invert_cdf(lag_cdf, pl+tau, z);toc
tic; norm(eval_cdf(lag_cdf, pl+tau, r)-z), toc

figure
plot(xs, fi, xs, Fi)
hold on
plot(lag_cdf.nodes, [0;dl.cdf_grid_unnormalised_right(:,1)]/dl.norm(1), 's')
ksdensity(r)


