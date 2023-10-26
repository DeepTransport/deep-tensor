d = 5; a = 0.5;
model = OU(d, a);

nsteps = 1E4;
tic;
out = NUTS(@(x) eval_potential_nuts(model,x), randn(d,1), nsteps);
toc

norm(cov(out.samples') - C, 'fro')/norm(C, 'fro')
mean(out.samples')

% this is a bad example for pCN, only showing the code is running
tic;
out = pCN(@(x) eval_potential_pcn(model,x), randn(d,1), nsteps, 1);
toc

norm(cov(out.samples') - C, 'fro')/norm(C, 'fro')
mean(out.samples')

%%%%%

