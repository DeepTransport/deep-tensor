
% setup the OU process
d = 3;
a = 0.7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = setup_ou_process(d, a);
func = @(x) eval_ou_process(data, x);

%%%%

debug_size = 1E4;
tmp = data.B\randn(d, debug_size);
tmp_f = feval(func, tmp);
deb = Debugger(tmp, tmp_f);
sample_x = data.B\randn(d, 1E3);

%%%%

% setup the reference polynomial

p = 50;
c1 = Chebyshev1st(p);
c2 = Chebyshev2nd(p);
le = Legendre(p);

bases{1} = ApproxBases(le, BoundedDomain([-4,4]), d);
bases{2} = ApproxBases(c1, LogarithmicMapping(4), d);

%

opt = SparseOption();
opt.tol = 1e-2;
opt.init_total_degree = 5;
%opt.indexset = 'hyperbolic';
opt.max_dim_basis = 2e4;
opt.max_sample_size = 2e5;
opt.enrich_degree = 1;
opt.init_sample_size = 2;
opt.enrich_sample_size = 1; 
opt.adaptation_rule = 'margin';


for i = 1:length(bases)
    tic;
    irts{i} = SparseSIRT(func, bases{i}, opt, 'debug', deb);
    toc
end

% sample
z = rand(d, 1E4);
for i = 1:length(bases)
    figure
    tic;[r,f] = eval_irt(irts{i}, z);toc
    tic;z0 = eval_rt(irts{i}, r);toc
    disp(' ')
    disp(['transform eror: ' num2str(norm(z-z0, 'fro'))])
    disp(['potential eror: ' num2str(norm(log(data.norm)+func(r) - f))])
    disp(['pdf eror: ' num2str(norm(exp(-log(data.norm)-func(r)) - exp(-f)))])
    disp(['cov eror: ' num2str(norm(data.C - cov(r'))/norm(data.C))])
    disp(' ')
    %
    subplot(2,2,1); plot(abs(log(data.norm)+func(r) - f), '.'); title('actual potential function vs fft')
    subplot(2,2,2); plot(log(data.norm)+func(r), f, '.');
    subplot(2,2,3); plot(data.C - cov(r')); title('actual covariance vs sample covariance')
end

