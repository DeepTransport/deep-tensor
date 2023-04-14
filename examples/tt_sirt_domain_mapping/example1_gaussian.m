
% setup the OU process
d = 20;
a = 0.5;

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

doms{1} = BoundedDomain([-5,5]);
bases{1} = ApproxBases(Lagrangep(5,8), doms{1}, d);

doms{2} = AlgebraicMapping(4);
bases{2} = ApproxBases(Chebyshev1st(40), doms{2}, d);

doms{3} = AlgebraicMapping(4);
bases{3} = ApproxBases(Lagrangep(5,8), doms{3}, d);

doms{4} = LogarithmicMapping(4);
bases{4} = ApproxBases(Chebyshev2nd(40), doms{4}, d);

doms{5} = LogarithmicMapping(4);
bases{5} = ApproxBases(Legendre(40), doms{5}, d);

opt = FTTOption('tt_method', 'random', 'als_tol', 1E-4, 'local_tol', 1E-20, ...
    'init_rank', 20, 'max_rank', 20, 'kick_rank', 0, 'max_als', 2);

for i = 1:length(bases)
    tic;
    irts{i} = TTSIRT(func, bases{i}, opt, 'debug', deb, 'samples', sample_x);
    toc
end

for i = 1:length(bases)
    tic;
    irts{i} = marginalise(irts{i}, -1);
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

