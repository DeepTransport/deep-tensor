
% setup the OU process
d = 20;
a = 0.5;
%p = randperm(d,d);
%pt(p) = 1:length(p);

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
%{
doms{1} = BoundedDomain([-5,5]);
bases{1} = ApproxBases(Lagrangep(5,8), doms{1}, d);

doms{2} = AlgebraicMapping(4);
bases{2} = ApproxBases(Chebyshev1st(40), doms{2}, d);

doms{3} = AlgebraicMapping(4);
bases{3} = ApproxBases(Lagrangep(5,8), doms{3}, d);

doms{4} = LogarithmicMapping(4);
bases{4} = ApproxBases(Chebyshev1st(40), doms{4}, d);

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
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% should test ind = 1, ind = 1:(d-1) for > 0
% should test ind = d, ind = 2:d for < 0
% sample
z = rand(d, 1E4);
for i = 1:length(bases)
    figure;
    % test 1
    ind  = 1:8;
    if irts{i}.int_dir ~= 1, irts{i} = marginalise(irts{i}, 1); end
    tic;[r,p] = eval_irt(irts{i}, z(ind,:));toc
    fx = eval_pdf(irts{i}, r);
    tic;z0 = eval_rt(irts{i}, r);toc
    disp(' ')
    disp(['transform eror: ' num2str(norm(z(ind,:) - z0))])
    disp(['density eror: ' num2str(norm(exp(-p) - fx))])
    fe   = eval_ou_process_marginal(data, ind, r);
    disp(['approx eror: ' num2str(norm(exp(-p) - exp(-fe)))])
    disp(' ')
    %
    subplot(2,3,1);plot(abs(fe - p)/max(abs(fe)), '.');
    subplot(2,3,2);plot(fe , p, '.');
    title('actual poential function value vs fft')
    subplot(2,3,3);plot(data.C(ind, ind) - cov(r'))
    
    % test 2
    ind  = d:-1:15;
    if irts{i}.int_dir ~= -1, irts{i} = marginalise(irts{i}, -1); end
    tic;[r,p] = eval_irt(irts{i}, z(ind,:));toc
    fx = eval_pdf(irts{i}, r);
    tic;z0 = eval_rt(irts{i}, r);toc
    disp(' ')
    disp(['transform eror: ' num2str(norm(z(ind,:) - z0))])
    disp(['density eror: ' num2str(norm(exp(-p) - fx))])
    %
    fe   = eval_ou_process_marginal(data, ind, r);
    disp(['approx eror: ' num2str(norm(exp(-p) - exp(-fe)))])
    disp(' ')
    subplot(2,3,4);plot(abs(fe - p)/max(abs(fe)), '.');
    subplot(2,3,5);plot(fe , p, '.');
    title('actual potential function value vs fft')
    subplot(2,3,6);plot(data.C(ind, ind) - cov(r'))
end



