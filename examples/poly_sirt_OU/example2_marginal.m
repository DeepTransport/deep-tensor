
% setup the OU process
d = 4;
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

p = 50;
c1 = Chebyshev1st(p);
c2 = Chebyshev2nd(p);
le = Legendre(p);

bases{1} = ApproxBases(le, BoundedDomain([-4,4]), d);
bases{2} = ApproxBases(c1, LogarithmicMapping(4), d);

%
%{
opt = SparseOption();
opt.tol = 1e-2;
opt.init_total_degree = 3;
opt.max_dim_basis = 2e3;
opt.max_sample_size = 2e4;
opt.enrich_degree = 1;
opt.init_sample_size = 3;
opt.enrich_sample_size = 1;
opt.adaptation_rule = 'reducedmargin';


for i = 1:length(bases)
    tic;
    irts{i} = SparseSIRT(func, bases{i}, opt, 'debug', deb);
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
    ind  = 1:1;
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
    ind  = d:-1:2;
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



