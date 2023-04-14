
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

dom = BoundedDomain([-5,5]);
bases{1} = ApproxBases(Legendre(40), dom, d);
bases{2} = ApproxBases(Fourier(20), dom, d);
bases{3} = ApproxBases(Lagrange1(40), dom, d);
bases{4} = ApproxBases(Lagrangep(5,8), dom, d);
%bases{5} = ApproxBases(Hermite(10), UnboundedDomain(), d);

opts{1} = FTTOption('tt_method', 'amen', ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 4);
opts{2} = FTTOption('tt_method', 'random', ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 4);

for i = 1:length(bases)
    for j = 1:2
        tic;
        irts{i,j} = TTSIRT(func, bases{i}, opts{j}, 'debug', deb, 'samples', sample_x);
        toc
    end
    tmp = round(irts{i,1}.approx, 1E-2);
    irts{i,3} = TTSIRT(func, tmp, 'debug', deb, 'samples', sample_x);
end

for i = 1:length(bases)
    for j = 1:3
        tic;
        irts{i,j} = marginalise(irts{i,j}, -1);
        toc
    end
end

%{
for i = 1:4
    for j = 1:3
        irts{i,j}.tau = irts{i,j}.l2_err.^2;
    end
end
%}

% sample
z = rand(d, 1E4);
for i = 1:length(bases)
    for j = 1:3
        figure
        tic;[r,f] = eval_irt(irts{i,j}, z);toc
        tic;z0 = eval_rt(irts{i,j}, r);toc
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
end

