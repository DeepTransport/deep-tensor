
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
dom = BoundedDomain([-5,5]);
bases{1} = ApproxBases(Legendre(40), dom, d);
bases{2} = ApproxBases(Fourier(20), dom, d);
bases{3} = ApproxBases(Lagrange1(40), dom, d);
bases{4} = ApproxBases(Lagrangep(5,8), dom, d);

opts{1} = FTTOption('tt_method', 'amen', 'max_als', 5, ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 5);
opts{2} = FTTOption('tt_method', 'random', 'max_als', 5, ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 5);


for i = 1:length(bases)
    for j = 1:2
        tic;
        irts{i,j} = TTSIRT(func, bases{i}, opts{j}, 'debug', deb, 'samples', sample_x);
        toc
    end
    tmp = round(irts{i,1}, 1E-2);
    irts{i,3} = TTSIRT(func, tmp, 'debug', deb, 'samples', sample_x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% should test ind = 1, ind = 1:(d-1) for > 0
% should test ind = d, ind = 2:d for < 0
% sample
z = rand(d, 1E4);
for i = 1:4
    for j = 1:3
        figure;
        % test 1
        ind  = 1:8;
        if irts{i,j}.int_dir ~= 1, irts{i,j} = marginalise(irts{i,j}, 1); end
        tic;[r,p] = eval_irt(irts{i,j}, z(ind,:));toc
        fx = eval_pdf(irts{i,j}, r);
        tic;z0 = eval_rt(irts{i,j}, r);toc
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
        if irts{i,j}.int_dir ~= -1, irts{i,j} = marginalise(irts{i,j}, -1); end
        tic;[r,p] = eval_irt(irts{i,j}, z(ind,:));toc
        fx = eval_pdf(irts{i,j}, r);
        tic;z0 = eval_rt(irts{i,j}, r);toc
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
end



