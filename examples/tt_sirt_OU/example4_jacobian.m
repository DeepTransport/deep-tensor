
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

opts{1} = FTTOption('tt_method', 'amen', ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 4);
opts{2} = FTTOption('tt_method', 'random', ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 4);


for i = 1:4
    for j = 1:2
        tic;
        irts{i,j} = TTSIRT(func, bases{i}, opts{j}, 'debug', deb, 'samples', sample_x);
        toc
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = rand(d, 1E2);
for i = 1:4
    for j = 1:2
        irts{i,j} = set_defensive(irts{i,j}, 5E-2);
        debug_jac(irts{i,j}, z, 1);
        debug_jac(irts{i,j}, z, -1);
    end
end

