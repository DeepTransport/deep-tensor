
% setup the OU process
d = 20;
a = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = OU(d, a);
func = @(x) eval_potential(data, x);

%%%%

debug_size = 1E4;
debug_x = data.B\randn(d, debug_size);
sample_x = data.B\randn(d, 1E3);
deb = InputData(sample_x, debug_x);

%%%%

% setup the reference polynomial

dom = BoundedDomain([-5,5]);
bases{1} = ApproxBases(Legendre(40), dom, d);
bases{2} = ApproxBases(Fourier(20), dom, d);
bases{3} = ApproxBases(Lagrange1(40), dom, d);
bases{4} = ApproxBases(Lagrangep(5,8), dom, d);
%bases{5} = ApproxBases(Hermite(10), UnboundedDomain(), d);

opts{1} = TTOption('tt_method', 'amen', ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 1);
opts{2} = TTOption('tt_method', 'random', ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 1);

irts = cell(length(bases),3);
for i = 1:length(bases)
    for j = 1:2
        tic;
        deb.count.i = 0;
        irts{i,j} = TTSIRT(func, bases{i}, opts{j}, 'var', deb);
        toc
    end
    tmp = round(irts{i,1}.approx, 1E-2);
    deb.count.i = 0;
    irts{i,3} = TTSIRT(func, tmp, 'var', deb);
end

for i = 1:length(bases)
    for j = 1:3
        tic;
        irts{i,j} = marginalise(irts{i,j}, -1);
        toc
    end
end