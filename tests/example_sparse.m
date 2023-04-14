
% setup the reference
tol = 1E-10;
dom = BoundedDomain([tol,1]);
d = 3;

base = ApproxBases(Legendre(20), dom, d);

%%%%

func1 = @(x) sqrt(1./sum(1E-5+x.^2,1));
func2 = @(x) [sqrt(1./sum(1E-5+x.^2,1)); sqrt(1./sum(1E-2+x.^2,1))];
func3 = @(x) [(1 + sum(x,1)).^(-d-1); exp( - sum(abs(x - 0.5), 1)); cos( sum(x.^2,1) )];

func = func1;

%%%%

debug_size = 1E4;
debug_x = sample_measure(base, debug_size);
sample_size = 1E4;
sample_x = sample_measure(base, sample_size);


opt1 = SparseOption('method', 'least_square');
opt2 = SparseOption('method', 'least_square', 'adapt', true);
opt3 = SparseOption('method', 'magic_points');

spfun{1} = SparseFun(func, base, opt1, 'debug_x', debug_x);
spfun{2} = SparseFun(func, base, opt1, 'sample_x', sample_x, 'debug_x', debug_x);
spfun{3} = SparseFun(func, base, opt3, 'debug_x', debug_x);
spfun{4} = SparseFun(func, base, opt3, 'sample_x', sample_x, 'debug_x', debug_x);
%
spfun{5} = SparseFun(func, base, opt2, 'debug_x', debug_x);
spfun{6} = SparseFun(func, base, opt2, 'sample_x', sample_x, 'debug_x', debug_x);

tic; exact = func(debug_x); toc
for i = 1:6
    tic; app{i} = eval(spfun{i}, debug_x); toc
end

for i = 1:6
    figure
    plot(exact(:) - app{i}(:), '.')
end

x = debug_x(:,10);
for i = 1:6
    [gx,fx] = grad(spfun{i}, x);
    f = eval(spfun{i}, x);
    tol = 1E-5;
    fp = zeros(size(gx));
    fm = zeros(size(gx));
    for ii = 1:length(x)
        xp = x;
        xp(ii) = xp(ii)+tol;
        xm = x;
        xm(ii) = xm(ii)-tol;
        fp(ii) = eval(spfun{i}, xp);
        fm(ii) = eval(spfun{i}, xm);
    end
    disp(norm((fp-fm)/(2*tol) - gx))
end


