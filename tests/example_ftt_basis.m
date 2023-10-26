
% setup the reference
tol = 1E-10;
dom = BoundedDomain([0,1]);
d = 10;

bases1 = ApproxBases(Legendre(20), dom, d);
bases2 = ApproxBases(Lagrange1(50), dom, d);
bases3 = ApproxBases(Lagrangep(5,4), dom, d);

debug_size = 1E4;
debug_x = sample_measure(bases1, debug_size);
deb = InputData([],debug_x);

%%%%

func1 = @(x) sqrt(1./sum(1E-5+x.^2,1));
func2 = @(x) [sqrt(1./sum(1E-5+x.^2,1)); sqrt(1./sum(1E-2+x.^2,1))];
func3 = @(x) [(1 + sum(x,1)).^(-d-1); exp( - sum(abs(x - 0.5), 1)); cos( sum(x.^2,1) )];

func = func3;

%%%%

opt1 = TTOption('max_als', 5, 'als_tol', 1E-8, 'local_tol', 1E-20, 'kick_rank', 2, 'init_rank', 6, 'max_rank', 12);
opt2 = TTOption('tt_method', 'random', 'max_als', 5, 'als_tol', 1E-8, 'local_tol', 1E-20, 'kick_rank', 2, 'init_rank', 6, 'max_rank', 12);

bases = {bases1, bases2, bases3};

for i = 1:3
    f = @(z) ApproxFun.feval_reference(bases{i}, func, z);
    tic;
    ftt{i,1} = TTFun(f, bases{i}, opt1, 'var', deb);
    ftt{i,2} = round(ftt{i,1}, 1E-4);
    ftt{i,3} = TTFun(f, ftt{i,1}, 'var', deb);
    ftt{i,4} = TTFun(f, ftt{i,2}, 'var', deb);
    
    ftt{i,5} = TTFun(f, bases{i}, opt2, 'var', deb);
    ftt{i,6} = round(ftt{i,5}, 1E-4);
    ftt{i,7} = TTFun(f, ftt{i,5}, 'var', deb);
    ftt{i,8} = TTFun(f, ftt{i,6}, 'var', deb);
    toc
end

%{
if TTFun.direction > 0
    plot_ftt_1d(TTFun, func, d)
else
    plot_ftt_1d(TTFun, func, 1)
end
%}

%{
figure
[xx, yy] = meshgrid(linspace(0, 1, 100), linspace(0, 1, 100));
f = reshape(sqrt(1./(xx(:).^2 + yy(:).^2 + 0.001^2)), 100, 100);
surf(f)
%}
%
tic; exact = func(debug_x); toc
for i = 1:3
    for j = 1:8
        tic; approx{i,j} = eval(ftt{i,j}, debug_x); toc
    end
end

for i = 1:3
    figure
    for j = 1:8
        plot(exact(:) - approx{i,j}(:), '.')
        hold on
    end
    legend('amen', 'amen rounded', 'amen rebuild', 'amen rounded rebuild', 'rand', 'rand rounded' , 'random rebuild', 'random rounded rebuild')
end

