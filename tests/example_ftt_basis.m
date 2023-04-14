
% setup the reference
tol = 1E-10;
dom = BoundedDomain([-1,1]);
d = 10;

bases1 = ApproxBases(Legendre(20), dom, d);
bases2 = ApproxBases(Lagrange1(50), dom, d);
bases3 = ApproxBases(Lagrangep(5,4), dom, d);

%%%%

func1 = @(x) sqrt(1./sum(1E-5+x.^2,1));
func2 = @(x) [sqrt(1./sum(1E-5+x.^2,1)); sqrt(1./sum(1E-2+x.^2,1))];
func3 = @(x) [(1 + sum(x,1)).^(-d-1); exp( - sum(abs(x - 0.5), 1)); cos( sum(x.^2,1) )];

func = @(x) func2( (x-tol)*(1-tol)*2 -1 ); % rescale to the reference [-1, 1] domain

%%%%


debug_size = 1E4;
debug_x = rand(d, debug_size)*2-1;

deb = Debugger(debug_x);

opt1 = FTTOption('max_als', 6, 'als_tol', 1E-8, 'local_tol', 1E-20, 'kick_rank', 2, 'init_rank', 6, 'max_rank', 12);
opt2 = FTTOption('tt_method', 'random', 'max_als', 5, 'als_tol', 1E-8, 'local_tol', 1E-20, 'kick_rank', 2, 'init_rank', 6, 'max_rank', 12);

bases = {bases1, bases2, bases3};
for i = 1:3
    tic;
    ftt{i,1} = FTT(func, bases{i}, opt1, 'debug', deb);
    ftt{i,2} = round(ftt{i,1}, 1E-4);
    ftt{i,3} = FTT(func, ftt{i,1}, 'debug', deb);
    ftt{i,4} = FTT(func, ftt{i,2}, 'debug', deb);
    
    ftt{i,5} = FTT(func, bases{i}, opt2, 'debug', deb);
    ftt{i,6} = round(ftt{i,5}, 1E-4);
    ftt{i,7} = FTT(func, ftt{i,5}, 'debug', deb);
    ftt{i,8} = FTT(func, ftt{i,6}, 'debug', deb);
    toc
end

%{
if FTT.direction > 0
    plot_ftt_1d(FTT, func, d)
else
    plot_ftt_1d(FTT, func, 1)
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
        tic; approx{i,j} = eval_reference(ftt{i,j}, debug_x); toc
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

