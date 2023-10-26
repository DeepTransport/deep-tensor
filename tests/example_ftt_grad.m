
% setup the reference
tol = 1E-10;
dom = BoundedDomain([tol,1]);
d = 10;

bases1 = ApproxBases(Legendre(20), dom, d);
bases2 = ApproxBases(Lagrange1(50), dom, d);
bases3 = ApproxBases(Lagrangep(5,4), dom, d);

debug_size = 1E4;
debug_x = sample_measure(bases1, debug_size);
deb = InputData([],debug_x);

%%%%

func = @(x) sqrt(1./sum(1E-5+x.^2,1));

%%%%

opt1 = TTOption('max_als', 4, 'als_tol', 1E-8, 'local_tol', 1E-10, 'kick_rank', 2, 'init_rank', 6, 'max_rank', 12);
opt2 = TTOption('max_als', 5, 'als_tol', 1E-8, 'local_tol', 1E-10, 'kick_rank', 2, 'init_rank', 6, 'max_rank', 12);

bases = {bases1, bases2, bases3};
for i = 1:3
    f = @(z) ApproxFun.feval_reference(bases{i}, func, z);
    tic;
    ftt{i,1} = TTFun(f, bases{i}, opt1, 'var', deb);
    ftt{i,2} = round(ftt{i,1}, 1E-4);
    
    ftt{i,3} = TTFun(f, bases{i}, opt2, 'var', deb);
    ftt{i,4} = round(ftt{i,3}, 1E-4);
    toc
end

x = debug_x(:,1);
for i = 1:3
    for j= 1:4
        gx = grad(ftt{i,j}, x);
        tol = 1E-4;
        fp = zeros(size(gx));
        fm = zeros(size(gx));
        for ii = 1:length(x)
            xp = x;
            xp(ii) = xp(ii)+tol;
            xm = x;
            xm(ii) = xm(ii)-tol;
            fp(ii) = eval(ftt{i,j}, xp);
            fm(ii) = eval(ftt{i,j}, xm);
        end
        disp(norm((fp-fm)/(2*tol) - gx))
    end
end


tol = 1E-4;
fpe = zeros(size(gx));
fme = zeros(size(gx));
for ii = 1:length(x)
    xp = x;
    xp(ii) = xp(ii)+tol;
    xm = x;
    xm(ii) = xm(ii)-tol;
    fpe(ii) = func(xp);
    fme(ii) = func(xm);
end
gxe = (fp-fm)/(2*tol);