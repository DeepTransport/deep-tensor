
% setup the reference
tol = 1E-10;
dom = BoundedDomain([tol,1]);
d = 10;

bases1 = ApproxBases(Legendre(20), dom, d);
bases2 = ApproxBases(Lagrange1(50), dom, d);
bases3 = ApproxBases(Lagrangep(5,4), dom, d);
%%%%


funcin = @(x) sqrt(1./sum(1E-5+x.^2,1));


func = @(x) funcin( (x-tol)*(1-tol)*2 - 1 ); % rescale to the reference [-1, 1] domain

%%%%


debug_size = 1E4;
debug_x = rand(d, debug_size)*2-1;
deb = Debugger(debug_x);

opt1 = FTTOption('max_als', 4, 'als_tol', 1E-8, 'local_tol', 1E-10, 'kick_rank', 2, 'init_rank', 6, 'max_rank', 12);
opt2 = FTTOption('max_als', 5, 'als_tol', 1E-8, 'local_tol', 1E-10, 'kick_rank', 2, 'init_rank', 6, 'max_rank', 12);

bases = {bases1, bases2, bases3};
for i = 1:3
    tic;
    ftt{i,1} = FTT(func, bases{i}, opt1, 'debug', deb);
    ftt{i,2} = round(ftt{i,1}, 1E-4);
    
    ftt{i,3} = FTT(func, bases{i}, opt2, 'debug', deb);
    ftt{i,4} = round(ftt{i,3}, 1E-4);
    toc
end

x = debug_x(:,1);
for i = 1:3
    for j= 1:4
        gx = grad_reference(ftt{i,j}, x);
        tol = 1E-4;
        fp = zeros(size(gx));
        fm = zeros(size(gx));
        for ii = 1:length(x)
            xp = x;
            xp(ii) = xp(ii)+tol;
            xm = x;
            xm(ii) = xm(ii)-tol;
            fp(ii) = eval_reference(ftt{i,j}, xp);
            fm(ii) = eval_reference(ftt{i,j}, xm);
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