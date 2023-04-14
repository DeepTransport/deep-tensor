
sig = 0.3;
%dat = [3.940392199905546; 4.403271869259551]; % 2 bananas
dat = [3; 5]; % 2 bananas
fun = @(z,arg) fun_banana(z, dat, sig, 1);

% Gaussian
ref   = @(u) erfinv(u*2-1)*sqrt(2);

red    = '#D95319';
blue   = '#0072BD';
purple = '#7E2F8E';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% target

n  = 100;
xs = linspace(-2, 2, n);
ys = linspace(-1.5, 2.5, n);
[xx,yy] = meshgrid(xs, ys);
xts = [xx(:), yy(:)]';
const = 64/n^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = 2;
dom = BoundedDomain([-4,4]);
diag = GaussReference(0, 1, dom);
%base = ApproxBases(Lagrangep(5,20), dom, d);
base = ApproxBases(Fourier(30), dom, d);
temp = Tempering1('betas', [0.1, 1]);

opt1 = FTTOption('tt_method', 'amen', ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'init_rank', 12, 'max_rank', 12, 'max_als', 4);

dirt = TTDIRT(fun,base,temp,diag,opt1);


opt2 = FTTOption('tt_method', 'amen', ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'init_rank', 24, 'max_rank', 26, 'max_als', 4);
sirt1 = TTSIRT(fun,base,opt2);

%{
base2 = ApproxBases(Lagrangep(5,20), dom, d);
sirt2 = TTSIRT(fun,base2,opt);
%}

n  = 100;
rxs = linspace(-2.5, 2.5, n);
rys = linspace(-2, 2, n);
[xx,yy] = meshgrid(rxs, rys);
rts = [xx(:), yy(:)]';


figure
[mllkd, mlp] = fun(xts);
bf = exp(-mllkd-mlp);
bf = bf/(sum(bf(:))*const);

C = contour(xs, ys, reshape(bf(:), n, n), 8, 'linewidth', 1);
C = C';
save(['contour_double_banana.dat'], 'C', '-ascii', '-tabs' );


figure
fk = eval_pdf(sirt1, xts);
C = contour(xs, ys, reshape(fk(:), n, n), 8, 'linewidth', 1);
C = C';
save(['contour_double_banana_sirt.dat'], 'C', '-ascii', '-tabs' );
    
%{
figure
fk = eval_pdf(sirt2, rts);
contour(rxs, rys, reshape(fk(:), n, n), 8, 'linewidth', 1)
xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('tt', 'interpreter', 'latex', 'fontsize', 20)
%}  

figure
[x,logf] = eval_irt(dirt, rts, 1);
[mllkd, mlp] = fun(x);
logfz = log_joint_pdf(dirt.ref, rts);

bf = exp(-mllkd*(dirt.bridge.betas(2)-dirt.bridge.betas(1))+logfz);

C = contour(rxs, rys, reshape(bf(:), n, n), 8, 'linewidth', 1);
C = C';
save(['contour_double_banana_pull.dat'], 'C', '-ascii', '-tabs' );


figure
rf = eval_potential(dirt, xts);
rf = exp(-rf);

C = contour(xs, ys, reshape(rf(:), n, n), 8, 'linewidth', 1);
C = C';
save(['contour_double_banana_dirt.dat'], 'C', '-ascii', '-tabs' );


