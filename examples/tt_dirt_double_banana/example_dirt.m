
sig = 0.3;
%dat = [3.940392199905546; 4.403271869259551]; % 2 bananas
dat = [3; 5]; % 2 bananas
model = DoubleBanana(sig, dat);

fun = @(z) eval_potential_dirt(model, z);

% Gaussian
ref   = @(u) erfinv(u*2-1)*sqrt(2);

red    = '#D95319';
blue   = '#0072BD';
purple = '#7E2F8E';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% target

n  = 100;
xs = linspace(-4, 4, n);
ys = linspace(-4, 4, n);
[xx,yy] = meshgrid(xs, ys);
xts = [xx(:), yy(:)]';
const = 64/n^2;

[mllkd, mlp] = fun(xts);
bf = exp(-mllkd-mlp);
rf = exp(-0.5*sum(xts.^2, 1));

%{
h = figure(1);
set(h, 'position', [100, 100, 800, 400])
subplot(1,2,1)
contour(xs, ys, reshape(rf(:), n, n), 5, 'color', blue, 'linewidth', 2)
xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('$\mu$', 'interpreter', 'latex', 'fontsize', 20)
subplot(1,2,2)
contour(xs, ys, reshape(bf(:), n, n), 5, 'linewidth', 2)
xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('$\pi$', 'interpreter', 'latex', 'fontsize', 20)
colormap default
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% use more than one ALS with AMEN

% with uniform reference, the boundary layer may cause trouble
%diag = UniformReference(); 

% with Gaussian reference, the Legendre with single ALS iteration may not
% work, overall the Legendre basis does not work nicely with Gaussian ref
% 
d = 2;
dom = BoundedDomain([-4,4]);
diag = GaussReference(0, 1, dom);
%diag = UniformReference(BoundedDomain([-1,1]));
bases{1} = {ApproxBases(Legendre(50), dom, d), ApproxBases(Legendre(50), dom, d)};
bases{2} = ApproxBases(Fourier(30), dom, d);
bases{3} = {ApproxBases(Lagrange1(50), dom, d), ApproxBases(Lagrange1(100), dom, d)};
bases{4} = ApproxBases(Lagrangep(5,10), dom, d);

temp = Tempering1('min_beta', 1E-3, 'ess_tol', 0.3);

sirt_opt = TTOption('tt_method', 'random', 'max_als', 1, 'als_tol', 1E-8, 'local_tol', 1E-2, 'kick_rank', 2, 'init_rank', 40, 'max_rank', 50);

dirt_opt1 = DIRTOption('method', 'Aratio'); %, 'dhell_tol', 1E-3);
dirt_opt2 = DIRTOption('method', 'Eratio'); %, 'dhell_tol', 1E-3);

airt = TTDIRT(fun, bases{2}, temp, diag, sirt_opt, dirt_opt1);
eirt = TTDIRT(fun, bases{4}, temp, diag, sirt_opt, dirt_opt2);

n  = 100;
rxs = linspace(domain_left(diag.domain), domain_right(diag.domain), n);
rys = linspace(domain_left(diag.domain), domain_right(diag.domain), n);
[xx,yy] = meshgrid(rxs, rys);
rts = [xx(:), yy(:)]';

for k = 1:num_layers(airt)
    figure('position', [100, 100, 1200, 800])
    
    rf = eval_potential(airt, xts, k);
    rf = exp(-rf);
    [mllkd, mlp] = fun(xts);
    bf = exp(-mllkd*airt.bridge.betas(k)-mlp);
    bf = bf/(sum(bf(:))*const);
    
    subplot(2,2,1)
    contour(xs, ys, reshape(rf(:), n, n), 8, 'linewidth', 1)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$\hat\pi$', 'interpreter', 'latex', 'fontsize', 20)
    
    subplot(2,2,2)
    contour(xs, ys, reshape(bf(:), n, n), 8, 'linewidth', 1)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$\pi$', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    
    fk = eval_pdf(airt.irts{k}, rts);
    subplot(2,2,3)
    contour(rxs, rys, reshape(fk(:), n, n), 8, 'linewidth', 1)
    xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('tt', 'interpreter', 'latex', 'fontsize', 20)
    
    if k > 1
        [x,logf] = eval_irt(airt, rts, k-1);
        [mllkd, mlp] = fun(x);
        logfz = log_joint_pdf(diag, rts);
        
        bf = exp(-mllkd*(airt.bridge.betas(k)-airt.bridge.betas(k-1))+logfz);
        
        subplot(2,2,4)
        contour(rxs, rys, reshape(bf(:), n, n), 8, 'linewidth', 1)
        xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
        ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
        set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
        title('a. ratio fun', 'interpreter', 'latex', 'fontsize', 20)
    end
end


for k = 1:num_layers(eirt)
    figure('position', [100, 100, 1200, 800])
    
    rf = eval_potential(eirt, xts, k);
    rf = exp(-rf);
    [mllkd, mlp] = fun(xts);
    bf = exp(-mllkd*eirt.bridge.betas(k)-mlp);
    bf = bf/(sum(bf(:))*const);
    
    subplot(2,2,1)
    contour(xs, ys, reshape(rf(:), n, n), 8, 'linewidth', 1)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$\hat\pi$', 'interpreter', 'latex', 'fontsize', 20)
    
    subplot(2,2,2)
    contour(xs, ys, reshape(bf(:), n, n), 8, 'linewidth', 1)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$\pi$', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    
    fk = eval_pdf(eirt.irts{k}, rts);
    subplot(2,2,3)
    contour(rxs, rys, reshape(fk(:), n, n), 8, 'linewidth', 1)
    xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('tt', 'interpreter', 'latex', 'fontsize', 20)
    
    if k > 1
        [x,logf] = eval_irt(eirt, rts, k-1);
        [mllkd, mlp] = fun(x);
        logfz = log_joint_pdf(diag, rts);
        
        bf = exp(-mllkd*eirt.bridge.betas(k)-mlp+logf+logfz);
        
        subplot(2,2,4)
        contour(rxs, rys, reshape(bf(:), n, n), 8, 'linewidth', 1)
        xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
        ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
        set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
        title('e. ratio fun', 'interpreter', 'latex', 'fontsize', 20)
    end
end
