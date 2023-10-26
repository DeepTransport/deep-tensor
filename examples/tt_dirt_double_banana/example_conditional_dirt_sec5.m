close all

% define the joint density
sig = 0.3;
model = ConditionalBanana(sig);

dom  = BoundedDomain([-5,5]);
diag = GaussReference(0, 1, dom);
base = ApproxBases(Lagrangep(2,30), dom, 3);
% opt = FTToption('max_als', 2, 'als_tol', 1E-8, 'local_tol', 1E-5, 'kick_rank', 3, 'init_rank', 30, 'max_rank', 50);
% irt = DIRT(fun, 3, poly, refmap, opt, 'min_beta', 1E-4, 'ess_tol', 0.5, 'method', 'Aratio');

temp = Tempering1('min_beta', 1E-3, 'ess_tol', 0.5);

sirt_opt = TTOption('max_als', 2, 'als_tol', 1E-8, 'local_tol', 1E-5, 'kick_rank', 3, 'init_rank', 30, 'max_rank', 50);
dirt_opt = DIRTOption('method', 'Aratio');

if ~exist('irt')
    irt = TTDIRT(@(x)model.eval_potential_joint(x), base, temp, diag, sirt_opt, dirt_opt);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data = [-2, -1, 0, 1, 2, 3];

data = 0;

red    = '#D95319';
blue   = '#0072BD';
purple = '#7E2F8E';
grey = [0.7, 0.7, 0.7];

n  = 100;
xs = linspace(-4, 4, n);
ys = linspace(-4, 4, n);
[xx,yy] = meshgrid(xs, ys);
xts = [xx(:), yy(:)]';
%
rxs = linspace(domain_left(diag.domain), domain_right(diag.domain), n);
rys = linspace(domain_left(diag.domain), domain_right(diag.domain), n);
[xx,yy] = meshgrid(rxs, rys);
rts = [xx(:), yy(:)]';

for ii = 1:length(data)
    dat = data(ii); % data
    ry = eval_rt(irt, dat); % reference data sample
    
    % true conditional, unnormalised
    [mllkd,mlp] = model.eval_potential_joint([repmat(dat,1,size(xts,2));xts]);
    bf = exp(-mllkd-mlp);
    % conditional irt density in target space, unnormalised
    rf = eval_potential(irt, [repmat(dat,1,size(xts,2));xts]);
    
    figure('position', [100, 100, 1200, 600])
    subplot(2,4,1)
    contour(xs, ys, reshape(bf(:), n, n), 10, 'linewidth', 1, 'Color', blue)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('$\pi$', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
    
    subplot(2,4,5)
    contour(xs, ys, reshape(exp(-rf(:)), n, n), 10, 'linewidth', 1, 'Color', blue)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('$\hat\pi$', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
    
    % pullback density of the true conditinal
    mlf = pullback(irt, @(x)model.eval_potential_joint(x), [repmat(ry,1,size(rts,2));rts]);
    % the following does the same
    % mlf = model.pullback_potential_nuts(irt, ry, rts);
    %
    
    r = random(irt.ref, 2, 256); % sample the reference measure
    subplot(2,4,2)
    contour(rxs, rys, reshape(exp(-mlf(:)), n, n), 20, 'linewidth', 1, 'Color', blue)
    hold on
    plot(r(1,:), r(2,:), '.', 'MarkerSize', 10, 'Color', purple)
    xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('$T^\sharp \pi$, random', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-3, 3, -3, 3])
    
    xx = eval_irt(irt, [repmat(ry,1,size(r,2));r]);
    xx = xx(2:3,:);
    subplot(2,4,6)
    contour(xs, ys, reshape(bf(:), n, n), 10, 'linewidth', 1, 'Color', blue)
    hold on
    plot(xx(1,:), xx(2,:), '.', 'MarkerSize', 10, 'Color', purple)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('$\pi$', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
    
    
    r = sobol(irt.ref, 2, 256); % QMC sample the reference measure
    subplot(2,4,3)
    contour(rxs, rys, reshape(exp(-mlf(:)), n, n), 20, 'linewidth', 1, 'Color', blue)
    hold on
    plot(r(1,:), r(2,:), '.', 'MarkerSize', 10, 'Color', purple)
    xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('$T^\sharp \pi$, Sobol', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-3, 3, -3, 3])
    
    xx = eval_irt(irt, [repmat(ry,1,size(r,2));r]);
    xx = xx(2:3,:);
    subplot(2,4,7)
    contour(xs, ys, reshape(bf(:), n, n), 10, 'linewidth', 1, 'Color', blue)
    hold on
    plot(xx(1,:), xx(2,:), '.', 'MarkerSize', 10, 'Color', purple)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('$\pi$', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
    
    % regular grid in the reference measure
    [xx,yy] = meshgrid(linspace(-3,3,16), linspace(-3,3,16));
    r = [xx(:), yy(:)]';
    subplot(2,4,4)
    contour(rxs, rys, reshape(exp(-mlf(:)), n, n), 20, 'linewidth', 1, 'Color', blue)
    hold on
    plot(r(1,:), r(2,:), '.', 'MarkerSize', 10, 'Color', purple)
    xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('$T^\sharp \pi$, grid', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-3, 3, -3, 3])
    
    xx = eval_irt(irt, [repmat(ry,1,size(r,2));r]);
    xx = xx(2:3,:);
    subplot(2,4,8)
    contour(xs, ys, reshape(bf(:), n, n), 10, 'linewidth', 1, 'Color', blue)
    hold on
    plot(xx(1,:), xx(2,:), '.', 'MarkerSize', 10, 'Color', purple)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('$\pi$', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
    

    nsteps = 2^10;
    init = [0.9; 0];
    tic;
    out1 = NUTS(@(x)model.eval_potential_conditional(dat,x), init, nsteps);
    xx1 = out1.samples;
    toc
    %
    tic;
    out2 = NUTS(@(z)model.pullback_potential_nuts(irt,ry,z), init, nsteps);
    xs2 = eval_irt(irt, [repmat(ry,1,size(out2.samples,2));out2.samples]);
    xx2 = xs2(2:3,:);
    toc
    tic; % Metropolized independence sampler
    out3 = pCN(@(z)model.pullback_potential_pcn(irt,ry,z), init, nsteps, log(2));
    xs3 = eval_irt(irt, [repmat(ry,1,size(out3.samples,2));out3.samples]);
    xx3 = xs3(2:3,:);
    toc
    tic; % negatively correlated pCN
    out4 = pCN(@(z)model.pullback_potential_pcn(irt,ry,z), init, nsteps, log(10));
    xs4 = eval_irt(irt, [repmat(ry,1,size(out4.samples,2));out4.samples]);
    xx4 = xs4(2:3,:);
    toc   
    
    msteps = 1:5;
    figure('position', [100, 100, 1200, 600])
    subplot(4,4,[1,5])
    contour(xs, ys, reshape(bf(:), n, n), 10, 'linewidth', .5, 'Color', grey)
    hold on
    plot(xx1(1,:), xx1(2,:), '.', 'Color', blue)
    plot(xx1(1,msteps), xx1(2,msteps), '-o', 'linewidth', 1, 'Color', purple, 'MarkerFaceColor', purple)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('NUTS', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
        
    subplot(4,4,[2,6])
    contour(xs, ys, reshape(bf(:), n, n), 10, 'linewidth', .5, 'Color', grey)
    hold on
    plot(xx2(1,:), xx2(2,:), '.', 'Color', blue)
    plot(xx2(1,msteps), xx2(2,msteps), '-o', 'linewidth', 1, 'Color', purple, 'MarkerFaceColor', purple)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('preconditioned NUTS', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
    
    subplot(4,4,[3,7])
    contour(xs, ys, reshape(bf(:), n, n), 10, 'linewidth', .5, 'Color', grey)
    hold on
    plot(xx3(1,:), xx3(2,:), '.', 'Color', blue)
    plot(xx3(1,msteps), xx3(2,msteps), '-o', 'linewidth', 1, 'Color', purple, 'MarkerFaceColor', purple)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('Metropilis independence', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
    
    
    subplot(4,4,[4,8])
    contour(xs, ys, reshape(bf(:), n, n), 10, 'linewidth', .5, 'Color', grey)
    hold on
    plot(xx4(1,:), xx4(2,:), '.', 'Color', blue)
    plot(xx4(1,msteps), xx4(2,msteps), '-o', 'linewidth', 1, 'Color', purple, 'MarkerFaceColor', purple)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('preconditioned pCN', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
    
    subplot(4,4,9);  autocorr(xx1(1,:)); 
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    ylabel('$x_1$','interpreter', 'latex', 'fontsize', 16); 
    xlabel(''); title('')
    subplot(4,4,13); autocorr(xx2(2,:)); 
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    ylabel('$x_2$','interpreter', 'latex', 'fontsize', 16)
    xlabel('Lag','interpreter', 'latex', 'fontsize', 16)
    title('')
    %     
    subplot(4,4,10); autocorr(xx2(1,:));  
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    ylabel('$x_1$','interpreter', 'latex', 'fontsize', 16); 
    xlabel(''); title('')
    subplot(4,4,14); autocorr(xx2(2,:)); 
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    ylabel('$x_2$','interpreter', 'latex', 'fontsize', 16)
    xlabel('Lag','interpreter', 'latex', 'fontsize', 16)
    title('')
    %     
    subplot(4,4,11); autocorr(xx3(1,:)); 
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    ylabel('$x_1$','interpreter', 'latex', 'fontsize', 16); 
    xlabel(''); title('')
    subplot(4,4,15); autocorr(xx3(2,:)); 
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    ylabel('$x_2$','interpreter', 'latex', 'fontsize', 16)
    xlabel('Lag','interpreter', 'latex', 'fontsize', 16)
    title('')
    %
    subplot(4,4,12); autocorr(xx4(1,:)); 
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    ylabel('$x_1$','interpreter', 'latex', 'fontsize', 16); 
    xlabel('Lag','interpreter', 'latex', 'fontsize', 16)
    xlabel(''); title('')
    subplot(4,4,16); autocorr(xx4(2,:));
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    ylabel('$x_2$','interpreter', 'latex', 'fontsize', 16)
    xlabel('Lag','interpreter', 'latex', 'fontsize', 16)
    title('')
    
end
