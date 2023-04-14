% define the joint density
sig = 0.3;
fun = @(z) joint_banana(z,sig);

dom = BoundedDomain([-5,5]);
diag = GaussReference(0, 1, dom);
base = ApproxBases(Lagrangep(2,30), dom, 3);
% opt = FTToption('max_als', 2, 'als_tol', 1E-8, 'local_tol', 1E-5, 'kick_rank', 3, 'init_rank', 30, 'max_rank', 50);
% irt = DIRT(fun, 3, poly, refmap, opt, 'min_beta', 1E-4, 'ess_tol', 0.5, 'method', 'Aratio');

temp = Tempering1('min_beta', 1E-4, 'ess_tol', 0.5);

sirt_opt = FTTOption('max_als', 2, 'als_tol', 1E-8, 'local_tol', 1E-5, 'kick_rank', 3, 'init_rank', 30, 'max_rank', 50);
dirt_opt = DIRTOption('method', 'Aratio');

irt = TTDIRT(fun, base, temp, diag, sirt_opt, dirt_opt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data = [-2, -1, 0, 1, 2, 3];

data = 0;

red    = '#D95319';
blue   = '#0072BD';
purple = '#7E2F8E';

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
    ry = eval_rt(irt, dat); % reference sample
    
    % true conditional, unnormalised
    [mllkd,mlp] = fun([repmat(dat,1,size(xts,2));xts]);
    bf = exp(-mllkd-mlp);
    % conditional irt density in target space, unnormalised
    rf = eval_potential(irt, [repmat(dat,1,size(xts,2));xts]);
    
    figure('position', [100, 100, 800, 800])
    subplot(2,2,1)
    contour(xs, ys, reshape(bf(:), n, n), 5, 'linewidth', 2)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$\pi$', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    
    subplot(2,2,2)
    contour(xs, ys, reshape(exp(-rf(:)), n, n), 5, 'linewidth', 2)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$\hat\pi$', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    
    % conditional irt density in reference space, unnormalised
    [x,f] = eval_irt(irt, [repmat(ry,1,size(rts,2));rts]);
    % pullback density of the true conditinal
    mlf = pullback(irt, fun, [repmat(ry,1,size(rts,2));rts]);
    %
    subplot(2,2,3)
    contour(rxs, rys, reshape(exp(-f(:)), n, n), 5, 'linewidth', 2)
    xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$\hat\pi(u)$', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    
    subplot(2,2,4)
    contour(rxs, rys, reshape(exp(-mlf(:)), n, n), 5, 'linewidth', 2)
    xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$T^\sharp \hat\pi$', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    
    init = randn(2,1);
    tic;
    for irun=1:16
        out1 = NUTS(@(x) log_target(fun,dat,x), init, 2^12);
        xx1(:,:,irun) = out1.samples;
    end
    toc
    %
    tic;
    for irun=1:16
        out2 = NUTS(@(z) log_target_pullback_nuts(irt,fun,ry,z), init, 2^12);
        xs2 = eval_irt(irt, [repmat(ry,1,size(out2.samples,2));out2.samples]);
        xx2(:,:,irun) = xs2(2:3,:);
    end
    toc
    tic;
    for irun=1:16
        out3 = pCN(@(z) log_target_pullback_pcn(irt,fun,ry,z), init, 2^12, log(10));
        xs3 = eval_irt(irt, [repmat(ry,1,size(out3.samples,2));out3.samples]);
        xx3(:,:,irun) = xs3(2:3,:);
    end
    toc
    tic;
    for irun=1:16
        out4 = pCN(@(z) log_target_pullback_pcn(irt,fun,ry,z), init, 2^12, log(2));
        xs4 = eval_irt(irt, [repmat(ry,1,size(out4.samples,2));out4.samples]);
        xx4(:,:,irun) = xs4(2:3,:);
    end
    toc    
    
    % xx = reshape(xx, 2, [], 2^4); 
    % std([mean(xx,2); std(xx,0,2)],0,3)       % N = 2^12 x 2^4
    % xx3 (pCN with sigma=log(10))      xx4 (pCN with sigma=log(2))     xx1 (original NUTS)     xx2 (pNUTS)
    % 6.665e-03                         1.516e-02                       1.989e-02               1.307e-02 
    % 1.130e-02                         2.188e-02                       3.266e-01               2.650e-02
    % 7.976e-03                         7.248e-03                       4.358e-02               1.368e-02
    % 1.104e-02                         2.701e-02                       4.762e-02               1.371e-02

    
    % 7.803e-03                         8.179e-03                       2.215e-02
    % 1.295e-02                         1.198e-02                       2.446e-02
    % 7.378e-03                         6.455e-03                       2.011e-02
    % 1.048e-02                         8.109e-03                       3.255e-02    
    % 
    % std([mean(xx.*double(xx(2,:,:)>xx(1,:,:).^2), 2); std(xx.*double(xx(2,:,:)>xx(1,:,:).^2),0, 2)], 0,3)
    % (mean and std over the upper banana only)
    % xx3 (pCN with sigma=log(10))      xx4 (pCN with sigma=log(2))     xx1             xx2
    % 6.565e-03                         1.244e-02                       8.910e-03       9.190e-03
    % 1.007e-02                         2.069e-02                       2.331e-01       1.921e-02
    % 9.937e-03                         2.120e-02                       1.103e-01       1.354e-02
    % 1.199e-02                         3.476e-02                       6.272e-02       1.614e-02
    
    % 7.457e-03                         6.047e-03                       0
    % 8.309e-03                         7.818e-03                       0
    % 9.124e-03                         6.244e-03                       0
    % 1.308e-02                         7.932e-03                       0
    
    figure('position', [100, 100, 1200, 800])
    subplot(4,3,[1,4])
    plot(out1.samples(1,:), out1.samples(2,:), '.')
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('HMC', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    axis([-4, 4, -4, 4])
        
    subplot(4,3,[2,5])
    plot(xx2(1,:), xx2(2,:), '.')
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('HMC, preconditioned', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    axis([-4, 4, -4, 4])
    
    subplot(4,3,[3,6])
    plot(xx3(1,:), xx3(2,:), '.')
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('pCN, preconditioned', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    axis([-4, 4, -4, 4])
    
%     subplot(4,3,7);  autocorr(out1.samples(1,:)); title('$x_1$','interpreter', 'latex'); 
%     subplot(4,3,10); autocorr(out1.samples(2,:)); title('$x_2$','interpreter', 'latex')
%     
%     subplot(4,3,8);  autocorr(xx2(1,:)); title('$x_1$','interpreter', 'latex'); 
%     subplot(4,3,11); autocorr(xx2(2,:)); title('$x_2$','interpreter', 'latex')
%     
%     subplot(4,3,9);  autocorr(xx3(1,:)); title('$x_1$','interpreter', 'latex'); 
%     subplot(4,3,12); autocorr(xx3(2,:)); title('$x_2$','interpreter', 'latex')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%debug the gradient
x = [0,0,1]';
[~,~,g1,g2] = fun(x);
g = g1+g2;
tol = 1E-5;
fp = zeros(size(g));
fm = zeros(size(g));
for i = 1:3
    xp = x;
    xp(i) = xp(i)+tol;
    xm = x;
    xm(i) = xm(i)-tol;
    [f1,f2] = fun(xp);
    fp(i) = f1+f2;
    [f1,f2] = fun(xm);
    fm(i) = f1+f2;
end
norm((fp-fm)/(2*tol) - g)

x = [dat,0,1]';
z = eval_rt(irt, x);
%
[~,~,g] = eval_irt(irt,z);
tol = 1E-5;
fp = zeros(size(g));
fm = zeros(size(g));
for i = 1:3
    zp = z;
    zp(i) = zp(i)+tol;
    zm = z;
    zm(i) = zm(i)-tol;
    [~,fp(i)] = eval_irt(irt,zp);
    [~,fm(i)] = eval_irt(irt,zm);
end
norm((fp-fm)/(2*tol) - g)
%
ry = eval_rt(irt, dat);
z = z(2:3);
[f,g] = log_target_pullback(irt,fun,ry,z);
tol = 1E-8;
fp = zeros(size(g));
fm = zeros(size(g));
for i = 1:2
    zp = z;
    zp(i) = zp(i)+tol;
    zm = z;
    zm(i) = zm(i)-tol;
    fp(i) = log_target_pullback(irt,fun,ry,zp);
    fm(i) = log_target_pullback(irt,fun,ry,zm);
end
norm((fp-fm)/(2*tol) - g)
%}
