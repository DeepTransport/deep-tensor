clear all
%addpath(genpath('.'));

%% Generate the target = banana
d = 2;
xmin = -4;
xmax = 4;
nx = 100;

gridx = linspace(xmin,xmax,nx);
[gridxx,gridyy] = meshgrid(gridx,gridx);
X = [gridxx(:),gridyy(:)]';

sigma = .5;
a = 1;
func = @(x,mask) banana_dirt(x, sigma, a);
% plot
clf
fplot = reshape(func(X),nx,nx);
contourf(gridx,gridx,exp(-fplot))

%% L-layer DIRT
nlayer = 5;
betas = 10.^(-2:2/(nlayer-1):0)
temp1 = Tempering1(betas);
temp2 = Tempering1('min_beta', 1E-3, 'ess_tol', 0.6);

p = 30;
dom = BoundedDomain([xmin,xmax]);
base = ApproxBases(Legendre(p), dom, d);
diag = UniformReference(dom);
%diag = GaussReference(0, 1, dom);

opt = SparseOption();
opt.tol = 1e-2;
opt.init_total_degree = 3;
opt.max_dim_basis = 200;
opt.max_sample_size = 2e3;
opt.enrich_degree = 2;
opt.init_sample_size = 2;
opt.enrich_sample_size = 1; 
opt.display_iterations = false;
opt.adaptation_rule = 'margin'; %'reducedmargin';
opt.fast = true;

dirt_opt = DIRTOption('method', 'Eratio', 'defensive', 1E-4);

%eirt = SparseDIRT(func, base, diag, s, dirt_opt);%, 'betas', betas);
poly_dirt = SparseDIRT(func, base, temp1, diag, opt, dirt_opt); %, 'betas', betas);


tt_opt = FTTOption('tt_method', 'random', 'max_als', 1, 'als_tol', 1E-8, 'local_tol', 1E-2, 'kick_rank', 2, 'init_rank', 5, 'max_rank', 10);
tt_dirt = TTDIRT(func, base, temp1, diag, tt_opt, dirt_opt);

%% plot
n  = 100;
rxs = linspace(domain_left(diag.domain), domain_right(diag.domain), n);
rys = linspace(domain_left(diag.domain), domain_right(diag.domain), n);
[xx,yy] = meshgrid(rxs, rys);
rts = [xx(:), yy(:)]';
const = (domain_right(diag.domain)-domain_left(diag.domain)^2)/n^2;

eirt = poly_dirt;

for k = 1:num_layers(eirt)
    figure('position', [100, 100, 1200, 800])
    
    rf = eval_potential(eirt, rts, k);
    rf = exp(-rf);
    [mllkd, mlp] = func(rts);
    bf = exp(-mllkd*eirt.bridge.betas(k)-mlp);
    bf = bf/(sum(bf(:))*const);
    
    subplot(2,2,1)
    contour(rxs, rys, reshape(rf(:), n, n), 8, 'linewidth', 1)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$\hat\pi$', 'interpreter', 'latex', 'fontsize', 20)
    
    subplot(2,2,2)
    contour(rxs, rys, reshape(bf(:), n, n), 8, 'linewidth', 1)
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
        [x,mlogf] = eval_irt(eirt, rts, k-1);
        [mllkd, mlp] = func(x);
        logfz = log_joint_pdf(diag, rts);
        
        bf = exp(-mllkd*eirt.bridge.betas(k)-mlp+mlogf+logfz);
        
        subplot(2,2,4)
        contour(rxs, rys, reshape(bf(:), n, n), 8, 'linewidth', 1)
        xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
        ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
        set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
        title('e. ratio fun', 'interpreter', 'latex', 'fontsize', 20)
    end
end