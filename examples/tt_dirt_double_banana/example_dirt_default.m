
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = 2;
dom = BoundedDomain([-4,4]);
base = ApproxBases(Lagrange1(50), dom, d);
temp = Tempering1();


airt = TTDIRT(fun,base,temp);

n  = 100;
rxs = linspace(domain_left(dom), domain_right(dom), n);
rys = linspace(domain_left(dom), domain_right(dom), n);
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
        logfz = log_joint_pdf(airt.ref, rts);
        
        bf = exp(-mllkd*(airt.bridge.betas(k)-airt.bridge.betas(k-1))+logfz);
        
        subplot(2,2,4)
        contour(rxs, rys, reshape(bf(:), n, n), 8, 'linewidth', 1)
        xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
        ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
        set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
        title('a. ratio fun', 'interpreter', 'latex', 'fontsize', 20)
    end
end

