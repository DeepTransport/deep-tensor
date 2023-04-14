
% setup the OU process
d = 4;
a = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = setup_ou_process(d, a);
func = @(x) eval_ou_process(data, x);

%%%%

debug_size = 1E4;
tmp = data.B\randn(d, debug_size);
tmp_f = feval(func, tmp);
deb = Debugger(tmp, tmp_f);
sample_x = data.B\randn(d, 1E3);

%%%%

% setup the reference polynomial

p = 50;
c1 = Chebyshev1st(p);
c2 = Chebyshev2nd(p);
le = Legendre(p);

bases{1} = ApproxBases(le, BoundedDomain([-4,4]), d);
bases{2} = ApproxBases(c1, LogarithmicMapping(4), d);

%
%{
opt = SparseOption();
opt.tol = 1e-2;
opt.init_total_degree = 3;
opt.max_dim_basis = 2e3;
opt.max_sample_size = 2e4;
opt.enrich_degree = 1;
opt.init_sample_size = 3;
opt.enrich_sample_size = 1;
opt.adaptation_rule = 'reducedmargin';


for i = 1:length(bases)
    tic;
    irts{i} = SparseSIRT(func, bases{i}, opt, 'debug', deb);
    toc
end
%}

% test m = 1, m = d-1, m = 8
% sample
z = rand(d, 1E4);
for i = 1:length(bases)
    m = 2;
    ind1 = 1:m;
    ind2 = (m+1):d;
    xc1 = deb.samples(ind1,:);
    xc2 = deb.samples(ind2,:);
    % conditional sample
    zc1 = z(ind1, :);
    zc2 = z(ind2, :);
    %
    if irts{i}.int_dir ~= 1, irts{i} = marginalise(irts{i}, 1); end
    tic;[rc2,f] = eval_cirt(irts{i}, xc1, zc2);toc
    %
    figure;
    subplot(2,2,1);plot(abs(exp(-func2(xc1, rc2, 1)) - exp(-f))); title('actual function value vs fft')
    subplot(2,2,2);plot(data.C(ind2, ind2) - cov(rc2')); title('actual covariance vs sample covariance')
    %
    disp(' ')
    disp(['approx eror: ' num2str(norm(exp(-func2(xc1, rc2, 1)) - exp(-f)))])
    disp(['cov eror: ' num2str(norm(data.C(ind2, ind2) - cov(rc2'))/norm(data.C))])
    disp(' ')
    
    if irts{i}.int_dir ~= -1, irts{i} = marginalise(irts{i}, -1); end
    tic;[rc1,f] = eval_cirt(irts{i}, xc2, zc1);toc
    %
    subplot(2,2,3);plot(abs(exp(-func2(rc1, xc2, -1)) - exp(-f))); title('actual function value vs fft')
    subplot(2,2,4);plot(data.C(ind1, ind1) - cov(rc1')); title('actual covariance vs sample covariance')
    disp(' ')
    disp(['approx eror: ' num2str(norm(exp(-func2(rc1, xc2, -1)) - exp(-f)))])
    disp(['cov eror: ' num2str(norm(data.C(ind1, ind1) - cov(rc1'))/norm(data.C))])
    disp(' ')
end

for i = 1:4
    irts{i} = set_defensive(irts{i}, 1E-2);
    m = 2;
    indx = 1:m;
    indy = (m+1):d;
    x = deb.samples(indx,1);
    my = data.C(indy,indx)*(data.C(indx,indx)\x);
    Cy = data.C(indy,indy) - data.C(indy,indx)*(data.C(indx,indx)\data.C(indx,indy));
    % conditional sample
    zy = z(indy, :);
    %
    if irts{i}.int_dir ~= 1, irts{i} = marginalise(irts{i}, 1); end
    tic;[y,f] = eval_cirt(irts{i}, x, zy);toc
    tic;
    zx = eval_rt (irts{i}, x);
    fx = eval_potential(irts{i}, x);
    [r,fxy] = eval_irt(irts{i}, [repmat(zx, 1, size(zy,2)); zy]);
    y2 = r(indy,:);
    f2 = fxy - fx;
    toc
    disp([norm(y2-y, 'fro'), norm(f2-f)])
    
    %
    figure;
    subplot(2,3,1);plot(abs(exp(-func2(repmat(x,1,size(zy,2)), y, 1)) - exp(-f))); title('actual function value vs fft')
    subplot(2,3,2);plot(Cy - cov(y')); title('actual covariance vs sample covariance')
    subplot(2,3,3);plot(my - mean(y,2)); title('actual mean vs sample mean')
    
    disp(' ')
    disp(['approx eror: ' num2str(norm(exp(-func2(repmat(x,1,size(zy,2)), y, 1)) - exp(-f)))])
    disp(['cov eror: ' num2str(norm(Cy - cov(y'))/norm(data.C))])
    disp(' ')
    
    m = 2;
    indy = 1:m;
    indx = (m+1):d;
    x = deb.samples(indx,1);
    my = data.C(indy,indx)*(data.C(indx,indx)\x);
    Cy = data.C(indy,indy) - data.C(indy,indx)*(data.C(indx,indx)\data.C(indx,indy));
    % conditional sample
    zy = z(indy, :);
    %
    if irts{i}.int_dir ~= -1, irts{i} = marginalise(irts{i}, -1); end
    tic;[y,f] = eval_cirt(irts{i}, x, zy);toc
    %
    subplot(2,3,4);plot(abs(exp(-func2(y, repmat(x,1,size(zy,2)), -1)) - exp(-f))); title('actual function value vs fft')
    subplot(2,3,5);plot(Cy - cov(y')); title('actual covariance vs sample covariance')
    subplot(2,3,6);plot(my - mean(y,2)); title('actual mean vs sample mean')
    
    disp(' ')
    disp(['approx eror: ' num2str(norm(exp(-func2(y, repmat(x,1,size(zy,2)), -1)) - exp(-f)))])
    disp(['cov eror: ' num2str(norm(Cy - cov(y'))/norm(data.C))])
    disp(' ')
end
