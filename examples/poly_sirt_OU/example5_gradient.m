

z = rand(d, 123);

for i = 1:length(bases)
    [r,~,g] = eval_irt(irts{i}, z);
    tol = 1e-8;
    fp = zeros(size(g));
    fm = zeros(size(g));
    for j = 1:d
        rp = r;
        rp(j,:) = rp(j,:) + tol;
        fp(j,:) = eval_potential(irts{i}, rp);
        
        rm = r;
        rm(j,:) = rm(j,:) - tol;
        fm(j,:) = eval_potential(irts{i}, rm);
    end
    gd = ( fp - fm ) / (tol*2);
    norm(g-gd,'fro')
end

for i = 1:length(bases)
    [r,f,g] = eval_irt_reference(irts{i}, z);
    tol = 1e-8;
    fp = zeros(size(g));
    fm = zeros(size(g));
    for j = 1:d
        rp = r;
        rp(j,:) = rp(j,:) + tol;
        fp(j,:) = eval_potential_reference(irts{i}, rp);
        
        rm = r;
        rm(j,:) = rm(j,:) - tol;
        fm(j,:) = eval_potential_reference(irts{i}, rm);
    end
    gd = ( fp - fm ) / (tol*2);
    norm(g-gd,'fro')
end

i = 2;
r = eval_irt_reference(irts{i}, z);
g = eval_measure_potential_reference_grad(irts{i}.approx.base, r);
tol = 1e-8;
fp = zeros(size(g));
fm = zeros(size(g));
for j = 1:d
    rp = r;
    rp(j,:) = rp(j,:) + tol;
    fp(j,:) = eval_measure_potential_reference(irts{i}.approx.base, rp);
    
    rm = r;
    rm(j,:) = rm(j,:) - tol;
    fm(j,:) = eval_measure_potential_reference(irts{i}.approx.base, rm);
end
gd = ( fp - fm ) / (tol*2);
norm(g-gd,'fro')


%{
for i = 1:length(bases)
    [r,f,g] = eval_irt_reference(irts{i}, z);
    J = eval_rt_jac_reference(irts{i}, r, z);
    for j = 1:size(z,2)
        ind = (j-1)*d + (1:d);
        % J(:,ind) % dzdx, dfdx/dzdx
        g(:,j) = J(:,ind)'\g(:,j);
    end
    
    tol = 1e-4;
    fp = zeros(size(g));
    fm = zeros(size(g));
    for j = 1:d
        zp = z;
        zp(j,:) = zp(j,:) + tol;
        [~,fp(j,:)] = eval_irt_reference(irts{i}, zp);
        
        zm = z;
        zm(j,:) = zm(j,:) - tol;
        [~,fm(j,:)] = eval_irt_reference(irts{i}, zm);
    end
    gd = ( fp - fm ) / (tol*2);
    norm(g-gd,'fro')
end
%}