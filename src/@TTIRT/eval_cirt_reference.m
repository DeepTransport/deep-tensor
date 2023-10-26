function [r,f] = eval_cirt_reference(obj, x, z)
%
%NEED TO BE UPDATED
%
% Evaluate the inverse of the conditional squared Rosenblatt transport 
% Y | X = R^{-1}(Z, X), where X is given, (X, Y) jointly follow the target 
% distribution represented by SIRT, Z is uniform. 
%   [Y,f] = EVAL_CIRT_REFERENCE(irt, X, Z)
%
%   Z - uniform random variables, m x n
%   X - given variables following the marginal of the target distribution
%   Y - random variable drawn from Y | X
%   f - pdf of Y | X
%
% The conditioning depends on irt.int_dir.
%   * With >0 (marginalised from right to left), X = (X_1, ..., X_k) is the
%     the first k coordinates in the joint and Y = (X_{k+1}, ... X_d)
%   * With <0 (marginalised from left to right), X = (X_m, ..., X_d) is the
%     last (d-m+1) coordinates in the joint and Y = (X_1, ... X_{m-1})

d = ndims(obj.approx);
nr = size(z,2);
dr = size(z,1);
nx = size(x,2);
dx = size(x,1);

if dx == 0 || dr == 0
    error('dimension of x or the dimension of z should be nonzero')
end

if dr + dx ~= d
    error('dimension of x and the dimension of z mismatch the dimension of the joint')
end

if nx ~= nr
    if nx == 1
        x = repmat(x, 1, nr);
        nx = nr;
    else
        error('number of x and the number of z mismatch')
    end
end

r = zeros(size(z));

% first compute the marginal
if obj.int_dir > 0
    %from the last dimension, marginalise to the first
    %order of the samples (x, r)
    %first evaluate the marginal for x
    frl = ones(nr,1);
    for k = 1:(dx-1)
        rkm = size(obj.approx.data.cores{k}, 1);
        % evaluate the updated basis function
        T2  = TTSIRT.eval_oned_core_213(obj.approx.base.oneds{k}, obj.approx.data.cores{k}, x(k,:));
        jj  = reshape(reshape(1:rkm*nx, rkm, nx)', [], 1);
        ii  = repmat((1:nx)', 1, rkm);
        B   = sparse(ii(:), jj(:), frl(:), nx, rkm*nx);
        frl = B*T2;
        %frl(isnan(frl)) = 0;
        %
    end
    k   = dx;
    rkm = size(obj.approx.data.cores{k}, 1);
    % evaluate the updated basis function
    % for the marginal
    T2  = TTSIRT.eval_oned_core_213(obj.approx.base.oneds{k}, obj.ys{k}, x(k,:));
    jj  = reshape(reshape(1:rkm*nx, rkm, nx)', [], 1);
    ii  = repmat((1:nx)', 1, rkm);
    B   = sparse(ii(:), jj(:), frl(:), nx, rkm*nx);
    frl_m   = B*T2;
    fm  = sum(frl_m.^2, 2)';
    %
    %for the frl, evaluate the updated basis function
    T2  = TTSIRT.eval_oned_core_213(obj.approx.base.oneds{k}, obj.approx.data.cores{k}, x(k,:));
    jj  = reshape(reshape(1:rkm*nx, rkm, nx)', [], 1);
    ii  = repmat((1:nx)', 1, rkm);
    B   = sparse(ii(:), jj(:), frl(:), nx, rkm*nx);
    frl = B*T2;
    %then generate the conditional samples
    for j = 1:dr
        k = dx+j;
        rkm = size(obj.approx.data.cores{k}, 1);
        %nc  = cardinal(obj.oned_cdfs{k});
        nn  = length(obj.oned_cdfs{k}.nodes);
        
        %bk  = compute_basis(obj.approx.base.oneds{k}, obj.oned_cdfs{k});
        T1  = reshape(TTSIRT.eval_oned_core_213(obj.approx.base.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), rkm, []);
        pk  = reshape(sum(reshape(frl*T1, nr*nn, []).^2, 2), nr, nn)';
        %
        r(j,:)  = invert_cdf(obj.oned_cdfs{k}, pk+obj.tau, z(j,:));
        %
        % evaluate the updated basis function
        T2  = TTSIRT.eval_oned_core_213(obj.approx.base.oneds{k}, obj.approx.data.cores{k}, r(j,:));
        jj  = reshape(reshape(1:rkm*nr, rkm, nr)', [], 1);
        ii  = repmat((1:nr)', 1, rkm);
        B   = sparse(ii(:), jj(:), frl(:), nr, rkm*nr);
        frl = B*T2;
        %frl(isnan(frl)) = 0;
    end
    f = frl'.^2;
    mlogw = eval_measure_potential_reference(obj.approx.base, r, (1:dr)+dx);
    f = log(fm+obj.tau) - log(f+obj.tau) + mlogw;
else
    %from the first dimension, marginalise to the last
    %order of the samples (r, x)
    %first evaluate the marginal for x
    frg = ones(1,nr);
    for j = dx:-1:2
        k   = dr + j;
        rk  = size(obj.approx.data.cores{k}, 3);
        % evaluate the updated basis function
        T2  = TTSIRT.eval_oned_core_231(obj.approx.base.oneds{k}, obj.approx.data.cores{k}, x(j,:));
        ii  = reshape(1:rk*nx, [], 1);
        jj  = reshape(repmat(1:nx, rk, 1), [], 1);
        B   = sparse(ii, jj, frg(:), rk*nx, nx);
        frg = T2'*B;
    end
    k   = dr+1;
    rk  = size(obj.approx.data.cores{k}, 3);
    % evaluate the updated basis function
    
    % for marginal
    T2  = TTSIRT.eval_oned_core_231(obj.approx.base.oneds{k}, obj.ys{k}, x(1,:));
    ii  = reshape(1:rk*nx, [], 1);
    jj  = reshape(repmat(1:nx, rk, 1), [], 1);
    B   = sparse(ii, jj, frg(:), rk*nx, nx);
    frg_m   = T2'*B;
    fm  = sum(frg_m.^2, 1);
    %
    % for frg
    T2  = TTSIRT.eval_oned_core_231(obj.approx.base.oneds{k}, obj.approx.data.cores{k}, x(1,:));
    ii  = reshape(1:rk*nx, [], 1);
    jj  = reshape(repmat(1:nx, rk, 1), [], 1);
    B   = sparse(ii, jj, frg(:), rk*nx, nx);
    frg = T2'*B;
    %then generate the conditional samples
    for k = dr:-1:1
        rk  = size(obj.approx.data.cores{k}, 3);
        %nc  = cardinal(obj.oned_cdfs{k});
        nn  = length(obj.oned_cdfs{k}.nodes);
        T1  = reshape(TTSIRT.eval_oned_core_213(obj.approx.base.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), [], rk);
        pk  = reshape(sum(reshape(T1*frg, [], nn*nr).^2, 1), nn, nr);
        %
        r(k,:)  = invert_cdf(obj.oned_cdfs{k}, pk+obj.tau, z(k,:));
        %
        % evaluate the updated basis function
        T2  = TTSIRT.eval_oned_core_231(obj.approx.base.oneds{k}, obj.approx.data.cores{k}, r(k,:));
        ii  = reshape(1:rk*nr, [], 1);
        jj  = reshape(repmat(1:nr, rk, 1), [], 1);
        B   = sparse(ii, jj, frg(:), rk*nr, nr);
        frg = T2'*B;
    end
    f = frg.^2;
    mlogw = eval_measure_potential_reference(obj.approx.base, r, dr:-1:1);
    f = log(fm+obj.tau) - log(f+obj.tau) + mlogw;
end


end

