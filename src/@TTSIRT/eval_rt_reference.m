function z = eval_rt_reference(obj, r)
% Evaluate the squared Rosenblatt transport Z = R(X), where Z is the uniform
% random variable and X is the target random variable.
%   Z = EVAL_RT_REFERENCE(irt, X)
%
%   Z - uniform random variables, d x n
%   X - random variable drawn form the pdf defined by SIRT

d = ndims(obj.approx);
[dr,n] = size(r);
z = zeros(dr,n);

if obj.int_dir > 0 % from left to right
    frl = ones(n,1);
    for k = 1:dr
        rkm = size(obj.approx.cores{k}, 1);
        %nc  = cardinal(obj.oned_cdfs{k});
        nn  = length(obj.oned_cdfs{k}.nodes);
        %
        T1  = reshape( TTSIRT.eval_oned_core_213(obj.approx.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), rkm, []);
        pk  = reshape(sum(reshape(frl*T1, n*nn, []).^2, 2), n, nn)'; % squared TT
        %
        % OLD: z(k,:)  = eval_cdf(obj.oned_cdfs{k}, pk, obj.tau*fr_ref, obj.oned_refs{k}, r(k,:));
        z(k,:) = eval_cdf(obj.oned_cdfs{k}, pk+obj.tau, r(k,:));
        %
        % evaluate the updated basis function
        T2  = TTSIRT.eval_oned_core_213(obj.approx.oneds{k}, obj.approx.cores{k}, r(k,:));
        jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
        ii  = repmat((1:n)', 1, rkm);
        B   = sparse(ii(:), jj(:), frl(:), n, rkm*n);
        frl = B*T2;
        %frl(isnan(frl)) = 0;
    end
else % from right to left
    ie = (d-dr)+1;
    frg = ones(1,n);
    for k = d:-1:ie
        ck  = k-ie+1;
        rk  = size(obj.approx.cores{k}, 3);
        %nc  = cardinal(obj.oned_cdfs{k});
        nn  = length(obj.oned_cdfs{k}.nodes);
        T1  = reshape( TTSIRT.eval_oned_core_213(obj.approx.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), [], rk);
        pk  = reshape(sum(reshape(T1*frg, [], nn*n).^2, 1), nn, n);
        %
        % OLD: z(ck,:)  = eval_cdf(obj.oned_cdfs{k}, pk, obj.tau*fr_ref, obj.oned_refs{k}, r(ck,:));
        z(ck,:)  = eval_cdf(obj.oned_cdfs{k}, pk+obj.tau, r(ck,:));
        %
        % evaluate the updated basis function
        T2  = TTSIRT.eval_oned_core_231(obj.approx.oneds{k}, obj.approx.cores{k}, r(ck,:));
        ii  = reshape(1:rk*n, [], 1);
        jj  = reshape(repmat(1:n, rk, 1), [], 1);
        B   = sparse(ii, jj, frg(:), rk*n, n);
        frg = T2'*B;
        %frg(isnan(frg)) = 0;
    end
end

end

%{
function z = eval_rt_reference(obj, r)
% Evaluate the squared Rosenblatt transport Z = R(X), where Z is the uniform
% random variable and X is the target random variable.
%   Z = EVAL_RT_REFERENCE(irt, X)
%
%   Z - uniform random variables, d x n
%   X - random variable drawn form the pdf defined by SIRT

d = ndims(obj.approx);
[dr,n] = size(r);
z = zeros(dr,n);

if obj.int_dir > 0 % from left to right
    frl = ones(n,1);
    fr_ref = ones(1,n);
    for k = 1:dr
        rkm = size(obj.approx.cores{k}, 1);
        nc  = cardinal(obj.oned_cdfs{k});
        %
        T1  = reshape( TTSIRT.eval_oned_core_213(obj.approx.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), rkm, []);
        pk  = reshape(sum(reshape(frl*T1, n*nc, []).^2, 2), n, nc)'; % squared TT
        %
        % OLD: z(k,:)  = eval_cdf(obj.oned_cdfs{k}, pk, obj.tau*fr_ref, obj.oned_refs{k}, r(k,:));
        z(k,:) = eval_cdf(obj.oned_cdfs{k}, pk+obj.tau, r(k,:));
        %
        % evaluate the updated basis function
        T2  = TTSIRT.eval_oned_core_213(obj.approx.oneds{k}, obj.approx.cores{k}, r(k,:));
        jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
        ii  = repmat((1:n)', 1, rkm);
        B   = sparse(ii(:), jj(:), frl(:), n, rkm*n);
        frl = B*T2;
        %frl(isnan(frl)) = 0;
        %
        fk_ref = eval_pdf(obj.oned_refs{k}, r(k,:));
        fr_ref = fr_ref.*fk_ref;
    end
else % from right to left
    ie = (d-dr)+1;
    frg = ones(1,n);
    fr_ref = ones(1,n);
    for k = d:-1:ie
        ck  = k-ie+1;
        rk  = size(obj.approx.cores{k}, 3);
        nc  = cardinal(obj.oned_cdfs{k});
        T1  = reshape( TTSIRT.eval_oned_core_213(obj.approx.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), [], rk);
        pk  = reshape(sum(reshape(T1*frg, [], nc*n).^2, 1), nc, n);
        %
        z(ck,:)  = eval_cdf(obj.oned_cdfs{k}, pk, obj.tau*fr_ref, obj.oned_refs{k}, r(ck,:));
        %
        % evaluate the updated basis function
        T2  = TTSIRT.eval_oned_core_231(obj.approx.oneds{k}, obj.approx.cores{k}, r(ck,:));
        ii  = reshape(1:rk*n, [], 1);
        jj  = reshape(repmat(1:n, rk, 1), [], 1);
        B   = sparse(ii, jj, frg(:), rk*n, n);
        frg = T2'*B;
        %frg(isnan(frg)) = 0;
        %
        fk_ref = eval_pdf(obj.oned_refs{k}, r(ck,:));
        fr_ref = fr_ref.*fk_ref;
    end
end

end
%}