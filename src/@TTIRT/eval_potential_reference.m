function fx = eval_potential_reference(obj, x)
% Evaluate the normalised (marginal) pdf represented by squared FTT.
%   f = EVAL_POTENTIAL_REFERENCE(irt, x)
%
%   x - input variables in the reference coordinates
%   f - potential function at x

% map to the reference coordinate

d = ndims(obj.approx);
[dx,n] = size(x);
fradon = ones(1,n);

if obj.int_dir > 0 % from left to right
    frl = ones(n,1);
    for k = 1:dx
        rkm = size(obj.approx.data.cores{k}, 1);
        %
        pk  = abs(frl*obj.ys{k})';
        %
        tmp = eval_pdf_radon(obj.oned_cdfs{k}, pk+obj.tau, x(k,:));
        fradon = fradon.*reshape(tmp,size(fradon));
        %
        % evaluate the updated basis function
        T2  = TTSIRT.eval_oned_core_213(obj.approx.base.oneds{k}, obj.approx.data.cores{k}, x(k,:));
        jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
        ii  = repmat((1:n)', 1, rkm);
        B   = sparse(ii(:), jj(:), frl(:), n, rkm*n);
        frl = B*T2;
    end
    mlogw = eval_measure_potential_reference(obj.approx.base, x, 1:dx);
    fx = - log(fradon+obj.tau) + mlogw;
else % from right to left
    ie = (d-dx)+1;
    frg = ones(1,n);
    for k = d:-1:ie
        ck  = k-ie+1;
        rk  = size(obj.approx.data.cores{k}, 3);
        %
        pk  = abs(obj.ys{k}*frg);
        %
        tmp = eval_pdf_radon(obj.oned_cdfs{k}, pk+obj.tau, x(ck,:));
        fradon = fradon.*reshape(tmp,size(fradon));
        %
        % evaluate the updated basis function
        T2  = TTSIRT.eval_oned_core_231(obj.approx.base.oneds{k}, obj.approx.data.cores{k}, x(ck,:));
        ii  = reshape(1:rk*n, [], 1);
        jj  = reshape(repmat(1:n, rk, 1), [], 1);
        B   = sparse(ii, jj, frg(:), rk*n, n);
        frg = T2'*B;
    end
    mlogw = eval_measure_potential_reference(obj.approx.base, x, d:-1:ie);
    fx = - log(fradon+obj.tau) + mlogw;
end


end

