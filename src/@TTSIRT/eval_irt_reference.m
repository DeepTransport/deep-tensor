function [r,f,g] = eval_irt_reference(obj, z)
% Evaluate the inverse of the squared Rosenblatt transport X = R^{-1}(Z),
% where X is the target random variable and Z is uniform.
%   [X,f,g] = EVAL_IRT_REFERENCE(irt, Z)
%
%   Z - uniform random variables, d x n
%   X - random variable drawn form the pdf defined by SIRT
%   f - potential function at X
%   g - gradient of potential function at X

d = ndims(obj.approx);
[dz,n] = size(z);
r = zeros(dz,n);

if nargout < 3
    if obj.int_dir > 0 % from left to right
        frl = ones(n,1);
        for k = 1:dz
            rkm = size(obj.approx.data.cores{k}, 1);
            %nc  = cardinal(obj.oned_cdfs{k});
            nn  = length(obj.oned_cdfs{k}.nodes);
            T1  = reshape( TTSIRT.eval_oned_core_213(obj.approx.base.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), rkm, []);
            pk  = reshape(sum(reshape(frl*T1, n*nn, []).^2, 2), n, nn)';
            %
            %OLD: r(k,:)  = invert_cdf(obj.oned_cdfs{k}, pk, obj.tau*fr_ref, obj.oned_refs{k}, z(k,:));
            r(k,:)  = invert_cdf(obj.oned_cdfs{k}, pk+obj.tau, z(k,:));
            %
            % evaluate the updated basis function
            T2  = TTSIRT.eval_oned_core_213(obj.approx.base.oneds{k}, obj.approx.data.cores{k}, r(k,:));
            jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
            ii  = repmat((1:n)', 1, rkm);
            B   = sparse(ii(:), jj(:), frl(:), n, rkm*n);
            frl = B*T2;
            frl(isnan(frl)) = 0;
        end
        if dz < d
            f = sum((frl*obj.ms{dz+1}).^2, 2)';
        else
            f = frl'.^2;
        end
        mlogw = eval_measure_potential_reference(obj.approx.base, r, 1:dz);
        f = log(obj.z) - log(f+obj.tau) + mlogw;
    else % from right to left
        frg = ones(1,n);
        ie  = (d-dz)+1;
        for k = d:-1:ie
            ck  = k-ie+1;
            rk  = size(obj.approx.data.cores{k}, 3);
            %nc  = cardinal(obj.oned_cdfs{k});
            nn  = length(obj.oned_cdfs{k}.nodes);
            T1  = reshape( TTSIRT.eval_oned_core_213(obj.approx.base.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), [], rk);
            pk  = reshape(sum(reshape(T1*frg, [], nn*n).^2, 1), nn, n);
            %
            %OLD: r(ck,:)  = invert_cdf(obj.oned_cdfs{k}, pk, obj.tau*fr_ref, obj.oned_refs{k}, z(ck,:));
            r(ck,:)  = invert_cdf(obj.oned_cdfs{k}, pk+obj.tau, z(ck,:));
            %
            % evaluate the updated basis function
            T2  = TTSIRT.eval_oned_core_231(obj.approx.base.oneds{k}, obj.approx.data.cores{k}, r(ck,:));
            ii  = reshape(1:rk*n, [], 1);
            jj  = reshape(repmat(1:n, rk, 1), [], 1);
            B   = sparse(ii, jj, frg(:), rk*n, n);
            frg = T2'*B;
            frg(isnan(frg)) = 0;
        end
        if dz < d
            f = sum((obj.ms{ie-1}*frg).^2, 1);
        else
            f = frg.^2;
        end
        mlogw = eval_measure_potential_reference(obj.approx.base, r, d:-1:ie);
        f = log(obj.z) - log(f+obj.tau) + mlogw;
    end
else 
    if dz == d
        fls = cell(d,1);
        frs = cell(d,1);
        if obj.int_dir > 0 % from left to right
            %
            frl = ones(n,1);
            for k = 1:dz
                rkm = size(obj.approx.data.cores{k}, 1);
                nc  = cardinal(obj.oned_cdfs{k});
                %nn  = length(obj.oned_cdfs{k}.nodes);
                T1  = reshape( TTSIRT.eval_oned_core_213(obj.approx.base.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), rkm, []);
                pk  = reshape(sum(reshape(frl*T1, n*nc, []).^2, 2), n, nc)';
                %
                r(k,:)  = invert_cdf(obj.oned_cdfs{k}, pk+obj.tau, z(k,:));
                %
                % evaluate the updated basis function
                T2  = TTSIRT.eval_oned_core_213(obj.approx.base.oneds{k}, obj.approx.data.cores{k}, r(k,:));
                jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
                ii  = repmat((1:n)', 1, rkm);
                B   = sparse(ii(:), jj(:), frl(:), n, rkm*n);
                frl = B*T2;
                frl(isnan(frl)) = 0;
                fls{k} = frl;
            end
            %
            frg = ones(1,n);
            for k = dz:-1:1
                rk  = size(obj.approx.data.cores{k}, 3);
                %
                % evaluate the updated basis function
                T2  = TTSIRT.eval_oned_core_231(obj.approx.base.oneds{k}, obj.approx.data.cores{k}, r(k,:));
                ii  = reshape(1:rk*n, [], 1);
                jj  = reshape(repmat(1:n, rk, 1), [], 1);
                B   = sparse(ii, jj, frg(:), rk*n, n);
                frg = T2'*B;
                frs{k} = frg;
            end
            %
            fsqrt = frl';
        else % from right to left
            frg = ones(1,n);
            for k = dz:-1:1
                rk  = size(obj.approx.data.cores{k}, 3);
                nc  = cardinal(obj.oned_cdfs{k});
                %nn  = length(obj.oned_cdfs{k}.nodes);
                T1  = reshape( TTSIRT.eval_oned_core_213(obj.approx.base.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), [], rk);
                pk  = reshape(sum(reshape(T1*frg, [], nc*n).^2, 1), nc, n);
                %
                r(k,:)  = invert_cdf(obj.oned_cdfs{k}, pk+obj.tau, z(k,:));
                %
                % evaluate the updated basis function
                T2  = TTSIRT.eval_oned_core_231(obj.approx.base.oneds{k}, obj.approx.data.cores{k}, r(k,:));
                ii  = reshape(1:rk*n, [], 1);
                jj  = reshape(repmat(1:n, rk, 1), [], 1);
                B   = sparse(ii, jj, frg(:), rk*n, n);
                frg = T2'*B;
                frg(isnan(frg)) = 0;
                frs{k} = frg;
            end
            %
            frl = ones(n,1);
            for k = 1:dz
                rkm = size(obj.approx.data.cores{k}, 1);
                % evaluate the updated basis function
                T2  = TTSIRT.eval_oned_core_213(obj.approx.base.oneds{k}, obj.approx.data.cores{k}, r(k,:));
                jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
                ii  = repmat((1:n)', 1, rkm);
                B   = sparse(ii(:), jj(:), frl(:), n, rkm*n);
                frl = B*T2;
                fls{k} = frl;
            end
            %
            fsqrt = frg;
        end
        % assemble gradient
        g = zeros(dz,n);
        for k = 1:dz
            rkm = size(obj.approx.data.cores{k}, 1);
            %
            % (rjm nx) by rj
            D = TTSIRT.eval_oned_core_213_deri(obj.approx.base.oneds{k}, obj.approx.data.cores{k}, r(k,:));
            if k == 1
                g(k,:) = sum(D'.*frs{k+1}, 1);
            elseif k == d
                jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
                ii  = repmat((1:n)', 1, rkm);
                Bl  = sparse(ii(:), jj(:), fls{k-1}(:), n, rkm*n);
                g(k,:) =(Bl*D)';
            else
                jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
                ii  = repmat((1:n)', 1, rkm);
                Bl  = sparse(ii(:), jj(:), fls{k-1}(:), n, rkm*n);
                g(k,:) = sum((Bl*D)'.*frs{k+1}, 1);
            end
        end
        %
        mlogw = eval_measure_potential_reference(obj.approx.base, r);
        gmlogw = eval_measure_potential_reference_grad(obj.approx.base, r);
        %
        f = obj.tau + fsqrt.^2;
        %
        g = 2*(g.*fsqrt);
        g = -g./f + gmlogw; % -g/f
        g(isnan(g)) = 0;
        g(isinf(g)) = 0;
        %
        % normalize at the end
        f = log(obj.z) - log(f) + mlogw;
    else
        warning('grad is not implemented for marginals')
    end
end


end
