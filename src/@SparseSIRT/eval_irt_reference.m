function [r,f,g] = eval_irt_reference(obj, z)
% Evaluate the inverse of the squared Rosenblatt transport X = R^{-1}(Z),
% where X is the target random variable and Z is uniform.
%   [X,f,g] = EVAL_IRT(irt, Z)
%
%   Z - uniform random variables, d x n
%   X - random variable drawn form the pdf defined by SIRT
%   f - pdf of X
%   g - gradient of log pdf at X

d = ndims(obj.approx);
nk = cardinal(obj.approx);
[dz,n] = size(z);
r = zeros(dz,n);

nbatch = ceil(n/obj.batch_size);
start_i = 1:obj.batch_size:n;
end_i = start_i+(obj.batch_size-1);
end_i(end) = n;

if nargout < 3
    brl = repmat(obj.approx.data(:)', n, 1); % n x nk
    for k = 1:dz
        nc = cardinal(obj.oned_cdfs{obj.order(k)});
        % for each i in (1..nk), outter product of bc(:,i) and brl(:,i)
        % tmp is nc x (n x nk)
        for batch = 1:nbatch
            ind = start_i(batch):1:end_i(batch);
            ni = numel(ind);
            tmp = reshape(repmat(obj.bases_cdf_nodes{obj.order(k)},ni,1), nc, ni*nk).*reshape(brl(ind,:),1,ni*nk);
            if k < d
                pk = reshape( sum((reshape(tmp,nc*ni,nk)*obj.cum_P{obj.order(k)}).^2,2), nc, ni);
            else
                pk = reshape( sum(reshape(tmp,nc*ni,nk),2).^2, nc, ni);
            end
            r(obj.order(k),ind) = invert_cdf(obj.oned_cdfs{obj.order(k)}, pk+obj.tau, z(obj.order(k),ind));
        end
        %
        bk = eval_basis(obj.approx.oneds{obj.order(k)}, r(obj.order(k),:));
        brl = brl.*bk(:,obj.approx.indices.array(:,obj.order(k))+1);
    end
    if dz < d
        f = sum((reshape(brl,n,nk)*obj.cum_P{obj.order(dz+1)}).^2,2);
    else
        f = sum(reshape(brl,n,nk),2).^2;
    end
    mlogw = eval_measure_potential_reference(obj.approx, r, 1:dz);
    f = log(obj.z) - log(f(:)'+obj.tau) + mlogw;
else
    if dz == d
        g = zeros(d,n);
        fls = cell(d,1);
        frs = cell(d,1);
        block_basis = cell(d,1);
        %
        brl = repmat(obj.approx.data(:)', n, 1); % n x nk
        for k = 1:d
            nc = cardinal(obj.oned_cdfs{obj.order(k)});
            % for each i in (1..nk), outter product of bc(:,i) and brl(:,i)
            % tmp is nc x (n x nk)
            for batch = 1:nbatch
                ind = start_i(batch):1:end_i(batch);
                ni = numel(ind);
                tmp = reshape(repmat(obj.bases_cdf_nodes{obj.order(k)},ni,1), nc, ni*nk).*reshape(brl(ind,:),1,ni*nk);
                if k < d
                    pk = reshape( sum((reshape(tmp,nc*ni,nk)*obj.cum_P{obj.order(k)}).^2,2), nc, ni);
                else
                    pk = reshape( sum(reshape(tmp,nc*ni,nk),2).^2, nc, ni);
                end
                r(obj.order(k),ind) = invert_cdf(obj.oned_cdfs{obj.order(k)}, pk+obj.tau, z(obj.order(k),ind));
            end
            B = eval_basis(obj.approx.oneds{obj.order(k)}, r(obj.order(k),:));
            block_basis{obj.order(k)} = B(:,obj.approx.indices.array(:,obj.order(k))+1);
            %
            brl = brl.*block_basis{obj.order(k)};
            fls{obj.order(k)} = brl;
        end
        fsqrt = sum(reshape(brl,n,nk),2)';
        %
        frs{obj.order(d)} = block_basis{obj.order(d)};
        for k = (d-1):-1:2
            frs{obj.order(k)} = frs{obj.order(k+1)}.*block_basis{obj.order(k)};
        end
        % assemble gradient
        for k = 1:d
            dB = eval_basis_deri(obj.approx.oneds{obj.order(k)}, r(obj.order(k),:));
            D = dB(:,obj.approx.indices.array(:,obj.order(k))+1);
            if k == 1
                g(obj.order(k),:) = ( (D.*frs{obj.order(k+1)})*obj.approx.data(:) )';
            elseif k == d
                g(obj.order(k),:) = sum(D.*fls{obj.order(k-1)}, 2)';
            else
                g(obj.order(k),:) = sum(D.*fls{obj.order(k-1)}.*frs{obj.order(k+1)}, 2)';
            end
        end
        %
        mlogw = eval_measure_potential_reference(obj.approx, r);
        gmlogw = eval_measure_potential_reference_grad(obj.approx, r);
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
