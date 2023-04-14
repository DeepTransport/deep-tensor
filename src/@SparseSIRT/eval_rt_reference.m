function z = eval_rt_reference(obj, r)
% Evaluate the squared Rosenblatt transport Z = R(X), where Z is the uniform
% random variable and X is the target random variable.
%   Z = EVAL_RT(irt, X)
%
%   Z - uniform random variables, d x n
%   X - random variable drawn form the pdf defined by SIRT

d = ndims(obj.approx);
nk = cardinal(obj.approx);
[dr,n] = size(r);
z = zeros(dr,n);

nbatch = ceil(n/obj.batch_size);
start_i = 1:obj.batch_size:n;
end_i = start_i+(obj.batch_size-1);
end_i(end) = n;

brl = repmat(obj.approx.data(:)', n, 1); % n x nk
for k = 1:dr
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
        z(obj.order(k),ind) = eval_cdf(obj.oned_cdfs{obj.order(k)}, pk+obj.tau, r(obj.order(k),ind));
    end
    %
    bk = eval_basis(obj.approx.oneds{obj.order(k)}, r(obj.order(k),:));
    brl = brl.*bk(:,obj.approx.indices.array(:,obj.order(k))+1);
end

end
