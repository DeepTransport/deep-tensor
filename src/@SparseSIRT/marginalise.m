function obj = marginalise(obj, order)
% Marginalise the pdf represented by FTT dimension by dimension.
%   irt = MARGINALISE(irt, order)
%
%   int_dir - The direction of the marginalisation
%             >0: marginalise to the first from the last dimension 
%             <0: marginalise to the last from the first dimension 


if nargin == 1
    order = 1;
end

d = ndims(obj.approx);
nc = cardinal(obj.approx);

if numel(order) == 1
    if order > 0
        obj.order = 1:d;
    else
        obj.order = d:-1:1;
    end
elseif numel(order) == d
    obj.order = order;
else
    error('order should either be a scalar to specify the direction or a d dimensional vector')
end
    
% the multi indices: obj.est_fun.basis.indices.array; m x d
% m = size(obj.est_fun.basis.indices.array, 1);
obj.cum_P = cell(d-1, 1);
for k = 1:(d-1)
    alpha = obj.order((k+1):d);
    Kalpha = obj.approx.data.I.array(:,alpha);
    [~,~,iuni] = unique(Kalpha, 'rows');
    %
    obj.cum_P{obj.order(k)} = sparse((1:nc)',iuni(:),ones(nc,1),nc,numel(unique(iuni)) );
end
%obj.cum_P{obj.order(d)} = speye(nc,nc);

for k = 1:d
    ind = obj.approx.data.I.array(:,k);
    % evaluate polynomial basis on all cdf nodes
    bk = eval_basis(obj.approx.base.oneds{k}, obj.oned_cdfs{k}.nodes(:));
    obj.bases_cdf_nodes{k} = bk(:,ind+1); % cdf nodes for each indices, nc x nk
end

obj.fun_z = sum(obj.approx.data.coeff.^2);
obj.z = obj.fun_z + obj.tau;

end
