function J = eval_rt_jac_reference_inner(obj, r, z)
% Evaluate the jacobian of the squared Rosenblatt transport Z = R(X), where
% Z is the uniform random variable and X is the target random variable.
%   J = EVAL_RT_JAC_REFERENCE(irt, X, Z)
%
%   X - random variable drawn form the pdf defined by SIRT
%   Z - uniform random variables, d x n
%   J - Jacobian, d x (d x n), each d x d block is the Jabocian for X(:,j)

d = ndims(obj.approx);
n = size(r,2);
nk = cardinal(obj.approx);
if size(r,1) ~= d
    error('Input variable must be d x n')
end
J = zeros(d,n*d);

block_basis = cell(d,1);
block_basis_deri = cell(d,1);
block_mlogw = cell(d,1);
for k = 1:d
    B = eval_basis(obj.approx.base.oneds{k}, r(k,:));
    block_basis{k} = B(:,obj.approx.data.I.array(:,k)+1);
    dB = eval_basis_deri(obj.approx.base.oneds{k}, r(k,:));
    block_basis_deri{k} = dB(:,obj.approx.data.I.array(:,k)+1);
    logw = eval_log_measure(obj.approx.base.oneds{k}, r(k,:));
    block_mlogw{k} = -logw(:);
end

F = cell(d,1); % accumulated basis function, with coefficients
G = cell(d,1); % sum( G^.2 ) is the marginal
pk_left = cell(d,1);
Fm = cell(d,1); % accumulated ftt

%
F{obj.order(1)} = block_basis{obj.order(1)}.*obj.approx.data.coeff(:)';
for k = 2:d
    F{obj.order(k)} = F{obj.order(k-1)}.*block_basis{obj.order(k)};
end
for k = 1:d-1
    G{obj.order(k)} = F{obj.order(k)}*obj.cum_P{obj.order(k)};
    Fm{obj.order(k)} = sum(G{obj.order(k)}.^2, 2) + obj.tau;
end
G{obj.order(d)} = sum(F{obj.order(d)},2);
Fm{obj.order(d)} = G{obj.order(d)}.^2 + obj.tau;
%
for k = 2:d
    nc = cardinal(obj.oned_cdfs{obj.order(k)});
    tmp = reshape(repmat(obj.bases_cdf_nodes{obj.order(k)},n,1), nc, n*nk).*reshape(F{obj.order(k-1)},1,n*nk);
    if k < d
        pk_left{obj.order(k)} = reshape(tmp,nc*n,nk)*obj.cum_P{obj.order(k)};
    else
        pk_left{obj.order(k)} = sum(reshape(tmp,nc*n,nk),2);
    end
end
%
% j is the index of the coordinate of differentiation
% d(T_k)/d(x_j), k:row, j: col
for j = 1:d
    ind = (1:d:n*d) + (obj.order(j)-1); %dzdx(ind)
    %
    % the diagonal
    if j == 1
        J(obj.order(j),ind) = Fm{obj.order(j)}/obj.z;
    else
        J(obj.order(j),ind) = Fm{obj.order(j)}./Fm{obj.order(j-1)};
    end
    J(obj.order(j),ind) = J(obj.order(j),ind).*exp(-block_mlogw{obj.order(j)}');
    %
    if j < d % skip the (d,d) element
        % derivative of the FTT
        drl = block_basis_deri{obj.order(j)};
        if j > 1
            drl = F{obj.order(j-1)}.*drl;
        else
            drl = drl.*obj.approx.data.coeff(:)';
        end
        mrl = drl*obj.cum_P{obj.order(j)};
        % first sub, the second term, for the d(j+1)/dj term
        J(obj.order(j+1),ind) = J(obj.order(j+1),ind) - 2*sum(G{obj.order(j)}.*mrl,2)'.*z(obj.order(j+1),:);
        %
        for k = (j+1):d
            nc = cardinal(obj.oned_cdfs{obj.order(k)});
            tmp = reshape(repmat(obj.bases_cdf_nodes{obj.order(k)},n,1), nc, n*nk).*reshape(drl,1,n*nk);
            if k < d
                pk = reshape( sum(pk_left{obj.order(k)}.*(reshape(tmp,nc*n,nk)*obj.cum_P{obj.order(k)}),2), nc, n);
            else
                pk = reshape( pk_left{obj.order(k)}.*sum(reshape(tmp,nc*n,nk),2), nc, n);
            end
            if obj.approx.base.oneds{obj.order(k)}.constant_weight
                tmp = eval_measure(obj.approx.base.oneds{obj.order(k)}, obj.oned_cdfs{obj.order(k)}.nodes);
                pk = pk.*tmp(:);
            end
            %
            J(obj.order(k),ind) = J(obj.order(k),ind) + 2*reshape(eval_int_deri(obj.oned_cdfs{obj.order(k)}, pk, r(obj.order(k),:)), 1, []);
            %
            if k < d
                drl = drl.*block_basis{obj.order(k)};
                % the second term, for the d(k+1)/dj term
                mrl = drl*obj.cum_P{obj.order(k)};
                J(obj.order(k+1),ind) = J(obj.order(k+1),ind) - 2*sum(G{obj.order(k)}.*mrl,2)'.*z(obj.order(k+1),:);
            end
            J(obj.order(k),ind) = J(obj.order(k),ind)./Fm{obj.order(k-1)}';
        end
    end
    
    %
end

end
