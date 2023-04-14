function J = eval_rt_jac_reference(obj, r, z)
% Evaluate the jacobian of the squared Rosenblatt transport Z = R(X), where 
% Z is the uniform random variable and X is the target random variable. 
%   J = EVAL_RT_JAC_REFERENCE(irt, X, Z)
%
%   X - random variable drawn form the pdf defined by SIRT
%   Z - uniform random variables, d x n
%   J - Jacobian, d x (d x n), each d x d block is the Jabocian for X(:,j)

d = ndims(obj.approx);
n = size(r,2);
if size(r,1) ~= d
    error('Input variable must be d x n')
end
J = zeros(d,n*d);

if obj.int_dir > 0 % from left to right
    % function value of each block
    block_ftt = cell(d,1);
    block_mar = cell(d,1);
    mlogw = cell(d,1);
    T = cell(d,1);
    block_ftt_d = cell(d,1);
    for k = 1:d
        block_ftt{k} = TTSIRT.eval_oned_core_213(obj.approx.oneds{k}, obj.approx.cores{k}, r(k,:));
        block_mar{k} = TTSIRT.eval_oned_core_213(obj.approx.oneds{k}, obj.ys{k}, r(k,:));
        %
        rkm  = size(obj.approx.cores{k}, 1);
        T{k} = reshape(TTSIRT.eval_oned_core_213(obj.approx.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), rkm, []);
        %
        block_ftt_d{k} = TTSIRT.eval_oned_core_213_deri(obj.approx.oneds{k}, obj.approx.cores{k}, r(k,:));
        %
        logw = eval_log_measure(obj.approx.oneds{k}, r(k,:));
        mlogw{k} = -logw;
    end
    %
    F = cell(d,1); % accumulated ftt
    G = cell(d,1); % sum( G^.2 ) is the marginal, each G{k} is nxr
    F{1} = block_ftt{1};
    G{1} = block_mar{1};
    for k = 2:d
        rkm = size(obj.approx.cores{k}, 1);
        jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
        ii  = repmat((1:n)', 1, rkm);
        B   = sparse(ii(:), jj(:), F{k-1}(:), n, rkm*n);
        F{k} = B*block_ftt{k};
        G{k} = B*block_mar{k};
    end
    %
    Fm = cell(d,1); % accumulated ftt
    for k = 1:d
        Fm{k} = sum(G{k}.^2, 2)' + obj.tau;
    end
    %
    % j is the index of the coordinate of differentiation
    % d(T_k)/d(x_j), k:row, j: col
    for j = 1:d
        ind = (1:d:n*d) + (j-1);
        %
        % the diagonal
        if j == 1
            J(1,ind) = Fm{1}/obj.z;
        else
            J(j,ind) = Fm{j}./Fm{j-1};
        end
        J(j,ind) = J(j,ind).*exp(-mlogw{j});
        %
        if j < d % skip the (d,d) element
            % derivative of the FTT
            drl = block_ftt_d{j};
            %
            % derivative of the FTT, 2nd term, for the d(j+1)/dj term
            mrl = TTSIRT.eval_oned_core_213_deri(obj.approx.oneds{j}, obj.ys{j}, r(j,:));
            if j > 1
                rjm = size(obj.approx.cores{j}, 1);
                jj  = reshape(reshape(1:rjm*n, rjm, n)', [], 1);
                ii  = repmat((1:n)', 1, rjm);
                B   = sparse(ii(:), jj(:), F{j-1}(:), n, rjm*n);
                drl = B*drl;
                mrl = B*mrl;
            end
            % first sub, the second term, for the d(j+1)/dj term
            J(j+1,ind) = J(j+1,ind) - 2*sum(G{j}.*mrl,2)'.*z(j+1,:);
            %
            for k = (j+1):d 
                % accumulate the j-th block and ealuate the integral
                rkm = size(obj.approx.cores{k}, 1);
                %nc  = cardinal(obj.oned_cdfs{k});
                nn  = length(obj.oned_cdfs{k}.nodes);
                pk  = reshape(sum( reshape((F{k-1}*T{k}).*(drl*T{k}), n*nn, []), 2), n, nn)';
                % the first term
                if obj.approx.oneds{k}.constant_weight
                    tmp = eval_measure(obj.approx.oneds{k}, obj.oned_cdfs{k}.nodes);
                    pk = pk.*tmp(:);
                end
                J(k,ind) = J(k,ind) + 2*reshape(eval_int_deri(obj.oned_cdfs{k}, pk, r(k,:)), 1, []);
                %
                if k < d
                    %
                    jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
                    ii  = repmat((1:n)', 1, rkm);
                    B   = sparse(ii(:), jj(:), drl(:), n, rkm*n);
                    % acumulate
                    drl = B*block_ftt{k};
                    % the second term, for the d(k+1)/dj term
                    mrl = B*block_mar{k};
                    J(k+1,ind) = J(k+1,ind) - 2*sum(G{k}.*mrl,2)'.*z(k+1,:);
                end
                J(k,ind) = J(k,ind)./Fm{k-1};
            end
        end
        %
    end
    
else % from right to left
    % function value of each block
    block_ftt = cell(d,1);
    block_mar = cell(d,1);
    mlogw = cell(d,1);
    T = cell(d,1);
    block_ftt_d = cell(d,1);
    for k = 1:d
        block_ftt{k} = TTSIRT.eval_oned_core_231(obj.approx.oneds{k}, obj.approx.cores{k}, r(k,:));
        block_mar{k} = TTSIRT.eval_oned_core_231(obj.approx.oneds{k}, obj.ys{k}, r(k,:));
        rk   = size(obj.approx.cores{k}, 3);
        T{k} = reshape(TTSIRT.eval_oned_core_213(obj.approx.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), [], rk);
        block_ftt_d{k} = TTSIRT.eval_oned_core_231_deri(obj.approx.oneds{k}, obj.approx.cores{k}, r(k,:));
        %
        logw = eval_log_measure(obj.approx.oneds{k}, r(k,:));
        mlogw{k} = -logw;
    end
    %
    F = cell(d,1); % accumulated ftt
    G = cell(d,1); % sum( G^.2 ) is the marginal, each G{k} is nxr
    F{d} = block_ftt{d}';
    G{d} = block_mar{d}';
    for k = d-1:-1:1
        rk  = size(obj.approx.cores{k}, 3);
        ii  = reshape(1:rk*n, [], 1);
        jj  = reshape(repmat(1:n, rk, 1), [], 1);
        B   = sparse(ii, jj, F{k+1}(:), rk*n, n);
        F{k} = block_ftt{k}'*B;
        G{k} = block_mar{k}'*B;
    end
    %
    Fm = cell(d,1); % accumulated ftt
    for k = 1:d
        Fm{k} = sum(G{k}.^2, 1) + obj.tau;
    end
    %
    % j is the index of the coordinate of differentiation
    for j = d:-1:1
        ind = (1:d:n*d) + (j-1);
        %
        % the diagonal
        if j == d
            J(d,ind) = Fm{d}/obj.z;
        else
            J(j,ind) = Fm{j}./Fm{j+1};
        end
        J(j,ind) = J(j,ind).*exp(-mlogw{j});
        %
        if j > 1 % skip the (1,1) element
            % derivative of the FTT
            drg = block_ftt_d{j}';
            % derivative of the FTT, 2nd term, for the d(j+1)/dj term
            mrg = TTSIRT.eval_oned_core_231_deri(obj.approx.oneds{j}, obj.ys{j}, r(j,:))';
            if j < d
                rj  = size(obj.approx.cores{j}, 3);
                ii  = reshape(1:rj*n, [], 1);
                jj  = reshape(repmat(1:n, rj, 1), [], 1);
                B   = sparse(ii, jj, F{j+1}(:), rj*n, n);
                drg = drg*B;
                mrg = mrg*B;
            end
            % first super, the second term, for the d(j-1)/dj term
            J(j-1,ind) = J(j-1,ind) - 2*sum(G{j}.*mrg,1).*z(j-1,:);
            %
            for k = (j-1):-1:1
                % accumulate the j-th block and evaluate the integral
                rk  = size(obj.approx.cores{k}, 3);
                %nc  = cardinal(obj.oned_cdfs{k});
                nn  = length(obj.oned_cdfs{k}.nodes);
                pk  = reshape(sum(reshape((T{k}*F{k+1}).*(T{k}*drg), [], nn*n), 1), nn, n);
                if obj.approx.oneds{k}.constant_weight
                    tmp = eval_measure(obj.approx.oneds{k}, obj.oned_cdfs{k}.nodes);
                    pk = pk.*tmp(:);
                end
                % the first term
                J(k,ind) = J(k,ind) + 2*reshape(eval_int_deri(obj.oned_cdfs{k}, pk, r(k,:)), 1, []);
                %
                if k > 1
                    ii  = reshape(1:rk*n, [], 1);
                    jj  = reshape(repmat(1:n, rk, 1), [], 1);
                    B   = sparse(ii, jj, drg(:), rk*n, n);
                    % acumulate
                    drg = block_ftt{k}'*B;
                    % the second term, for the d(k+1)/dj term
                    mrg = block_mar{k}'*B;
                    J(k-1,ind) = J(k-1,ind) - 2*sum(G{k}.*mrg,1).*z(k-1,:);
                end
                J(k,ind) = J(k,ind)./Fm{k+1};
            end
        end
    end
end

end
