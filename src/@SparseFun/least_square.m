function obj = least_square(obj,func)

% Implement Alg 4 of "Adaptive approximation by optimal weighted least 
% squares methods". Giovanni Migliorati
% The code is modified based on the AdaptiveSparseTensorAlgorithm of the
% ApproximationToolbox

d = ndims(obj.base);
max_index = MultiIndices(reshape(cardinals(obj.base)-1,1,d));

% sample_factor: average number of samples per basis function
sample_factor = ceil(obj.opt.init_sample_size);
factor_add = ceil(obj.opt.enrich_sample_size);
basis_adaptation = true;
force_enrich_sample = false;

if cardinal(obj.data.I) == 0
    switch obj.opt.indexset
        case {'hyperbolic'}
            obj.data.I = SparseTools.hyperbolic_cross(d, obj.opt.init_total_degree);
        otherwise
            obj.data.I = SparseTools.total_degree(d, obj.opt.init_total_degree);
    end
    rem = find((~(obj.data.I<=max_index)));
    obj.data.I = removeIndices(obj.data.I,rem);
    obj.data.x = random(obj, obj.data.I, sample_factor);
    obj.data.y = func(obj.data.x);
    obj.data.y = obj.data.y';
    %
    obj.data.A = eval_basis (obj, obj.data.I, obj.data.x);
    obj.data.w = eval_weight(obj, obj.data.I, obj.data.x);
else
    if isempty(obj.data.x)
        obj.data.x = random(obj, obj.data.I, sample_factor);
        obj.data.y = func(obj.data.x);
        obj.data.y = obj.data.y';
    end
    obj.data.A = eval_basis (obj, obj.data.I, obj.data.x);
    obj.data.w = eval_weight(obj, obj.data.I, obj.data.x);
end
Iold = obj.data.I;

if isweighted(obj)
    obj.data.A = obj.data.A.*sqrt(obj.data.w);
    obj.data.y = obj.data.y.*sqrt(obj.data.w);
end
[obj.data.coeff,obj.data.qA,obj.data.rA] = SparseTools.ls_solve(obj.data.A,obj.data.y);
[obj.data.err,obj.data.opt_err] = SparseTools.ls_cross_validate(obj.data);
obj.l2_err = rel_error(obj);

if isdebug(obj.var)
    fprintf('+-----------+------------+------------+------------+------------+\n');
    fprintf('| Dim Basis | Nb Samples |  Dot prod. |  CV error  |  L2 error  |\n');
    fprintf('+-----------+------------+------------+------------+------------+\n');
else
    fprintf('+-----------+------------+------------+------------+\n');
    fprintf('| Dim Basis | Nb Samples |  Dot prod. |  CV error  |\n');
    fprintf('+-----------+------------+------------+------------+\n');
end
                
%
while (norm(obj.data.err) > obj.opt.tol) && (size(obj.data.x,2)<obj.opt.max_sample_size) ...
        && (cardinal(obj.data.I)<obj.opt.max_dim_basis) && basis_adaptation
    % Adaptive sampling on fixed basis
    % first add samples to the enriched basis functions
    if isempty(Iold)
        Idiff = obj.data.I;
    else
        Idiff = removeIndices(obj.data.I,Iold);
    end
    if cardinal(Idiff) > 0
        xadd = random(obj, Idiff, sample_factor);
        obj.data.x = [obj.data.x, xadd];
        yadd = func(xadd);
        %
        Aadd = eval_basis (obj, obj.data.I, xadd);
        wadd = eval_weight(obj, obj.data.I, xadd);
        if isweighted(obj)
            obj.data.y = [obj.data.y; yadd'.*sqrt(wadd)];
            obj.data.w = [obj.data.w(:); wadd(:)];
            Aadd = Aadd.*sqrt(wadd);
        else
            obj.data.y = [obj.data.y; yadd'];
        end
        obj.data.A = [obj.data.A; Aadd];
        if obj.opt.fast
            [obj.data.coeff,obj.data.qA,obj.data.rA] = SparseTools.ls_add_rows(obj.data.qA, obj.data.rA, Aadd, obj.data.y);
        else
            [obj.data.coeff,obj.data.qA,obj.data.rA] = SparseTools.ls_solve(obj.data.A,obj.data.y);
        end
        [obj.data.err,obj.data.opt_err] = SparseTools.ls_cross_validate(obj.data);
        obj.l2_err = rel_error(obj);
    end
    
    % adaptive sampling
    % err_stagn = Inf;
    flag = (obj.data.opt_err > obj.opt.opt_tol) || force_enrich_sample;
    while flag && norm(obj.data.err) > obj.opt.tol && size(obj.data.x,2) < obj.opt.max_sample_size
        %disp('adapt sampling')
        %err_old = err;
        %opt_err_old = opt_err;
        force_enrich_sample = false;
        %
        sample_factor = sample_factor + factor_add;
        %
        xadd = random(obj, obj.data.I, factor_add);
        obj.data.x = [obj.data.x, xadd];
        yadd = func(xadd);
        %
        Aadd = eval_basis (obj, obj.data.I, xadd);
        wadd = eval_weight(obj, obj.data.I, xadd);
        if isweighted(obj)
            obj.data.y = [obj.data.y; yadd'.*sqrt(wadd)];
            obj.data.w = [obj.data.w(:); wadd(:)];
            Aadd = Aadd.*sqrt(wadd);
        else
            obj.data.y = [obj.data.y; yadd'];
        end
        obj.data.A = [obj.data.A; Aadd];
        if obj.opt.fast
            [obj.data.coeff,obj.data.qA,obj.data.rA] = SparseTools.ls_add_rows(obj.data.qA, obj.data.rA, Aadd, obj.data.y);
        else
            [obj.data.coeff,obj.data.qA,obj.data.rA] = SparseTools.ls_solve(obj.data.A,obj.data.y);
        end
        [obj.data.err,obj.data.opt_err] = SparseTools.ls_cross_validate(obj.data);
        %{
        err_stagn = norm(err-err_old)/obj.data.err;
        opt_err_stagn = norm(opt_err-opt_err_old)/obj.data.opt_err;
        if obj.opt.display_iterations
            fprintf('|           | %10d | %4.4e |\n',size(obj.data.x,2),obj.data.err);
        end
        %}
        flag = (obj.data.opt_err > obj.opt.opt_tol) || force_enrich_sample;
    end
    Iold = obj.data.I;
    %
    obj.l2_err = rel_error(obj);
    if obj.opt.display_iterations
        if isdebug(obj.var)
            fprintf('|           | %10d | %4.4e | %4.4e | %4.4e |\n',size(obj.data.x,2),obj.data.opt_err,obj.data.err,obj.l2_err);
        else
            fprintf('|           | %10d | %4.4e | %4.4e |\n',size(obj.data.x,2),obj.data.opt_err,obj.data.err);
        end
    end
    % Adaptative basis with fixed sample, should we change weights? no!
    m = size(obj.data.x,2);
    while (norm(obj.data.err) > obj.opt.tol) && (obj.l2_err > obj.opt.tol)
        err_old = obj.data.err;
        Itest = obj.data.I;
        for kk = 1:obj.opt.enrich_degree
            switch lower(obj.opt.adaptation_rule)
                case 'margin'
                    Iadd = getMargin(Itest);
                case 'reducedmargin'
                    Iadd = getReducedMargin(Itest);
            end
            %
            rem = find((~(Iadd<=max_index)));
            Iadd = removeIndices(Iadd,rem);
            %
            if cardinal(Iadd)==0
                basis_adaptation = false;
                break;
            end
            %
            Inew = Itest.addIndices(Iadd);
            if cardinal(Inew) > m
                force_enrich_sample = true;
                warning('dimension of candidate basis + dimension of existing basis > number of samples');
                break;
            else
                Itest = Inew;
            end
        end
        if force_enrich_sample
            break;
        end
        Iadd = removeIndices(Itest,obj.data.I);
        Aadd = eval_basis(obj, Iadd, obj.data.x);
        if isweighted(obj)
            Aadd = Aadd.*sqrt(obj.data.w);
        end
        tmp_coeff = SparseTools.ls_add_cols(obj.data.qA, obj.data.rA, Aadd, obj.data.y);
        %
        [~,loc] = ismember(Iadd.array,Itest.array,'rows');
        c_marg = tmp_coeff(loc,:);
        norm_a_marg = sqrt(sum(c_marg.^2,2));
        switch lower(obj.opt.adaptation_rule)
            case 'margin'
                env = envelope(Iadd,norm_a_marg);
                [~,ind] = sort(env,'descend');
            case 'reducedmargin'
                [~,ind] = sort(norm_a_marg,'descend');
        end
        %size(norm_a_marg)
        %size(ind)
        if ~isempty(ind)
            energy = cumsum(norm_a_marg(ind).^2);
            %
            rep = find( energy >= energy(end)*obj.opt.bulk_parameter, 1, 'first' );
            Iadd.array = Iadd.array(ind(1:rep),:);
        end
        obj.data.I = addIndices(obj.data.I,Iadd);
        Aadd = eval_basis(obj, Iadd, obj.data.x);
        if isweighted(obj)
            Aadd = Aadd.*sqrt(obj.data.w);
        end
        obj.data.A = [obj.data.A, Aadd];
        [obj.data.coeff,obj.data.qA,obj.data.rA] = SparseTools.ls_add_cols(obj.data.qA, obj.data.rA, Aadd, obj.data.y);
        [obj.data.err,obj.data.opt_err] = SparseTools.ls_cross_validate(obj.data);
        %
        if obj.data.opt_err > obj.opt.opt_tol
            break
        end
        %
        err_stagn = norm(obj.data.err-err_old)/norm(obj.data.err);
        if (norm(obj.data.err) > obj.opt.overfit_tol*norm(err_old)) || (err_stagn <= obj.opt.stagnation_tol)
            break
        end
    end
    obj.l2_err = rel_error(obj);
    if obj.opt.display_iterations
        if isdebug(obj.var)
            fprintf('| %9d |            | %4.4e | %4.4e | %4.4e |\n',cardinal(obj.data.I),obj.data.opt_err,obj.data.err,obj.l2_err);
        else
            fprintf('| %9d |            | %4.4e | %4.4e |\n',cardinal(obj.data.I),obj.data.opt_err,obj.data.err);
        end
    end
    if obj.l2_err < obj.opt.tol
        break
    end
end
if obj.opt.display_iterations
    if isdebug(obj.var)
        fprintf('+-----------+------------+------------+------------+------------+\n');
    else
        fprintf('+-----------+------------+------------+------------+\n');
    end
end

if isdebug(obj.var)
    fprintf('| %9d | %10d | %4.4e | %4.4e | %4.4e |\n',cardinal(obj.data.I),size(obj.data.x,2),obj.data.opt_err,obj.data.err,obj.l2_err);
    fprintf('+-----------+------------+------------+------------+------------+\n');
else
    fprintf('| %9d | %10d | %4.4e | %4.4e |\n',cardinal(obj.data.I),size(obj.data.x,2),obj.data.opt_err,obj.data.err);
    fprintf('+-----------+------------+------------+------------+\n');
end

obj.n_eval = size(obj.data.x,2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

