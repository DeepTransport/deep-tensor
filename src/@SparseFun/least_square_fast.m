function obj = least_square_fast(obj,func,deb)

% Implement Alg 4 of "Adaptive approximation by optimal weighted least 
% squares methods". Giovanni Migliorati
% The code is modified based on the AdaptiveSparseTensorAlgorithm of the
% ApproximationToolbox

d = ndims(obj);
max_index = MultiIndices(reshape(cardinals(obj)-1,1,d));
ls = LeastSquare();

% sample_factor: average number of samples per basis function
sample_factor = ceil(obj.opt.init_sample_size);
factor_add = ceil(obj.opt.enrich_sample_size);
basis_adaptation = true;
force_enrich_sample = false;

if cardinal(obj.indices) == 0
    Iold = [];
    switch obj.opt.indexset
        case {'hyperbolic'}
            I = hyperbolic_cross(d, obj.opt.init_total_degree);
        otherwise
            I = total_degree(d, obj.opt.init_total_degree);
    end
    rem = find((~(I<=max_index)));
    I = removeIndices(I,rem);
    l2_err = Inf;
    x = random(obj.importance, I, sample_factor);
    y = func(x);
    %
    A = eval_basis(obj, I, x);
    w = eval_weight(obj.importance, I, x);
    w = w(:);
else
    I = obj.indices;
    l2_err = Inf;
    if isempty(obj.x)
        x = random(obj.importance, I, sample_factor);
        y = func(x);
        %
        A = eval_basis(obj, I, x);
        w = eval_weight(obj.importance, I, x);
        w = w(:);
    else
        x = obj.x;
        y = obj.y;
        A = eval_basis(obj, I, x);
        w = eval_weight(obj.importance, I, x);
        w = w(:);
    end
end

Iold = I;
ls = initialize(ls, A, y', w);
if ~isempty(deb)
    basis_at_z = eval_basis(obj, I, deb.samples);
    approx = reshape(basis_at_z*ls.x,[],size(deb.samples,2));
    l2_err = sqrt(mean((deb.f(:) - approx(:)).^2))/sqrt(mean(deb.f(:)));
end

if isempty(deb)
    fprintf('+-----------+------------+------------+------------+\n');
    fprintf('| Dim Basis | Nb Samples |  Dot prod. |  CV error  |\n');
    fprintf('+-----------+------------+------------+------------+\n');
else
    fprintf('+-----------+------------+------------+------------+------------+\n');
    fprintf('| Dim Basis | Nb Samples |  Dot prod. |  CV error  |  L2 error  |\n');
    fprintf('+-----------+------------+------------+------------+------------+\n');
end
                
%
while (norm(ls.err) > obj.opt.tol) && (size(x,2)<obj.opt.max_sample_size) && (cardinal(I)<obj.opt.max_dim_basis) && basis_adaptation
    % Adaptive sampling on fixed basis
    % first add samples to the enriched basis functions
    if isempty(Iold)
        Idiff = I;
    else
        Idiff = removeIndices(I,Iold);
    end
    if cardinal(Idiff) > 0
        xadd = random(obj.importance, Idiff, sample_factor);
        x = [x, xadd];
        yadd = func(xadd);
        %
        Aadd = eval_basis(obj, I, xadd);
        wadd = eval_weight(obj.importance, I, xadd);
        ls = add_rows(ls, Aadd, yadd', wadd);
        if ~isempty(deb)            
            basis_at_z = eval_basis(obj, I, deb.samples);
            approx = reshape(basis_at_z*ls.x,[],size(deb.samples,2));
            l2_err = sqrt(mean((deb.f(:) - approx(:)).^2))/sqrt(mean(deb.f(:)));
        end
    end
    
    % adaptive sampling
    % err_stagn = Inf;
    flag = (ls.opt_err > obj.opt.opt_tol) || force_enrich_sample;
    while flag && norm(ls.err) > obj.opt.tol && size(x,2) < obj.opt.max_sample_size
        %disp('adapt sampling')
        %err_old = err;
        %opt_err_old = opt_err;
        force_enrich_sample = false;
        %
        sample_factor = sample_factor + factor_add;
        %
        xadd = random(obj.importance, I, factor_add);
        x = [x, xadd];
        yadd = func(xadd);
        %
        Aadd = eval_basis(obj, I, xadd);
        wadd = eval_weight(obj.importance, I, xadd);
        ls = add_rows(ls, Aadd, yadd', wadd);
        %{
        err_stagn = norm(err-err_old)/ls.err;
        opt_err_stagn = norm(opt_err-opt_err_old)/ls.opt_err;
        if obj.opt.display_iterations
            fprintf('|           | %10d | %4.4e |\n',size(x,2),ls.err);
        end
        %}
        flag = (ls.opt_err > obj.opt.opt_tol) || force_enrich_sample;
    end
    Iold = I;
    %
    if ~isempty(deb)
        basis_at_z = eval_basis(obj, I, deb.samples);
        approx = reshape(basis_at_z*ls.x,[],size(deb.samples,2));
        l2_err = sqrt(mean((deb.f(:) - approx(:)).^2))/sqrt(mean(deb.f(:)));
    end
    if obj.opt.display_iterations
        if isempty(deb)
            fprintf('|           | %10d | %4.4e | %4.4e |\n',size(x,2),ls.opt_err,ls.err);
        else
            fprintf('|           | %10d | %4.4e | %4.4e | %4.4e |\n',size(x,2),ls.opt_err,ls.err,l2_err);
        end
    end
    % Adaptative basis with fixed sample, should we change weights? no!
    m = size(x,2);
    while (norm(ls.err) > obj.opt.tol) 
        err_old = ls.err;
        Itest = I;
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
        Iadd = removeIndices(Itest,I);
        Aadd = eval_basis(obj, Iadd, x);
        %Atest = [A,Aadd];
        %Atest2 = eval_basis(obj, Itest, x);
        %norm(Atest-Atest2)
        %
        %wtest = eval_weight(obj.importance, Itest, x);
        %wtest = wtest(:);
        coeff = add_cols_solve(ls, Aadd);
        %
        [~,loc] = ismember(Iadd.array,Itest.array,'rows');
        c_marg = coeff(loc,:);
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
        I = I.addIndices(Iadd);
        Aadd = eval_basis(obj, Iadd, x);
        ls = add_cols(ls, Aadd);
        if ls.opt_err > obj.opt.opt_tol
            break
        end
        %
        err_stagn = norm(ls.err-err_old)/norm(ls.err);
        if (norm(ls.err) > obj.opt.overfit_tol*norm(err_old)) || (err_stagn <= obj.opt.stagnation_tol)
            break
        end
    end
    if ~isempty(deb)
        basis_at_z = eval_basis(obj, I, deb.samples);
        approx = reshape(basis_at_z*ls.x,[],size(deb.samples,2));
        l2_err = sqrt(mean((deb.f(:) - approx(:)).^2))/sqrt(mean(deb.f(:)));
    end
    if obj.opt.display_iterations
        if isempty(deb)
            fprintf('| %9d |            | %4.4e | %4.4e |\n',cardinal(I),ls.opt_err,ls.err);
        else
            fprintf('| %9d |            | %4.4e | %4.4e | %4.4e |\n',cardinal(I),ls.opt_err,ls.err,l2_err);
        end
    end
    if (~isempty(deb) && l2_err < obj.opt.tol)
        break
    end
end
if obj.opt.display_iterations
    if isempty(deb)
        fprintf('+-----------+------------+------------+------------+\n');
    else
        fprintf('+-----------+------------+------------+------------+------------+\n');
    end
end

if isempty(deb)
    fprintf('| %9d | %10d | %4.4e | %4.4e |\n',cardinal(I),size(x,2),ls.opt_err,ls.err);
    fprintf('+-----------+------------+------------+------------+\n');
else
    fprintf('| %9d | %10d | %4.4e | %4.4e | %4.4e |\n',cardinal(I),size(x,2),ls.opt_err,ls.err,l2_err);
    fprintf('+-----------+------------+------------+------------+------------+\n');
end

obj.indices = I;
obj.data = ls.x;
obj.err = ls.err;
obj.x = x;
obj.ls = ls;
%obj.A = A;
%obj.w = w;
%obj.y = y;
obj.n_evals = size(x,2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

