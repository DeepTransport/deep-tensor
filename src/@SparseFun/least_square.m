function obj = least_square(obj,func,deb)

% Implement Alg 4 of "Adaptive approximation by optimal weighted least 
% squares methods". Giovanni Migliorati
% The code is modified based on the AdaptiveSparseTensorAlgorithm of the
% ApproximationToolbox

d = ndims(obj);
max_index = MultiIndices(reshape(cardinals(obj)-1,1,d));

if cardinal(obj.indices) == 0
    Iold = [];
    I = MultiIndices(zeros(1,d));
    for i = 1:obj.opt.init_total_degree
        Iadd = getReducedMargin(I);
        I = I.addIndices(Iadd);
    end
    rem = find((~(I<=max_index)));
    Iadd = removeIndices(I,rem);
    x = [];
    y = [];
    A = [];
    w = [];
    err = Inf;
    l2_err = Inf;
else
    I = obj.indices;
    l2_err = Inf;
    if isempty(obj.x)
        Iold = [];
        x = [];
        y = [];
        A = [];
        w = [];
        err = Inf;
    else
        Iold = I;
        x = obj.x;
        y = obj.y;
        A = eval_basis(obj, I, x);
        [coeff,err] = SparseFun.ls_solve(A,y',[]);
        if ~isempty(deb)            
            basis_at_z = eval_basis(obj, I, deb.samples);
            approx = reshape(basis_at_z*coeff,[],size(deb.samples,2));
            l2_err = sqrt(mean((deb.f(:) - approx(:)).^2))/sqrt(mean(deb.f(:)));
        end
    end
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
                
% sample_factor: average number of samples per basis function
sample_factor = obj.opt.init_sample_size;
factor_add = obj.opt.enrich_sample_size;
basis_adaptation = true;
%
while (norm(err) > obj.opt.tol) && (size(x,2)<obj.opt.max_sample_size) && (cardinal(I)<obj.opt.max_dim_basis) && basis_adaptation
    % Adaptive sampling on fixed basis
    % first add samples to the enriched basis functions
    if isempty(Iold)
        Idiff = I;
    else
        Idiff = removeIndices(I,Iold);
    end
    if cardinal(Idiff) > 0
        xadd = sample_measure_reference(obj, ceil(cardinal(Idiff)*sample_factor));
        x = [x, xadd];
        yadd = func(xadd);
        y = [y, yadd];
        %
        Aadd = eval_basis(obj, I, xadd);
        A = [A;Aadd];
        [coeff,err] = SparseFun.ls_solve(A,y',[]);
        if ~isempty(deb)            
            basis_at_z = eval_basis(obj, I, deb.samples);
            approx = reshape(basis_at_z*coeff,[],size(deb.samples,2));
            l2_err = sqrt(mean((deb.f(:) - approx(:)).^2))/sqrt(mean(deb.f(:)));
        end
    end
    opt_err = norm( A'*(A./size(x,2)) - eye(size(A,2)), 2);
    
    % adaptive sampling
    % err_stagn = Inf;
    while (opt_err > obj.opt.opt_tol) && norm(err) > obj.opt.tol && size(x,2) < obj.opt.max_sample_size
        %disp('adapt sampling')
        %err_old = err;
        %opt_err_old = opt_err;
        %
        sample_factor = sample_factor + factor_add;
        xadd = sample_measure_reference(obj, ceil(cardinal(I)*factor_add));
        x = [x, xadd];
        yadd = func(xadd);
        y = [y, yadd];
        %
        Aadd = eval_basis(obj, I, xadd);
        A = [A;Aadd];
        [coeff,err] = SparseFun.ls_solve(A,y',[]);
        %
        opt_err = norm( A'*(A./size(x,2)) - eye(size(A,2)), 2);
    end
    Iold = I;
    %
    if ~isempty(deb)
        basis_at_z = eval_basis(obj, I, deb.samples);
        approx = reshape(basis_at_z*coeff,[],size(deb.samples,2));
        l2_err = sqrt(mean((deb.f(:) - approx(:)).^2))/sqrt(mean(deb.f(:)));
    end
    if obj.opt.display_iterations
        if isempty(deb)
            fprintf('|           | %10d | %4.4e | %4.4e |\n',size(x,2),norm(opt_err),norm(err));
        else
            fprintf('|           | %10d | %4.4e | %4.4e | %4.4e |\n',size(x,2),norm(opt_err),norm(err),l2_err);
        end
    end
    % Adaptative basis with fixed sample, should we change weights? no!
    m = size(x,2);
    while (norm(err) > obj.opt.tol) 
        err_old = err;
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
                break;
            else
                Itest = Inew;
            end
        end
        Iadd = removeIndices(Itest,I);
        Aadd = eval_basis(obj, Iadd, x);
        Atest = [A,Aadd];
        %Atest2 = eval_basis(obj, Itest, x);
        %norm(Atest-Atest2)
        coeff = SparseFun.ls_solve(Atest,y',[]);
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
        A = [A,Aadd];
        %A2 = eval_basis(obj, I, x);
        %norm(A-A2)
        [coeff,err] = SparseFun.ls_solve(A,y',[]);
        opt_err = norm( A'*(A./size(x,2)) - eye(size(A,2)), 2);
        if opt_err > obj.opt.opt_tol
            break
        end
        %
        err_stagn = norm(err-err_old)/norm(err);
        if (norm(err) > obj.opt.overfit_tol*norm(err_old)) || (err_stagn <= obj.opt.stagnation_tol)
            break
        end
    end    
    if ~isempty(deb)
        basis_at_z = eval_basis(obj, I, deb.samples);
        approx = reshape(basis_at_z*coeff,[],size(deb.samples,2));
        l2_err = sqrt(mean((deb.f(:) - approx(:)).^2))/sqrt(mean(deb.f(:)));
    end
    if obj.opt.display_iterations
        if isempty(deb)
            fprintf('| %9d |            | %4.4e | %4.4e |\n',cardinal(I),norm(opt_err),norm(err));
        else
            fprintf('| %9d |            | %4.4e | %4.4e | %4.4e |\n',cardinal(I),norm(opt_err),norm(err),l2_err);
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
    fprintf('| %9d | %10d | %4.4e | %4.4e |\n',cardinal(I),size(x,2),norm(opt_err),norm(err));
    fprintf('+-----------+------------+------------+------------+\n');
else
    fprintf('| %9d | %10d | %4.4e | %4.4e | %4.4e |\n',cardinal(I),size(x,2),norm(opt_err),norm(err),l2_err);
    fprintf('+-----------+------------+------------+------------+------------+\n');
end

obj.indices = I;
obj.data = coeff;
obj.err = err;
obj.A = A;
obj.x = x;
obj.w = [];
obj.y = y;
obj.n_evals = size(x,2);

end
