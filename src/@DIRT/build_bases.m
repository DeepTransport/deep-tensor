function [bases, d] = build_bases(arg, ref)

ref_dom = get_domain(ref);

if isa(arg, 'cell')
    % This should contain 2 cells: for levels 0 and 1
    if (numel(arg)>2)
        warning('bases cells 3:%d are not used and will be ignored', numel(arg));
        bases = arg(1:2);
    elseif (numel(arg)<2)
        warning('repeat the first base');
        bases = repmat(arg(1), 1, 2);
    else
        bases = arg;
    end
    for i = 1:2
        if (isa(bases{i}, 'ApproxBases'))
            oneds = bases{i}.oneds;
            if i == 1
                doms = bases{i}.oned_domains;
            else
                doms = ref_dom;
            end
            d = ndims(bases{i});
        else
            error('bases cells element should be ApproxBases');
        end
        bases{i} = ApproxBases(oneds,doms,d);
    end
elseif (isa(arg, 'ApproxBases')) 
    bases = cell(1,2);
    oneds = arg.oneds;
    d = ndims(arg);
    for i = 1:2
        if i == 1
            doms = arg.oned_domains;
        else
            doms = ref_dom;
        end
        bases{i} = ApproxBases(oneds,doms,d);
    end
elseif isa(arg, 'numeric') % only gives the dimension
    bases = cell(1,2);
    d = arg;
    if isscalar(arg) && (arg > 0)
        for i = 1:2
            if i == 1
                doms = DIRT.defaultDomain;
            else
                doms = ref_dom;
            end
            bases{i} = ApproxBases(DIRT.defaultPoly,doms,arg);
        end
    else
        error('dimension should be a positive scalar')
    end
else
    error('wrong type of argument')
end

end