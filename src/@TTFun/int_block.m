function ftt = int_block(obj, ind)
% Integrate a block of FTT cores. tt = INT_BLOCK(tt, ind)
%
%   ind - indices of coordinates will be integrated
%   tt  - output FTT after integration

ftt = obj;
d = ndims(ftt.base);
%
if length(ind) ~= length(unique(ind))
    disp('integration indices are not unique')
    ind = unique(ind);
end
%
if max(ind) > d || min(ind) < 1
    error('integration indices out of bounds')
end
%
if length(ind) == d
    error('integration over all indices, should use int for integration')
end
ind = reshape(ind, 1, []);
if ind(end) == d
    jr = d;
    for i = (length(ind)-1):-1:1
        j = ind(i);
        if j+1 == jr
            % if the index is continous, move to the left
            jr = j;
        else
            % otherwise stop searching, jr is the last stopping index
            % jr:d are continous block
            break;
        end
    end
    ind2 = jr:d;
    ind1 = setdiff(ind, ind2);
    ind2 = fliplr(ind2);
else
    ind1 = ind;
    ind2 = [];
end
%the left blocks
if ~isempty(ind1)
    for i = 1:length(ind1)
        k   = ind1(i);
        nx  = cardinal(ftt.base.oneds{k});
        
        rkm = size(ftt.data.cores{k}, 1);
        rk  = size(ftt.data.cores{k}, 3);
        tmp = integral(ftt.base.oneds{k}, reshape(permute(ftt.data.cores{k}, [2,1,3]), nx,rkm*rk));
        tmp = reshape(tmp, rkm, rk);
        ftt.data.cores{k+1} = reshape(tmp*reshape(ftt.data.cores{k+1}, rk, []), rkm, nx, []);
    end
end
% the right blocks
if ~isempty(ind2)
    for i = 1:length(ind2)
        k = ind2(i);
        nx  = cardinal(ftt.base.oneds{k});
        rkm = size(ftt.data.cores{k}, 1);
        rk  = size(ftt.data.cores{k}, 3); % rk should always be 1
        tmp = integral(ftt.base.oneds{k}, reshape(permute(ftt.data.cores{k}, [2,1,3]), nx,rkm*rk));
        tmp = reshape(tmp, rkm, rk);
        ftt.data.cores{k-1} = reshape(reshape(ftt.data.cores{k-1}, [], rkm)*tmp, [], nx, rk);
    end
end
% delete marginalised blocks
ftt.data.cores(ind) = [];
ftt.base = remove_bases(ftt.base, ind);
ftt.data = clean(ftt.data);

end