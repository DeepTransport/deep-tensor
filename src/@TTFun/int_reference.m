function z = int_reference(obj)
% Integrate the entire TT. z = INT(tt)

if obj.data.direction > 0
    % last dimenion
    figeqk  = 1;
    for k = ndims(obj.base):-1:1
        nx  = cardinal(obj.base.oneds{k});
        rkm = size(obj.data.cores{k}, 1);
        rk  = size(obj.data.cores{k}, 3);
        % push the previous dimenion into the current ftt by modifying the
        % coefficients, then ys{k} is used to replace the ftt core at the
        % k-th coordinate for obtaining the integrated function over the
        % last d-k dimensions, ys{k} is a rk-1 by nx matrix
        ys = reshape(reshape(obj.data.cores{k}, rkm*nx, rk)*figeqk, rkm, nx);
        figeqk = integral(obj.base.oneds{k}, ys')';
    end
    z = figeqk;
else
    % first dimenion
    fileqk = 1;
    for k = 1:ndims(obj.base)
        nx  = cardinal(obj.base.oneds{k});
        rkm = size(obj.data.cores{k}, 1);
        rk  = size(obj.data.cores{k}, 3);
        % push the previous dimenion into the current ftt
        % ys{k} is a nx by rk matrix
        ys = reshape(fileqk*reshape(obj.data.cores{k}, rkm, nx*rk), nx, rk);
        fileqk = integral(obj.base.oneds{k}, ys);
    end
    z = fileqk;
end
end