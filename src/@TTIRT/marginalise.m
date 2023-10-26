function obj = marginalise(obj, dir)
% Marginalise the pdf represented by FTT dimension by dimension.
%   irt = MARGINALISE(irt, int_dir)
%
%   int_dir - The direction of the marginalisation
%             >0: marginalise to the first from the last dimension 
%             <0: marginalise to the last from the first dimension 

if nargin > 1
    obj.int_dir = dir;
else
    obj.int_dir = 1;
end

d  = ndims(obj.approx);
obj.ys = cell(d, 1);
obj.ms = cell(d, 1);

if obj.int_dir > 0
    obj.order = 1:d;
else
    obj.order = d:-1:1;
end

if obj.int_dir > 0
    % ys{d} is built but shouldn't be used
    Ligeqk = 1;
    for k = d:-1:1
        nx  = cardinal(obj.approx.base.oneds{k});
        rkm = size(obj.approx.data.cores{k}, 1);
        rk  = size(obj.approx.data.cores{k}, 3);
        obj.ys{k} = reshape( reshape(obj.approx.data.cores{k}, rkm*nx, rk)*Ligeqk, rkm, nx);
        
        %disp(size(obj.approx.data.cores{k}))
        %disp(size(ys{k}))
        
        Ligeqk = integral(obj.approx.base.oneds{k}, obj.ys{k}')';
        obj.ms{k} = Ligeqk;
    end
    obj.fun_z = Ligeqk;
else
    % ys{1} is built but shouldn't be used
    Rileqk = 1;
    for k = 1:d
        nx  = cardinal(obj.approx.base.oneds{k});
        rkm = size(obj.approx.data.cores{k}, 1);
        rk  = size(obj.approx.data.cores{k}, 3);
        obj.ys{k} = reshape( Rileqk*reshape(obj.approx.data.cores{k}, rkm,nx*rk), nx, rk);
        
        %disp(size(obj.approx.data.cores{k}))
        %disp(size(ys{k}))
        
        Rileqk = integral(obj.approx.base.oneds{k}, obj.ys{k});
        obj.ms{k} = Rileqk;
    end
    obj.fun_z = Rileqk;
end

obj.z = obj.fun_z + obj.tau;

end