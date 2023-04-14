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
    % start with the last dim, upper triangular Chol of the mass matrix
    % ys{d} is built but shouldn't be used
    Ligeqk  = 1;
    for k = d:-1:1
        nx  = cardinal(obj.approx.oneds{k});
        rkm = size(obj.approx.cores{k}, 1);
        rk  = size(obj.approx.cores{k}, 3);
        % push the previous dimenion into the current ftt by modifying the
        % coefficients, then ys{k} is used to replace the ftt core at the
        % k-th coordinate for obtaining the integrated function over the
        % last d-k dimensions, ys{k} is a rk-1 by nx by rk tensor, as we
        % still need rk functions to evaluate the function
        obj.ys{k} = reshape( reshape(obj.approx.cores{k}, rkm*nx, rk)*Ligeqk, rkm, nx, rk);
        
        %disp(size(obj.approx.cores{k}))
        %disp(size(ys{k}))
        
        B = reshape( mass_r(obj.approx.oneds{k}, reshape(permute(obj.ys{k},[2,3,1]),nx,[])), [],rkm);
        [~, R] = qr(B,0);
        if size(R, 1) < size(R, 2)
            R = R(:,1:size(R, 1));
        end
        Ligeqk = R';
        obj.ms{k} = Ligeqk;
    end
    obj.fun_z = sum(sum(Ligeqk.^2, 1));
    %Ligeqk
else
    % start with the 1st dim, upper triangular Chol of the mass matrix
    % ys{1} is built but shouldn't be used
    Rileqk  = 1;
    for k = 1:d
        nx  = cardinal(obj.approx.oneds{k});
        rkm = size(obj.approx.cores{k}, 1);
        rk  = size(obj.approx.cores{k}, 3);
        % push the previous dimenion into the current ftt by modifying the
        % coefficients, then ys{k} is used to replace the ftt core at the
        % k-th coordinate for obtaining the integrated function over the
        % last d-k dimensions, ys{k} is a rk-1 by nx by rk tensor, as we
        % still need rk-1 functions to evaluate the function
        obj.ys{k} = reshape( Rileqk*reshape(obj.approx.cores{k}, rkm,nx*rk), [], nx, rk);
        
        %disp(size(obj.approx.cores{k}))
        %disp(size(ys{k}))
        
        B = reshape(mass_r(obj.approx.oneds{k},reshape(permute(obj.ys{k},[2,1,3]),nx,[])), [],rk);
        [~,Rileqk] = qr(B,0);
        %if size(Rileqk, 1) < size(Rileqk, 2)
        %   Rileqk = Rileqk(:,1:size(Rileqk, 1));
        %end
        
        %size(B)
        %size(Rileqk)
        obj.ms{k} = Rileqk;
    end
    obj.fun_z = sum(sum(Rileqk.^2, 1));
    %Rileqk
end

obj.z = obj.fun_z + obj.tau;

end