function [ind,B] = maxvol(H)
%Build the cross indices using maxvol
%
%Tiangang Cui, August, 2020

tol     = 1E-6;
maxiter = 200;
[m,n]   = size(H);

% initialisation using QR-DEIM
[~,~,e] = qr(H', 'vector');
ind = e(1:n);
B   = H/H(ind,:);

for k = 1:maxiter
    [mb,mi] = max(abs(B(:)));
    if ( mb <= 1 + tol )
        return
    end
    isub = floor((mi-0.5)/m) + 1; % col index of Bmax
    inew = mi - (isub-1)*m; % row index of Bmax
    %
    iold = ind(isub);
    %
    % we use H(inew,:) to replace H(iold,:)
    % 
    % H(ind_new,:) = H(ind_old,:) + e_sub * (H(inew,:) - H(iold,:))
    % rank-1 update of the inversion
    % write the above as 
    % inv( A + e * b^t ) = inv(A) - (inv(A) e) * (b^t inv(A)) / ( 1 + b^t inv(A) e )
    % we have H * inv(A) = B and b^t = (H(inew,:) - H(iold,:)), so
    % b^t inv(A) = (B(inew,:) - B(iold,:))
    % b^t inv(A) e + 1= (B(inew,i_sub) - B(iold,i_sub)) + 1
    % B(iold,i_sub) = 1 as it comes from B = H inv(H(ind,:)) and iold in ind
    %
    % Then, we have
    % H inv( A + e * b^t ) = H inv(A) - H (inv(A) e) * (B(inew,:) - B(iold,:)) / B(inew,isub)
    %                      = B - B(:,isub) * (B(inew,:) - B(iold,:)) / B(inew,isub)
    
    %disp([k, isub, inew, mb])
    
    B = B - B(:,isub) * ( (B(inew,:) - B(iold,:)) / B(inew,isub) );
    ind(isub) = inew;
    
    % debug
    %B0 = H/H(ind,:);
    %norm(B - B0, 'fro')
end
end