function [quadx,quadw] = quadcc_piecewise(poly, nquad, x)

log_order = ceil(log2(nquad-1));
nquad = 2^log_order + 1;

quadx = ones(nquad*poly.num_elems, length(x))*poly.grid(1);
quadw = zeros(nquad*poly.num_elems, length(x));

% index in which the boundary x falls into
% ei = 0: left ghost cell, ei = poly.num_elems+1: right ghost cell
ei = ceil( (x-poly.grid(1))./poly.elem_size );

cc = 0;
[z,w] = cc_rule(log_order);
for k = 1:poly.num_elems
    ind1 = (ei>k);  % integrate the whole element
    ind2 = (ei==k); % integrate to x
    n1 = sum(ind1);
    n2 = sum(ind2);
    ri = (1:nquad) + cc*nquad;
    cc = cc+1;
    if n1 > 0
        a = poly.grid(k);
        b = poly.grid(k+1);
        dxdz = (b-a)/2;
        mid = (b+a)/2;
        quadx(ri,ind1) = repmat(z*dxdz + mid, 1, n1);
        quadw(ri,ind1) = repmat(w*dxdz, 1, n1);
    end
    if n2 > 0
        a = poly.grid(k)*ones(1,n2);
        b = x(ind2);
        dxdz = (b-a)/2;
        mid = (b+a)/2;
        quadx(ri,ind2) = z.*dxdz + mid;
        quadw(ri,ind2) = w.*dxdz;
    end
end

end

%{

function [quadx,quadw] = quadcc_piecewise(poly, nquad, x)

tol = 1E-6;
log_order = ceil(log2(nquad-1));
nquad = 2^log_order + 1;

if poly.gs > tol
    quadx = ones(nquad*(poly.num_elems+2), length(x))*poly.grid(1);
    quadw = zeros(nquad*(poly.num_elems+2), length(x));
else
    quadx = ones(nquad*poly.num_elems, length(x))*poly.grid(1);
    quadw = zeros(nquad*poly.num_elems, length(x));
end

% index in which the boundary x falls into
% ei = 0: left ghost cell, ei = poly.num_elems+1: right ghost cell
ei = ceil( (x-poly.grid(1))./poly.elem_size );

cc = 0;
[z,w] = cc_rule(log_order);
for k = 0:(poly.num_elems+1)
    ind1 = (ei>k);  % integrate the whole element
    ind2 = (ei==k); % integrate to x
    n1 = sum(ind1);
    n2 = sum(ind2);
    if k == 0
        % ghost cell
        if poly.gs > tol
            ri = (1:nquad) + cc*nquad;
            cc = cc+1;
            if n1 > 0
                a = poly.domain(1);
                b = poly.grid(1);
                dxdz = (b-a)/2;
                mid = (b+a)/2;
                quadx(ri,ind1) = repmat(z*dxdz + mid, 1, n1);
                quadw(ri,ind1) = repmat(w*dxdz, 1, n1);
            end
            if n2 > 0
                a = poly.domain(1)*ones(1,n2);
                b = x(ind2);
                dxdz = (b-a)/2;
                mid = (b+a)/2;
                quadx(ri,ind2) = z.*dxdz + mid;
                quadw(ri,ind2) = w.*dxdz;
            end
        end
    elseif k == (poly.num_elems+1)
        % ghost cell
        if poly.gs > tol
            ri = (1:nquad) + cc*nquad;
            cc = cc+1;
            if n1 > 0
                a = poly.grid(end);
                b = poly.domain(2);
                dxdz = (b-a)/2;
                mid = (b+a)/2;
                quadx(ri,ind1) = repmat(z*dxdz + mid, 1, n1);
                quadw(ri,ind1) = repmat(w*dxdz, 1, n1);
            end
            if n2 > 0
                a = poly.grid(end)*ones(1,n2);
                b = x(ind2);
                dxdz = (b-a)/2;
                mid = (b+a)/2;
                quadx(ri,ind2) = z.*dxdz + mid;
                quadw(ri,ind2) = w.*dxdz;
            end
        end
    else
        ri = (1:nquad) + cc*nquad;
        cc = cc+1;
        if n1 > 0
            a = poly.grid(k);
            b = poly.grid(k+1);
            dxdz = (b-a)/2;
            mid = (b+a)/2;
            quadx(ri,ind1) = repmat(z*dxdz + mid, 1, n1);
            quadw(ri,ind1) = repmat(w*dxdz, 1, n1);
        end
        if n2 > 0
            a = poly.grid(k)*ones(1,n2);
            b = x(ind2);
            dxdz = (b-a)/2;
            mid = (b+a)/2;
            quadx(ri,ind2) = z.*dxdz + mid;
            quadw(ri,ind2) = w.*dxdz;
        end
    end
end

end

%}
