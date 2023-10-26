function [D]=lagr_diff(x, y)
%Differentiation matrix in the basis of Lagrange interpolating polynomials
%   function [D]=lagr_diff(x, y) Generates the differentiation matrix in the 
%   basis of the Lagrange interpolating polynomials on grid X, so that 
%   D(i,j) = L'_j(y_i)
%
%
% TT Toolbox 2.1, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------


n = numel(x);
M = numel(y);
D = zeros(M,n);

denom = ones(n,1);
for i=1:n % grid points in result
    for k=1:n % inner product
        if (k~=i)
            denom(i)=denom(i)*(x(i)-x(k));
        end
    end
end

for j=1:M % grid points in result
    for i=1:n % polynomial index
        curD=0;
        for p=1:n % sum
            if (p~=i)
                curnom = 1;
                for k=1:n % product
                    if (k~=i)&&(k~=p)
                        curnom=curnom*(y(j)-x(k));
                    end
                end
                curD=curD+curnom;
            end
        end
        D(j, i)=curD/denom(i);
    end
end

end