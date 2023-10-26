classdef SparseTools
    % SparseTools class: tools used for building sparse approximation.
    %
    % SparseTools Methods (Static):
    %   ls_solve - 
    %           Solves the least square system ||Ax -y||.
    %   ls_add_cols  - 
    %           Re-solves the least square system ||Ax -y|| after adding
    %           columns to A, which also increases the number of coefficients.
    %   ls_add_rows  - 
    %           Re-solves the least square system ||Ax -y|| after adding
    %           rows to A and y.
    %   ls_cross_validate -
    %           Estimates L2 err and projector error. 
    %   total_degree -     
    %           Generates a total-degree index set.
    %   hyperbolic_cross  -
    %           Generates a hyperblic-cross index set.
    %
    % See also SPARSEFUN
    
    methods (Static)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [coeff,qA,rA] = ls_solve(A,y)
            [qA,rA] = qr(A,0);
            coeff = rA\(qA'*y);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [coeff,qA,rA] = ls_add_cols(qA, rA, Aadd, y)
            R1 = qA'*Aadd;
            [Q2, R2] = qr(Aadd - qA*R1, 0);
            %
            rA = [rA, R1; zeros(size(Aadd,2),size(qA,2)), R2];
            qA = [qA, Q2];
            %
            coeff = rA\(qA'*y);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [coeff,qA,rA] = ls_add_rows(qA, rA, Aadd, y)
            n = size(qA,2);
            [Q1,R1] = qr([rA;Aadd],0);
            %
            %[m,n] = size(qA);
            %r = size(Aadd,1);
            % qA = [qA, zeros(m,r); zeros(r,n), eye(r)]*Q1(:,1:n);
            qA = [qA*Q1(1:n,1:n); Q1((n+1):end,1:n)];
            rA = R1(1:n,1:n);
            %
            coeff = rA\(qA'*y);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [err, opt_err] = ls_cross_validate(data)
            [m,n] = size(data.A);
            nout = size(data.y,2);
            second_moment = sum(data.y.^2,1)/size(data.y,1);
            % Compute the relative leave-one-out cross-validation error the fast
            % leave-one-out cross-validation procedure [Cawlet & Talbot 2004] based on
            % the Sherman-Morrison-Woodbury formula
            %
            err = zeros(1,nout);
            if m-1 < n
                warning('Not enough samples for cross-validation');
                err(:) = Inf;
                return
            end
            % Compute the predicted residuals using the Bartlett matrix inversion formula
            % (special case of the Sherman-Morrison-Woodbury formula)
            T = sum(data.qA.^2,2);
            del = (data.y-data.qA*(data.rA*data.coeff))./(1-T);
            err = sum(del.^2,1)/m;
            err = err./second_moment;
            %
            % Compute the relative cross-validation error
            %if rcond(data.rA) > 1E-16
            dra = diag(data.rA);
            if min(abs(dra))/max(abs(dra)) > 1E-16
                %invrA = inv(rA);
                %C = invrA*invrA';
                % (1-n/m)^(-1)*(1+trace(C))
                % Direct Eigenvalue Estimator [Chapelle, Vapnik & Bengio, 2002], [Blatman & Sudret, 2011]
                % -> accurate even when m is not >> n (corr ~ 1+2*n/m when m->Inf, as trace(C) ~ n/m when m->Inf)
                %corr2 = (m/(m-n))*(1+trace(C));
                %norm(sum(diag(rA).^(-1))^2 - trace(C))
                %norm(sum(diag(rA).^(-1)) - trace(invrA))
                %norm(corr-corr2)
                corr = (m/(m-n))*(1+sum(dra.^(-1))^2);
                if corr > 0
                    err = err*corr;
                end
            end
            err = sqrt(err);
            %
            opt_err = normest( data.rA'*(data.rA/size(data.A,1)) - eye(size(data.A,2)), 1e-2);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function I = total_degree(d, k)
            I = MultiIndices(zeros(1,d));
            for i = 1:k
                Iadd = getReducedMargin(I);
                I = I.addIndices(Iadd);
            end 
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function I = hyperbolic_cross(d, k, weights)
            
            if nargin == 2
                weights = ones(1,d);
            end
            array = SparseTools.hyperbolic_cross_recur(d, k, weights);
            
            I = MultiIndices(unique(array,'rows'));
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function I = hyperbolic_cross_recur(d, k, weights)
            
            if d == 1
                I = reshape(0:k, [], 1);
            else
                I_pre = SparseTools.hyperbolic_cross_recur(d-1,k,weights);
                n = size(I_pre,1);
                I = zeros(n*(k+1),d);
                %
                pos = 0;
                for j = 1:n
                    tmp = prod(max(I_pre(j,:), 1)./weights(1:d-1));
                    for i = 1:k+1
                        if tmp*(max(i-1,1)/weights(d)) <= k
                            pos = pos + 1;
                            I(pos,1:d-1) = I_pre(j,:);
                            I(pos,d) = i-1;
                        end
                    end
                end
                I = I(1:pos,:);
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end