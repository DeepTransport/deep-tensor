classdef Wrench
    
    % Linear elasticity problem on a wrench
    
    properties
        
        
        PDEmodel  % model for the PDE to solve
        Geometry  % geometry of the wrench
        
        MeshSize  % diameter of the largest element
        CenterOfElements
        
        dimParam
        Sigma
        Sigma12 % Square root Sigma such that Sigma = Sigma12 * Sigma12'
        
        
        K % Stiffness matrix associated with the Young Modulus E=1
        M % Mass matrix
        B % 
        
    end
    
    properties (Access = private)
        % Those properties are meant to speedup the code (precomputations)
        TR
        deltaA  % non-assembled stiffness matrix: nbElem x 1 cells containing matrices you have to sum up to get A
        FEM
        
    end
    
    methods
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function wrench = Wrench(hmax,flag_deltaA)
            
            % Definition of the model
            if nargin<1
                hmax = 0.032;
            end
            if nargin<2
                flag_deltaA = 2;
            end
            
            % Initialization
            wrench = initGeo(wrench);
            wrench = initPDEmodel(wrench);
            wrench = initMesh(wrench,hmax);
            wrench = initKM(wrench)
            
            % parameter of the covariance function
            %  (x,y) ->  sigma^alpha * exp[ -( ||x-y||_2/l0 )^2 ]
            alpha = 2;
            sigma = 1;
            l0 = 1.0;
            wrench = initYoungsModulusField(wrench,alpha,sigma,l0);
            
            if flag_deltaA == 1 
                wrench = wrench.initDeltaA();
            elseif flag_deltaA == 2
                wrench = wrench.initDeltaAparallel();
            end
            
        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function results = evalSol(wrench,logE)
            
            % logE : Element-wise log Young's modulus (vector of size
            %        wrench.dimParam)
            if nargin<2
                logE = wrench.logE();
            end
            E = log(1 + exp(logE));
%             E = exp(logE);
            
            % Hooke's tensor: the parameter "c" (function handle)\
            
%             YoungsModulus = @(x,y) ones(length(x),1) ; % Young modulus
            
            
            
            % Hooke tensor for plain-stress problems, see https://www.mathworks.com/help/pde/ug/3-d-linear-elasticity-equations-in-toolbox-form.html?searchHighlight=linear%2520elasticity&s_tid=doc_srchtitle#bumr56f-1
            eta = 0.3; % Poissons ratio
            hookeEta = [ 1/(1-eta^2) 0 0 eta/(1-eta^2) ;...
                0 1/(1+eta)/2 1/(1+eta)/2 0 ; ...
                0 1/(1+eta)/2 1/(1+eta)/2 0 ; ...
                eta/(1-eta^2) 0 0 1/(1-eta^2)] ;
            hookeEta = reshape(hookeEta,[2,2,2,2]);
            hookeEta = permute(hookeEta,[1 3 2 4]); % I know... it is crazy... but that works!
            hookeEta = reshape(hookeEta,[4,4]);
            
            hookeTensor = @(Ex) repmat(hookeEta(:),1,length(Ex)) * spdiags(Ex(:),0,length(Ex),length(Ex)) ;
            
            c = @(region,state) hookeTensor( E(wrench.TR.pointLocation(region.x(:),region.y(:))) );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Specify all the coefficients
            specifyCoefficients(wrench.PDEmodel,'m',0,'d',0,'c',c,'a',0,'f',[0;0]);
            
            
            results = solvepde(wrench.PDEmodel);
            
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function g = evalGradSol(wrench,logE,L)
            
            % Compute y = L*u(logE) and the gradient of y w.r.t. logE using
            % the adjoint
            
            
            %%% Solution y = L*u
            E = log(1 + exp(logE));          
            D = 1./(1 + exp(-logE));
%             E = exp(logE);
            nbElem = wrench.dimParam;
            
            A = wrench.deltaA{1}*E(1);
            b = wrench.FEM.Fc;
            for k=2:nbElem
                A = A + wrench.deltaA{k}*E(k);
            end
            
            u  = A\b;
            L = wrench.FEM.B' * L;
            % y = L'*u;
            
            
            %%% Gradient
            q = A\L;
            
            deltaAu = cellfun(@(A,D) D*A*u ,wrench.deltaA,num2cell(D),'Uniformoutput',0);
            deltaAu = [deltaAu{:}];
            
            g = - q'*deltaAu ;
            
            
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function [y,g,Xobs,L] = evalObs(wrench,logE,type)
            
            % logE : Element-wise log Young's modulus (vector of size
            %        wrench.dimParam)
            if nargin<2
                logE = wrench.logE();
            end
            if nargin<3
                type = 3;
            end
            
            
            results = wrench.evalSol(logE);
            
            
            if type == 1
                % Vertical displacement at the point
                Xobs = [-4.1 ; 0 ];
                
                [~,indObs] = min(sum((wrench.PDEmodel.Mesh.Nodes - Xobs).^2));
                nbNodes = size(wrench.PDEmodel.Mesh.Nodes,2);
                L = sparse(nbNodes,2);
                L(indObs,2)=1;
                L=L(:);
                
                Xobs = wrench.PDEmodel.Mesh.Nodes(:,indObs);
                
                y = L' * results.NodalSolution(:);
                if nargout>1
                    g = wrench.evalGradSol(logE,L);
                    g = g(:);
                end
                
            elseif type == 2
                % Von Mises stress at the point
                Xobs = [1.3 ; 0.3];
                
                indObs = wrench.TR.pointLocation(Xobs(1),Xobs(2));
                if nargout==1
                    y = wrench.VonMises(results,logE,indObs);
                    g=[];
                else
                    [y,g] = wrench.VonMises(results,logE,indObs);
                    g = g(:);
                end
                
            elseif type == 3
                % Displacement along the line where the force is applied
                
                indObs = wrench.PDEmodel.Mesh.findNodes('region','Edge',14);
                nbNodes = size(wrench.PDEmodel.Mesh.Nodes,2);
                L = sparse(2*nbNodes,length(indObs));
                for i=1:length(indObs)
                    ell = sparse(nbNodes,2);
                    ell(indObs(i),2)=1;% the 2 refers to the vertical axis
                    ell = ell(:);
                    L(:,i) = ell;
                end
                
                Xobs = wrench.PDEmodel.Mesh.Nodes(:,indObs);
                
                y = L' * results.NodalSolution(:);
                if nargout>1
                    g = wrench.evalGradSol(logE,L);
                end
                
            end
            
            % For node indexing of the homebrew stiffness matrix
            L = wrench.FEM.B' * L;
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function [U,A,F] = evalFEMSystem(wrench, coeffEall)
            % Assembly and Solution interface for POD/ALS-Cross
            
            Nsamples = size(coeffEall, 3);
            U = cell(1, Nsamples);
            A = cell(1, Nsamples);
            F = cell(1, Nsamples);
            
            for i=1:Nsamples
                
%                 % Cast random variables into logE field
%                 logE = wrench.Sigma12*reshape(logEparam(:,:,i), [], 1);
                
                %%% Solution y = L*u
%                 E = exp(logE);
                E = reshape(coeffEall(:,:,i), [], 1);
                nbElem = wrench.dimParam;
                
                A{i} = wrench.deltaA{1}*E(1);
                F{i} = wrench.FEM.Fc;
                for k=2:nbElem
                    A{i} = A{i} + wrench.deltaA{k}*E(k);
                end
                
                U{i}  = A{i}\F{i};
                fprintf('wrench solved problem %d out of %d\n', i, Nsamples);
            end
            
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function [VM,g] = VonMises(wrench,results,logE,indElem)
            
            
            if nargin<4
                nbElements = size(wrench.CenterOfElements,2);
                indElem = 1:nbElements;
            end
            
            eta = 0.3;
            hookeEta = [ 1/(1-eta^2) 0 0 eta/(1-eta^2) ;...
                0 1/(1+eta)/2 1/(1+eta)/2 0 ; ...
                eta/(1-eta^2) 0 0 1/(1-eta^2)] ;
            
            [p,~,t] = meshToPet(wrench.PDEmodel.Mesh);
            [u1x,u1y] = pdegrad(p,t,results.NodalSolution(:,1));
            [u2x,u2y] = pdegrad(p,t,results.NodalSolution(:,2));
            
            
            VM = zeros(length(indElem),1);
            i=1;
            for k=indElem
                
                sigma = exp(logE(k)) * hookeEta * [u1x(k) u1y(k) u2x(k) u2y(k)]';
                VM(i) = sqrt(sigma'*[ 1 1/2 0 ; 1/2 1 0 ; 0 0 3 ]*sigma);
                i=i+1;
            end
            
            % If the gradient of the VM stress at one point is needed 
            if nargout>1 && length(indElem)==1
                
                % Diff of [u1x(k) u1y(k) u2x(k) u2y(k)]
                % GRAD_gradu : matrix of 4 columns 
                
                indNode = t(1:3,indElem);
                nbNode = size(p,2);
                
                gradX = sparse(nbNode,1);
                gradY = sparse(nbNode,1);
                zz = sparse(nbNode,1);
                NodalField = sparse(nbNode,1);
                for i=1:3
                    NodalField(indNode(i)) = 1;
                    
                    [aa,bb] = pdegrad(p,t(:,indElem),NodalField);
                    gradX(indNode(i)) = aa;
                    gradY(indNode(i)) = bb;
                    
                    NodalField(indNode(i)) = 0;
                end
                
                grads = [gradX , gradY  , zz , zz;...
                    zz , zz ,   gradX  , gradY];
                
                GRAD_gradu = wrench.evalGradSol(logE,grads);
                
                
                % Diff of sigma = exp(logE(k)) * hookeEta * [u1x(k) u1y(k) u2x(k) u2y(k)]';
                % GRAD_sigma : matrix of 3 columns
                GRAD_sigma = exp(logE(indElem)) * hookeEta * GRAD_gradu;
                GRAD_sigma(:,indElem) = GRAD_sigma(:,indElem) + sigma(:);
                
                % Diff of VM(i) = sqrt(sigma'*[ 1 1/2 0 ; 1/2 1 0 ; 0 0 3 ]*sigma);
                g = 1/VM * GRAD_sigma' * [ 1 1/2 0 ; 1/2 1 0 ; 0 0 3 ]*sigma;
                
                
            end
            
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function logE = logE(wrench, xi_ext)
            % Draw the parameter (the log of the Young's Modulus)
            if (nargin<2) || (isempty(xi_ext))
                xi = randn(size(wrench.Sigma12,2),1);
                logE = wrench.Sigma12*xi;
            else
                xi = xi_ext;
                logE = xi; % ?
            end
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function plot(wrench,data)
            
            
            if length(data) == size(wrench.PDEmodel.Mesh.Nodes,2)
                % Nodal-based field
                pdemesh(wrench.PDEmodel,'XYData',data,'ColorMap','jet')
                axis equal
                
            elseif length(data) == size(wrench.PDEmodel.Mesh.Elements,2)
                % Element-based field : plot a piecewise constance field
                
                [p,~,t]= wrench.PDEmodel.Mesh.meshToPet();
                t(4,:)=[];
                
                x=p(1,:);
                y=p(2,:);
                P=[x(t(:));y(t(:))];
                T=reshape(1:size(P,2),[3 size(P,2)/3]);
                
                tmp=repmat(data(:)',3,1);
                h = trisurf(T',P(1,:),P(2,:),tmp(:));
                set(h,'edgecolor','none')
                
                grid off
                view(0,90)
                axis([-5 5 -3 2])
                axis equal
                
                colorbar
                caxis([min(data) max(data)])
            end
            
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function plotGeometry(wrench)
            
            pdegplot(wrench.Geometry,'EdgeLabels','on','FaceLabels','off')
            
            axis equal
        end % endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function plotMesh(wrench)
            
            
            pdemesh(wrench.PDEmodel,'ElementLabels','off','NodeLabels','off')
            
            
            axis equal
        end % endFunction
        %------------------------------------------------------------------
        function plotSol(wrench,u,VM)
            
            plotVM = 0;
            if nargin==1
                x = wrench.x();
                u = wrench.evalSol(x);
                u = u.NodalSolution;
                VM = wrench.VonMises(u, x);
                plotVM = 1;
            end
            if nargin==2 && isa(u,'double')
                if numel(u)==wrench.dimParam
                    x = u;
                    u = wrench.evalSol(x);
                    VM = wrench.VonMises(u, x);
                    u = u.NodalSolution;
                    plotVM = 1;
                else
                    u = reshape(wrench.B*u,[size(wrench.B,1)/2 2]);
                    plotVM = 0;
                end
                
            end
            
            
            if nargin==3
                plotVM = 1;
            end
            
            [p,e,t] = meshToPet(wrench.PDEmodel.Mesh);
            factor = 0.0005;
            pnew = p+u'*factor;
            
            if plotVM
                
                t(4,:)=[];
                
                x = pnew(1,:);
                y = pnew(2,:);
                P = [x(t(:));y(t(:))];
                T = reshape(1:size(P,2),[3 size(P,2)/3]);
                
                tmp = repmat(VM(:)',3,1);
                h = trisurf(T',P(1,:),P(2,:),tmp(:));
                set(h,'edgecolor','none')
                grid off
                
                view(0,90)
                colorbar
%                 caxis([min(VM) max(VM)])
                
                hold off
                
            else
                pdeplot(pnew,e,t)
                hold off
                
            end
%             axis([-5 5 -3 2])
            axis equal
            
%             set(gca,'YTickLabel',[],'XTickLabel',[],'ytick',[],'xtick',[]);
            
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
%     end % endMethods
%     methods (Access = private) % Set the initialization functions...
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function wrench = initGeo(wrench)
            
            l=3.5;
            e=0.4;
            r1=1.4;
            r2=1;
            c2= 0.85;
            r11=1.1;
            r22=0.6;
            
            r = (e^2+(l-1.4)^2-r1^2)/(2*r1-2*e);
            dd = l-r1/(r1+r)*(l-1.4);
            
            % Elementary geometries (Rectangles, Circles etc.)
            R1 = [3;4;-l*.8;-l*.8;l*.8;l*.8;-e;e;e;-e];
            C1 = [1;l;0;r1];
            C2 = [1;-l;0;r2];
            P1 = [2 6 0.45+l+r11*cos([0:5]/6*(2*pi) + pi/6 ) -0.26+r11*sin([0:5]/6*(2*pi) + pi/6 )]';
            P2 = [2 6 -l+r22*cos([0:5]/6*(2*pi) ) r22*sin([0:5]/6*(2*pi) )]';
            P0 = [3 4 1.4 1.4 dd dd  -1 1 1 -1]';
            P00= [2 4 -l/2 -l+0.5 -2.7 -l+0.5 0 c2 0 -c2]';
            E1 = [1 1.4 e+r r]';
            E2 = [1 1.4 -e-r r]';
            
            gd = [[R1;zeros(length(P1)-length(R1),1)] , ...
                [C1;zeros(length(P1)-length(C1),1)],...
                [C2;zeros(length(P1)-length(C2),1)],...
                P1,...
                P2,...
                [P0;zeros(length(P1)-length(P0),1)],...
                [P00;zeros(length(P1)-length(P00),1)],...
                [E1;zeros(length(P1)-length(E1),1)],...
                [E2;zeros(length(P1)-length(E2),1)]    ];
            
            % Assemble the geometry
            ns = char('R1','C1','C2','P1','P2','P0','P00','E1','E2');
            ns = ns';
            sf = '(R1+C1+C2+P0+P00)-P1-P2-E1-E2';
            [dl,bt] = decsg(gd,sf,ns);
            [dl,~] = csgdel(dl,bt);
            
            wrench.Geometry = dl;
            
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function wrench = initPDEmodel(wrench)
            
            % Initialize
            wrench.PDEmodel = createpde(2); % the "2" means we have a set of two equations
            geometryFromEdges(wrench.PDEmodel,wrench.Geometry);
            
            % Set the boundary conditions
            applyBoundaryCondition(wrench.PDEmodel,'dirichlet','Edge',[10,11],'u',[0;0]);
            applyBoundaryCondition(wrench.PDEmodel,'neumann','Edge',[14],'g',[0;-1]);
            
%             % Hooke's tensor: the parameter "c" (function handle)\
%             eta = 0.3; % Poissons ratio
%             E = @(x,y) ones(length(x),1) ; % Young modulus
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % Hooke tensor for plain-stress problems, see https://www.mathworks.com/help/pde/ug/3-d-linear-elasticity-equations-in-toolbox-form.html?searchHighlight=linear%2520elasticity&s_tid=doc_srchtitle#bumr56f-1
%             hookeEta = [ 1/(1-eta^2) 0 0 eta/(1-eta^2) ;...
%                 0 1/(1+eta)/2 1/(1+eta)/2 0 ; ...
%                 0 1/(1+eta)/2 1/(1+eta)/2 0 ; ...
%                 eta/(1-eta^2) 0 0 1/(1-eta^2)] ;
%             hookeEta = reshape(hookeEta,[2,2,2,2]);
%             hookeEta = permute(hookeEta,[1 3 2 4]); % I know... it is crazy... but that works!
%             hookeEta = reshape(hookeEta,[4,4]);
%             
%             hookeTensor = @(Ex) repmat(hookeEta(:),1,length(Ex)) * spdiags(Ex(:),0,length(Ex),length(Ex)) ;
%             
%             c = @(region,state) hookeTensor( E(region.x,region.y) );
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%             
%             % Specify all the coefficients
%             specifyCoefficients(wrench.PDEmodel,'m',0,'d',0,'c',c,'a',0,'f',[0;0]);
%             
            
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function wrench = initMesh(wrench,hmax)
            
            generateMesh(wrench.PDEmodel,'Hmax',hmax,'GeometricOrder','linear');
            
            wrench.TR = triangulation(wrench.PDEmodel.Mesh.Elements',...
                wrench.PDEmodel.Mesh.Nodes');
            
            wrench.MeshSize = hmax;
            wrench.CenterOfElements = incenter(wrench.TR)';
            
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function wrench = initYoungsModulusField(wrench,alpha,sigma,l0)
            
            % Distance matrix between the centers of the elements
            DXY = pdist(wrench.CenterOfElements','euclidean');
            DXY = squareform(DXY);
            
            % Creating the covariance matrix
            wrench.Sigma = sigma^2 *exp( -(DXY/l0).^alpha );
            wrench.dimParam = size(wrench.Sigma,1);
            
            % Choleski factorization of Sigma (for fast sampling)
            [U,D] = eig(wrench.Sigma);
            D=diag(D);
            ind = D<(10*eps*max(D));
            
            D(ind)=[];
            U(:,ind)=[];
            wrench.Sigma12 = U*diag(sqrt(D));
            % wrench.SigmaM12 = U*diag(1./sqrt(D));
            % wrench.SigmaInv= U*diag(1./D)*U';
            
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function wrench = initDeltaA(wrench)
            
            nbElem = wrench.dimParam;
            wrench.deltaA = cell(nbElem,1);
            
            E = zeros(nbElem,1);
            for k=1:nbElem
                
                E(k) = 1;
                
                eta = 0.3; % Poissons ratio
                hookeEta = [ 1/(1-eta^2) 0 0 eta/(1-eta^2) ;...
                    0 1/(1+eta)/2 1/(1+eta)/2 0 ; ...
                    0 1/(1+eta)/2 1/(1+eta)/2 0 ; ...
                    eta/(1-eta^2) 0 0 1/(1-eta^2)] ;
                hookeEta = reshape(hookeEta,[2,2,2,2]);
                hookeEta = permute(hookeEta,[1 3 2 4]); % I know... it is crazy... but that works!
                hookeEta = reshape(hookeEta,[4,4]);
                
                hookeTensor = @(Ex) repmat(hookeEta(:),1,length(Ex)) * spdiags(Ex(:),0,length(Ex),length(Ex)) ;
                c = @(region,state) hookeTensor( E(wrench.TR.pointLocation(region.x(:),region.y(:))) );
                
                specifyCoefficients(wrench.PDEmodel,'m',0,'d',0,'c',c,'a',0,'f',[0;0]);
                FEM = assembleFEMatrices(wrench.PDEmodel , 'nullspace');
                wrench.deltaA{k} = FEM.Kc;
                
                
                
                E(k)=0;
                
            end
            FEM.Kc = [];
            wrench.FEM = FEM;
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function wrench = initDeltaAparallel(wrench)
            
            nbElem = wrench.dimParam;
            wrench.deltaA = cell(nbElem,1);
            
            parfor k=1:nbElem
                E = zeros(nbElem,1);
                
                E(k) = 1;
                
                eta = 0.3; % Poissons ratio
                hookeEta = [ 1/(1-eta^2) 0 0 eta/(1-eta^2) ;...
                    0 1/(1+eta)/2 1/(1+eta)/2 0 ; ...
                    0 1/(1+eta)/2 1/(1+eta)/2 0 ; ...
                    eta/(1-eta^2) 0 0 1/(1-eta^2)] ;
                hookeEta = reshape(hookeEta,[2,2,2,2]);
                hookeEta = permute(hookeEta,[1 3 2 4]); % I know... it is crazy... but that works!
                hookeEta = reshape(hookeEta,[4,4]);
                
                hookeTensor = @(Ex) repmat(hookeEta(:),1,length(Ex)) * spdiags(Ex(:),0,length(Ex),length(Ex)) ;
                c = @(region,state) hookeTensor( E(wrench.TR.pointLocation(region.x(:),region.y(:))) );
                
                specifyCoefficients(wrench.PDEmodel,'m',0,'d',0,'c',c,'a',0,'f',[0;0]);
                FEM = assembleFEMatrices(wrench.PDEmodel , 'nullspace');
%                 tmpFEM{k} = FEM
                tempdeltaA{k} = FEM.Kc;
                
                
                E(k)=0;
                
            end
            
            E = zeros(nbElem,1);
            E(nbElem) = 1;

            eta = 0.3; % Poissons ratio
            hookeEta = [ 1/(1-eta^2) 0 0 eta/(1-eta^2) ;...
                0 1/(1+eta)/2 1/(1+eta)/2 0 ; ...
                0 1/(1+eta)/2 1/(1+eta)/2 0 ; ...
                eta/(1-eta^2) 0 0 1/(1-eta^2)] ;
            hookeEta = reshape(hookeEta,[2,2,2,2]);
            hookeEta = permute(hookeEta,[1 3 2 4]); % I know... it is crazy... but that works!
            hookeEta = reshape(hookeEta,[4,4]);

            hookeTensor = @(Ex) repmat(hookeEta(:),1,length(Ex)) * spdiags(Ex(:),0,length(Ex),length(Ex)) ;
            c = @(region,state) hookeTensor( E(wrench.TR.pointLocation(region.x(:),region.y(:))) );

            specifyCoefficients(wrench.PDEmodel,'m',0,'d',0,'c',c,'a',0,'f',[0;0]);
            FEM = assembleFEMatrices(wrench.PDEmodel , 'nullspace');


            wrench.deltaA = tempdeltaA';  
            FEM.Kc = [];
            wrench.FEM = FEM;
        end
        
        %------------------------------------------------------------------
        function wrench = initKM(wrench)
            
            % M
            specifyCoefficients(wrench.PDEmodel,'m',0,'d',0,'c',0,'a',1,'f',[0;0]);
            FEM = assembleFEMatrices(wrench.PDEmodel , 'nullspace');
            wrench.M = FEM.Kc;
            wrench.B = FEM.B;
            
            % K0
            specifyCoefficients(wrench.PDEmodel,'m',0,'d',0,'c',1,'a',0,'f',[0;0]);
            FEM = assembleFEMatrices(wrench.PDEmodel , 'nullspace');
            wrench.K = FEM.Kc;
%             E = ones(size(wrench.Sigma12,1),1);
%             eta = 0.0; % Poissons ratio
%             
%             % Hooke tensor for plain-stress problems, see https://www.mathworks.com/help/pde/ug/3-d-linear-elasticity-equations-in-toolbox-form.html?searchHighlight=linear%2520elasticity&s_tid=doc_srchtitle#bumr56f-1
%             hookeEta = [ 1/(1-eta^2) 0 0 eta/(1-eta^2) ;...
%                 0 1/(1+eta)/2 1/(1+eta)/2 0 ; ...
%                 0 1/(1+eta)/2 1/(1+eta)/2 0 ; ...
%                 eta/(1-eta^2) 0 0 1/(1-eta^2)] ;
%             hookeEta = reshape(hookeEta,[2,2,2,2]);
%             hookeEta = permute(hookeEta,[1 3 2 4]);
%             hookeEta = reshape(hookeEta,[4,4]);
%             hookeTensor = @(Ex) repmat(hookeEta(:),1,length(Ex)) * spdiags(Ex(:),0,length(Ex),length(Ex)) ;
%             c = @(region,state) hookeTensor( E(wrench.TR.pointLocation(region.x(:),region.y(:))) );
%             
%             specifyCoefficients(wrench.PDEmodel,'m',0,'d',0,'c',c,'a',0,'f',[0;0]);
%             FEM = assembleFEMatrices(wrench.PDEmodel , 'nullspace');
%             wrench.K0 = FEM.Kc;
            
            
        end
        %------------------------------------------------------------------

        %------------------------------------------------------------------
        %------------------------------------------------------------------
    end % endMethods
end % endClass






