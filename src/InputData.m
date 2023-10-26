classdef InputData
    % InputData class
    %
    % This class contains user-specified sampling data for initializing TT 
    % and debugging various approximations .
    %
    % InputData Properties:
    %   sample_x  - 
    %           Samples for initializing and enriching TT. Must be drawn
    %           from the domain of approximation.
    %   sample_z  - 
    %           Transformation of sample_x to the reference domain.
    %   count - A counter that indicate how much data have been used.
    %   debug_x - 
    %           Samples for estimating approximation errors. Must be drawn
    %           from the domain of approximation.
    %   debug_z - 
    %           Transformation of debug_x to the reference domain.
    %   debug_f - 
    %           Function evaluations at debug_z.
    %
    % InputData Methods:
    %   isdebug   - 
    %           If the debug samples are given.
    %   rel_error - 
    %           Computes the relative error.
    %   set_samples - 
    %           Sets the samples for TT, given the approximation base.
    %   get_samples - 
    %           Gets n samples for TT.
    %   set_debug   - 
    %           Sets debug samples and evaluate the target function at
    %           debug_z samples if needed. 
    %   isevaluated - 
    %           If the debug samples are evaluted.
    %   reset_counter - 
    %           Resets the counter to zero.
    
    properties
        sample_x
        sample_z
        debug_x
        debug_z
        debug_f
    end
        
    properties (Constant)
        count = Counter
    end

    methods
        function obj = InputData(sample_x, debug_x, debug_f)
            if nargin == 1
                obj.sample_x = sample_x;
                obj.debug_x = [];
                obj.debug_f = [];
            elseif nargin == 2
                obj.sample_x = sample_x;
                obj.debug_x = debug_x;
                obj.debug_f = [];
            elseif nargin == 3
                obj.sample_x = sample_x;
                obj.debug_x = debug_x;
                obj.debug_f = debug_f;
            else
                obj.sample_x = [];
                obj.debug_x = [];
                obj.debug_f = [];
            end
            obj.sample_z = [];
            obj.debug_z = [];
            obj.count.i = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %{
        function obj = set_debug_z(obj, z)
            obj.debug_z = z;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = set_debug_f(obj, f)
            obj.debug_f = f;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = set_sample_z(obj, z)
            obj.sample_z = z;
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = set_samples(obj, base, n)
            if isempty(obj.sample_x)
                disp('Generate initialization samples from the base measure');
                if nargin == 3     
                    obj.sample_z = sample_measure_reference(base, n);
                else
                    warning('There is no initialization samples available.\\Either provide sample set or specify a sample size.');
                end
            else
                if nargin == 3 && size(obj.sample_x,2) < n
                    disp('Not enough number of samples to initialise ftt')
                    obj.sample_z = sample_measure_reference(base, n);
                else
                    obj.sample_z = domain2reference(base, obj.sample_x);
                end
            end
            obj.count.i = 0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = get_samples(obj, n)
            nn = size(obj.sample_z,2);
            if (obj.count.i+n) > nn
                n1 = nn - obj.count.i + 1;
                n2 = n - n1;
                ind = [obj.count.i:nn, 1:n2];
                obj.count.i = n2;
                warning('Running out of samples. Start from beginning.')
            else
                ind = obj.count.i+(1:n);
                obj.count.i = obj.count.i+n;
            end
            z = obj.sample_z(:,ind);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function reset_counter(obj)
            obj.count.i = 0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = set_debug(obj, func, base)
            if ~isempty(obj.debug_x)
                obj.debug_z = domain2reference(base, obj.debug_x);
                if isempty(obj.debug_f)
                    obj.debug_f = feval(func, obj.debug_z);
                end
            else
                disp('Debug samples are not provided. Turn off evaluation-based cross validation.')
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function flag = isdebug(obj)
            flag = ~isempty(obj.debug_z);
            if isempty(obj.debug_z) && ~isempty(obj.debug_x)
                warning('Debug samples needs to be transformed into the reference domain.')
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function flag = isevaluated(obj)
            flag = ~isempty(obj.debug_f);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [l2, linf] = rel_error(obj, approx)
            if isdebug(obj)
                l2 = mean((obj.debug_f(:) - approx(:)).^2).^0.5/mean(obj.debug_f(:).^2).^0.5;
                linf = max(abs(obj.debug_f(:) - approx(:)))  / max(abs(obj.debug_f(:)));
            else
                l2 = inf;
                linf = inf;
            end
        end
    end
end