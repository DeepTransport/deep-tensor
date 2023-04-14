classdef Debugger
    properties
        samples
        f
    end
    
    methods
        function obj = Debugger(samples, f)
            if nargin == 1
                obj.samples = samples;
                obj.f = [];
            elseif nargin == 2
                obj.samples = samples;
                obj.f = f;
            else
                obj.samples = [];
                obj.f = [];
            end
        end
        
        function obj = set_samples(obj, samples)
            obj.samples = samples;
        end
        
        function obj = set_functions(obj, f)
            obj.f = f;
        end
        
        function flag = isempty(obj)
            flag = isempty(obj.samples);
        end
        
        function flag = isevaluated(obj)
            flag = ~isempty(obj.f);
        end
    end
end