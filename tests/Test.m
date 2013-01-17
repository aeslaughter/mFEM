classdef Test < handle
    %TEST A class for testing code
    %   Detailed explanation goes here
    
    properties
        name;
        compare_count = 1;
    end
    
    
    methods
        function obj = Test(name)
            obj.name = name;
        end
        
        function out = func(obj, fhandle)
            try
                disp(['TEST: ', obj.name]);
                out = feval(fhandle);
            catch err
                disp(['   FAIL: ', obj.name]);
                rethrow(err);
            end
        end

        function compare(obj,mat1,mat2,msg,varargin)
            
            opt.tol = 0;
            opt = gather_user_options(opt,varargin{:});
            
            mat1 = mat2 - mat1;
            mat2 = zeros(size(mat2)) + opt.tol;
            
            if ~all(mat1 <= mat2)
                result = 'FAILED';
            else
                result = 'PASSED';
            end
            
            disp(['  ',result,' Test ', num2str(obj.compare_count),...
                ': ', msg, ' (Tol = ', num2str(opt.tol),')']);  
            obj.compare_count = obj.compare_count + 1;
        end
        
    end
    
    methods (Access = private)
    end   
    
end

