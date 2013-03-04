classdef Test < handle
    %TEST A class for testing code
    %   Detailed explanation goes here
    
    properties
        name;
        results = [];
    end
    
    
    methods
        function obj = Test(varargin)
            if nargin == 1;
                in = varargin{1};
                idx = regexp(in,'(?<=_).*');
                obj.name = in(idx:end);
                disp(['Testing ', obj.name]);
            end
        end

        function str = compare(obj,mat1,mat2,msg,varargin)
            
            opt.tol = 0;
            opt.main = false;
            opt.showmax = false;
            opt = gatherUserOptions(opt,varargin{:});
            
            if ischar(mat1) && ischar(mat2)
                bool = strcmp(mat1, mat2);
            
            elseif opt.tol > 0;
                mat1 = mat2 - mat1;
                mat2 = zeros(size(mat2)) + opt.tol;
                bool = all(all(mat1 <= mat2));
            else
                bool = all(all(mat1 == mat2));
            end
               
            if ~opt.main && opt.tol ~= 0;
                msg = ['  ', msg, ' (Tol = ', num2str(opt.tol, '%3.2e\n'),')'];
            elseif ~opt.main
                msg = ['  ', msg];
            end
            
            obj.printResult(msg,all(bool));
            
        end
    end
    
    methods (Access = private)
        function printResult(obj, msg, TF)
            len = 82;
            str(1:len) = '.';
            if length(msg) > len - 6;
                error('Test:buildstr', 'Messages must be limited to %d characters.', len-8);
            end
            
            str(1:length(msg)) = msg;
            if TF;
                result = 'OK';
            else
                result = 'FAIL';
            end
                
            str(end-length(result)+1:end) = result;
            disp(str);
            obj.results(end+1) = TF;
        end
    end       
end

