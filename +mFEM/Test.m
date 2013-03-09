classdef Test < handle
    %TEST A class for testing code
    %   Detailed explanation goes here
    
    properties
        results = [];
        handle
        err = {};
        options = struct('throw',false,'name','','type','class',...
            'handle',[]);
    end
    
    properties (Hidden, Access = private)
        ticID;  
    end
    
    methods
        function obj = Test(varargin)
            
            obj.options = gatherUserOptions(obj.options,varargin{:});
            obj.handle = obj.options.handle;
            
            if strcmp(obj.options.type,'filename');
                in = varargin{1};
                idx = regexp(in,'(?<=_).*');
                obj.options.name = in(idx:end);
            end
            disp(['Testing ', obj.options.name,' Class:']);
            obj.ticID = tic;
        end

        function delete(obj)
            msg = sprintf(' Completed in %f sec.',toc(obj.ticID));
            disp(msg);    
        end
        
        function varargout = eval(obj,method,input,varargin)
            
            opt.throw = obj.options.throw;
            opt.nout = 0;
            opt = gatherUserOptions(opt,varargin{:});
            
            
           try
                if opt.nout == 0;
                    obj.handle.(method)(input{:});
                else
                    [varargout{1:opt.nout}] = obj.handle.(method)(input{:});
                end
                
                msg = sprintf('  %s: The %s method was evaluated',obj.options.name,method);
                obj.printResult(msg,true);
                
           catch err
                obj.err{end+1} = err;
                if opt.throw; rethrow(err); end
                msg = sprintf('  %s: The %s method produced an error: %s',obj.options.name,method,err.identifier);
                obj.printResult(msg,false);
           end
        end
        
        
        function caught(obj,err)
           msg = sprintf('  %s: Caught error: %s',obj.options.name,err.identifier);

            obj.printResult(msg,false); 
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
            len = 128;
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

