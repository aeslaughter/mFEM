classdef Test < handle
    %TEST A class for testing code
    %   The Test class provides some basic testing functions including an
    %   evaluation method that allows this class to be added to the access
    %   list of protected methods such that they may be tested.
    %
    % Syntax
    %   T = Test('PropertyName',PropertyValue,...);
    %
    % Description
    %   T = Test('PropertyName',PropertyValue,...) defines a class for
    %   testing class and function behavior.
    %
    % TEST Property Descriptions
    %   throw
    %       true | {false}
    %       A flag for throwing the errors recieved, by default a simple
    %       fail message is shown.
    %
    %   name
    %       char
    %       The name of the test, this is printed in the output as the test
    %       progresses.
    %
    %   type
    %       {'class'} | 'filename'
    %       A flag for determining if the test being run is being peformed
    %       on a class of a test function (see bin/test.m for the latter)
    %
    %   handle
    %       function or class handle
    %       A handle to be evalued with the eval method of this class,
    %       which is useful for testing private methods.
    %
    % See Also bin\test
    %
    %----------------------------------------------------------------------
    %  mFEM: A Parallel, Object-Oriented MATLAB Finite Element Library
    %  Copyright (C) 2013 Andrew E Slaughter
    % 
    %  This program is free software: you can redistribute it and/or modify
    %  it under the terms of the GNU General Public License as published by
    %  the Free Software Foundation, either version 3 of the License, or
    %  (at your option) any later version.
    % 
    %  This program is distributed in the hope that it will be useful,
    %  but WITHOUT ANY WARRANTY; without even the implied warranty of
    %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %  GNU General Public License for more details.
    % 
    %  You should have received a copy of the GNU General Public License
    %  along with this program. If not, see <http://www.gnu.org/licenses/>.
    %
    %  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
    %----------------------------------------------------------------------
    
    % Public properties 
    properties (Access = public)
        results = [];
        handle
        err = {};
        options = struct('throw',false,'name','','type','class',...
            'handle',[]);
    end
    
    % Private properties
    properties (Hidden, Access = private)
        ticID;  
    end
    
    % Public methods
    methods
        function obj = Test(varargin)
            % Class constructor
            
            % Gather the user options
            obj.options = gatherUserOptions(obj.options,varargin{:});
            obj.handle = obj.options.handle;
            
            % Special filename input (see bin\test)
            if strcmp(obj.options.type,'filename');
                in = varargin{1};
                idx = regexp(in,'(?<=_).*');
                obj.options.name = in(idx:end);
            end
            
            % Display the message and start timer
            disp(['Testing ', obj.options.name,' Class:']);
            obj.ticID = tic;
        end

        function delete(obj)
            % Class destructor
            
            % Display the test time
            msg = sprintf(' Completed in %f sec.',toc(obj.ticID));
            disp(msg);    
        end
        
        function varargout = eval(obj,method,varargin)
            %EVAL Performs evaluation of handle
            %   By adding mFEM.Test to the access of a class it enables
            %   this function to test private class methods.
            %
            % Syntax
            %   eval(method,...)
            %   [out,...] = eval(...)
            %
            % Description
            %   eval(method,...) evalutes the method for the
            %   previously defined handle and passes the data given in the
            %   additional inputs using MATLAB's eval.
            %
            %   [out,...] = eval(...) outputs data from the 
            %   evaluation performed, the user must specify the no. of
            %   outputs to be collected.
      
            % Gather the input
            nout = nargout;
            
            % Attempt to perform evaluation
            try
                if nout == 0;
                    obj.handle.(method)(varargin{:});
                else
                    [varargout{1:nout}] = obj.handle.(method)(varargin{:});
                end

                % Success
                msg = sprintf('  %s: The %s method was evaluated',obj.options.name,method);
                obj.printResult(msg,true);

            % Evaluation failed
            catch err
                obj.err{end+1} = err;
                if obj.options.throw; rethrow(err); end
                msg = sprintf('  %s: The %s method produced an error: %s',obj.options.name,method,err.identifier);
                obj.printResult(msg,false);
            end
        end
        
        function caught(obj,err)
            %CAUGHT A simple tool for handling try/catch
            %
            % Syntax
            %   caught(err)
            %
            % Description
            %   caught(err) takes the exception caught with a try/catch
            %   statement and passes a message through the test class.
         
            msg = sprintf('  %s: Caught error: %s',obj.options.name,err.identifier);
            obj.printResult(msg,false); 
            obj.err{end+1} = err;
            if obj.options.throw;
                rethrow(err);
            end
        end
        
        function compare(obj,mat1,mat2,msg,varargin)
            %COMPARE A tool for performing comparisions
            %
            % Syntax
            %   compare(val1,val2,msg)
            %   compare(...,'PropertyName',PropertyValue,...)
            %
            % Description
            %   compare(val1,val2,msg) compares the values in val1 and val2
            %   and displays the message. the message includes a OK or FALL
            %   depending on the results of the comparision. The values may
            %   be numeric or strings.
            %
            %   compare(...,'PropertyName',PropertyValue,...) allows for
            %   additional control given the propertys below.
            %
            % compare Property Descriptions
            %   tol
            %       numeric
            %       Set the comparision tolerance, the default is zero
            %
            %   main
            %       true | {false}
            %       Flag for excluding message, used by bin\test.m

            % Gather options
            opt.tol = 0;
            opt.main = false;
            opt = gatherUserOptions(opt,varargin{:});
            
            % String compare
            if ischar(mat1) && ischar(mat2)
                bool = strcmp(mat1, mat2);
            
            % Apply tolerance comparision
            elseif opt.tol > 0;
                mat1 = mat2 - mat1;
                mat2 = zeros(size(mat2)) + opt.tol;
                bool = all(all(mat1 <= mat2));
                
            % Standard comparision
            else
                bool = all(all(mat1 == mat2));
            end
               
            % Build message
            if ~opt.main && opt.tol ~= 0;
                msg = ['  ', msg, ' (Tol = ', num2str(opt.tol, '%3.2e\n'),')'];
            elseif ~opt.main
                msg = ['  ', msg];
            end
            
            % Display message
            obj.printResult(msg,all(bool)); 
        end
    end
    
    methods (Hidden,Access = private)
        function printResult(obj, msg, TF)
            % Prints the test results
            len = 92;
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

