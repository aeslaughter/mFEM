function err = test(varargin)
%TEST
%
% Syntax
%   test;
%   test(TestName)
%   test({TestName1, TestName2, ...})
%   test('examples');
%   test(..., 'PropertyName', PropertyValue, ...)

% User options
opt.tests = {};
opt.throw = false;

% Get a list of tests to perform
[func, options] = getTestFunctions(varargin{:});

% Gather the user options
opt = gatherUserOptions(opt, options{:});

% Perform tests
T = mFEM.Test();
err = struct([]);
for i = 1:length(func);
    try
        t = feval(func{i});
        T.compare(all(t.results),true,t.name,'-main');
        disp(' ');
        
    catch e
        if isempty(err);
            err = e;
        else
            err(end+1) = e;
        end
        msg = [func{i},'.m failed to execute, see error structure ', num2str(length(err))];
        T.compare(false,true,msg,'-main');
        disp(' ');
        if opt.throw;
            rethrow(e);
       end
    end
end

function [output, options] = getTestFunctions(varargin)

% Account for special case of 'examples'
examples = false;
if nargin >= 1 && ischar(varargin{1}) && strcmpi('examples',varargin{1});
    loc = fullfile(getpref('MFEM_PREF','ROOT_DIR'), 'tests', 'test_example*.m');
    x = dir(loc);  
    examples = true;
    
% Extract all the available tests
else
    loc = fullfile(getpref('MFEM_PREF','ROOT_DIR'), 'tests', 'test_*.m');
    x = dir(loc);
end

% Build a cell array of all functions available
for i = 1:length(x);
    [~, func{i}, ~] = fileparts(x(i).name);  
end

% Return the complete list if no inputs are given
if nargin == 0;
   options = {};
   output = func;
   return;
elseif examples
    output = func;
    options = varargin(2:end);
    return;
end

% Extract the first optional argument
if ischar(varargin{1}); 
    input1 = varargin(1);
else
    input1 = varargin{1};
end

% Loop through all inputs to the first optional argument
output = {};
for i = 1:length(input1)
   t = ['test_',input1{i}];
   if any(strcmp(t, func));
       output{end+1} = t;
   end
end
    
% Extract remaining options
if isempty(output)
    options = varargin;
    output = func;
else
    options = varargin(2:end);
end