function err = test(varargin)
%TEST

% User options
opt.tests = {};
opt.throw = false;
opt = gather_user_options(opt,varargin{:});

% Build a list of test functions to run
func = getTestFunctions(opt.tests);

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

function func = getTestFunctions(input)

% Extract all the available tests
loc = fullfile(getpref('MFEM_PREF','ROOT_DIR'),'tests','test_*.m');
x = dir(loc);

% Use all inputs
if isempty(input);
    for i = 1:length(x);
        [~, func{i}, ~] = fileparts(x(i).name);  
    end
    
% Search function based on inputs    
else
    k = 1;
    for i = 1:length(x);
        [~,in,~] = fileparts(x(i).name);  
        idx = regexp(in,'(?<=_).*');
        current =  in(idx:end);
        TF = any(strcmp(current, input));
        if TF;
            func{k} = in;
            k = k + 1;
        end
    end
end

