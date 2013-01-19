function test(varargin)
%RUN_TESTS

% Extract all the available tests
loc = fullfile(getpref('MFEM_PREF','ROOT_DIR'),'tests','test_*.m');
x = dir(loc);

% Use all inputs
if nargin == 0;
    for i = 1:length(x);
        [~,func{i},~] = fileparts(x(i).name);  
    end
    
% Search function based on inputs    
else
    k = 1;
    for i = 1:length(x);
        [~,in,~] = fileparts(x(i).name);  
        idx = regexp(in,'(?<=_).*');
        current =  in(idx:end);
        TF = any(strcmp(current, varargin));
        if TF;
            func{k} = in;
            k = k + 1;
        end
    end
end

func
% Perform tests
T = mFEM.Test();
for i = 1:length(func);
    t = feval(func{i});
    T.compare(all(t.results),true,t.name,'-main');
    disp(' ');
end