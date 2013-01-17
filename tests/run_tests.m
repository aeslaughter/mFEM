function run_tests
%RUN_TESTS

loc = fullfile(getpref('MFEM_PREF','ROOT_DIR'),'tests','test_*.m');
x = dir(loc);
disp(' ');
for i = 1:length(x);
    [~,func,~] = fileparts(x(i).name);
    feval(func);
    disp(' ');
end