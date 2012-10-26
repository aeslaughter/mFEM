%BUILD A tool for building the homework problems

function build(num)

% Run latex for the desired problem (returns the release directory)
release = run_latex(num);

% Define a list of m-files to copy to the release directory
switch num;
    case 1;
        mfiles = {};
end
   
% Copy the necessary MATLAB files for the problem
copyfiles(mfiles, release);














function release = run_latex(num)
%LATEX
%
% Syntax:
%   str = latex(num,soln)
%
% Variables:
%   num
%       scalar
%
%       A number indicating which homework assignment to build
%
%   soln
%       'true' | 'false'
%
%       A string indicating to inclusion of the solution
%
%

% Build format string
frmt = ['!pdflatex -jobname %s "\\input{setup.tex}\\setboolean{solution}{%s}',...
        '\\buildonly{%d}\\input{questions.tex}"'];
    
% Change the directory    
cd('latex');

% Build the loc, assignment, and soln strings
hw = sprintf('HW%d',num);
soln = [hw,'_soln'];
loc = ['..',filesep,hw];
release = [loc,filesep,'release'];

% Create the homework assignment and solution
for i = 1:2; % run twice to get references correct
    eval(sprintf(frmt, hw, 'false', num));
    eval(sprintf(frmt, soln, 'true', num));
end

% Create the HW directory if it does not exist
if ~exist(loc,'dir'); mkdir(loc); end
if ~exist(release,'dir'); mkdir(release); end

% Copy homework and solution pdfs to the correct directory directory
copyfile([hw,'.pdf'], release);
copyfile([soln,'.pdf'], loc);

% Clean up the files
delete('*.aux','*.log','*.out','*.pdf');

% Return to the original directory
cd('..');













