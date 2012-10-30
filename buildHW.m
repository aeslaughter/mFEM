%BUILD A tool for building the homework problems
%
% Syntax:
%   buildHW(num)
%
% Description:

function buildHW(num)

% Run latex for the desired problem (returns the release directory)
name = run_latex(num);

% Define files in the +mFEM directory that are typically included
files.mFEM = {'Element.m','Line2.m','FEmesh.m','System.m','Gauss.m'};

% Files in this list are added to the release/+mFEM directory (location
% should be the current assignment directory, e.g. assignment/HW2)
files.to_mFEM = {};

% Files in this list are added to the release/ directory (location
% should be the current assignment directory, e.g. assignment/HW2)
files.other = {};

% List of build files that are always include (location should be specified
% relative to the current directory)
files.build = {'info.xml',['doc',filesep,'html',filesep,'helptoc.xml'],'install.m'};

% Edit/append the lists of files to copy to the release directory, 
switch num;
    case 1; % HW1 does not rely on mFEM library
        files.mFEM = {};
        files.build = {};

    case 2;
        files.to_mFEM = {'Line3.m'};
        files.other = [files.other, 'tutorial.m'];
end
   
% Copy files to release directory
copy_to_release(name, files);

% Create a zip file 
loc = ['assignments',filesep,name,filesep,'release'];
cd(loc);
zip(['..',filesep,name],{'*.*',['+mFEM',filesep,'*.*']});
cd(['..',filesep,'..',filesep,'..']);

function hw = run_latex(num)
%LATEX
%
% Syntax:
%   name = latex(num,soln)
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
        '\\buildonly{%d}\\input{assignments.tex}"'];
    
% Change the directory    
cd(['assignments',filesep,'latex']);

% Build the loc, assignment, and soln strings
hw = sprintf('HW%d',num);
loc = ['..',filesep,hw,filesep];
soln = [hw,'_soln'];
soln_loc = ['..',filesep,hw,filesep,'soln'];
release = ['..',filesep,hw,filesep,'release'];

% Create the homework assignment and solution
for i = 1:2; % run twice to get references correct
    eval(sprintf(frmt, hw, 'false', num));
    eval(sprintf(frmt, soln, 'true', num));
end

% Create the HW directory if it does not exist
if ~exist(loc,'dir'); mkdir(loc); end
if ~exist(soln_loc,'dir'); mkdir(soln_loc); end
if ~exist(release,'dir'); mkdir(release); end

% Copy homework and solution pdfs to the correct directory directory
copyfile([hw,'.pdf'], release);
copyfile([soln,'.pdf'], soln_loc);

% Clean up the files
delete('*.aux','*.log','*.out','*.pdf');

% Return to the original directory
cd(['..',filesep,'..']);

function copy_to_release(name, files)
%COPY_TO_RELEASE

% Release directory
assign  = ['assignments',filesep,name,filesep];
release = [assign,'release',filesep];

% +mFEM package
for i = 1:length(files.mFEM);
    in = ['+mFEM',filesep,files.mFEM{i}];
    out = [release,'+mFEM',filesep];
    if ~exist(out,'dir'); mkdir(out); end
    copyfile(in, [out, files.mFEM{i}]);
end

% Additions to the +mFEM package
for i = 1:length(files.to_mFEM);
    in = files.to_mFEM{i};
    out = [release,'+mFEM',filesep];
    if ~exist(out,'dir'); mkdir(out); end
    [~,fname,ext] = fileparts(files.to_mFEM{i});
    copyfile(in, [out,fname,ext]);
end

% Other files
for i = 1:length(files.other);
    in = assign;
    out = [release, files.other{i}];
    copyfile(in, out);
end

% Build related files
for i = 1:length(files.build);
    in = files.build{i};
    [~,fname,ext] = fileparts(in);
    out = [release,fname,ext];
    copyfile(in, out);  
end
