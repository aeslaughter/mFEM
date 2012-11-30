%BUILDHW A tool for building the homework problems
%
% Syntax
%   buildHW(num)
%
% Description
%
%----------------------------------------------------------------------
%  mFEM: An Object-Oriented MATLAB Finite Element Library
%  Copyright (C) 2012 Andrew E Slaughter
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
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
%----------------------------------------------------------------------
    
function buildHW(num)

% Run latex for the desired problem (returns the release directory)
name = run_latex(num);

% Define files in the +mFEM directory that are typically included
files.mFEM = {'Element.m','FEmesh.m','System.m','Gauss.m','Matrix.m','handle_hide.m'};

% Define the elements to include
files.elem = {'Point.m','Line2.m'};

% Files in this list are added to the release/+mFEM/+elements directory 
% (location should be current assignment directory, e.g. assignment/HW2)
files.to_elements = {};

% Files in this list are added to the release/ directory (location
% should be the current assignment directory, e.g. assignment/HW2)
files.other = {};

% List of build files that are always include (location should be specified
% relative to the current directory; the files will be added in a relative 
% path)
files.build = {'tutorial.m','install.m', 'info.xml',...
    ['doc',filesep,'html',filesep,'helptoc.xml'],['doc',filesep,'main_page.m']};

% Edit/append the lists of files to copy to the release directory, 
switch num;
    case 1; % HW1 does not rely on mFEM library
        files.mFEM = {};
        files.build = {};

    case 2;
        files.to_elements = {'Line3.m'};
end
   
% Copy files to release directory
copy_to_release(name, files);

% Create a zip file 
loc = ['assignments',filesep,name,filesep,'release'];
if exist(loc,'dir');
    cd(loc);
    zip(['..',filesep,'..',filesep,name],{'*.*',['+mFEM',filesep,'*.*'],['+mFEM',filesep,'+elements',filesep,'*.*']});
    cd(['..',filesep,'..',filesep,'..']);
end

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

% Delete release files
if exist(release,'dir');
    rmdir(release,'s');
end

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
pause(1);
delete('*.aux','*.log','*.out','*.pdf','*.synctex');

% Return to the original directory
cd(['..',filesep,'..']);

function copy_to_release(name, files)
%COPY_TO_RELEASE

% Release directory
assign  = ['assignments', filesep, name, filesep];
release = [assign, 'release', filesep];

% +mFEM package
for i = 1:length(files.mFEM);
    in = ['+mFEM', filesep, files.mFEM{i}];
    out = [release, '+mFEM', filesep];
    if ~exist(out,'dir'); mkdir(out); end
    copyfile(in, [out, files.mFEM{i}],'f');
end

% The elements
for i = 1:length(files.elem);
    in = ['+mFEM', filesep,'+elements',filesep, files.elem{i}];
    out = [release, '+mFEM', filesep,'+elements',filesep];
    if ~exist(out,'dir'); mkdir(out); end
    copyfile(in, [out, files.elem{i}],'f');
end

% Additions to the +mFEM package
for i = 1:length(files.to_elements);
    in = [assign, files.to_elements{i}];
    out = [release, '+mFEM', filesep,'+elements',filesep];
    if ~exist(out,'dir'); mkdir(out); end
    [~,fname,ext] = fileparts(files.to_elements{i});
    copyfile(in, [out,fname,ext],'f');
end

% Other files
for i = 1:length(files.other);
    in = [assign, files.other{i}];
    out = [release, files.other{i}];
    copyfile(in, out, 'f');
end

% Build related files
for i = 1:length(files.build);
    in = files.build{i};
    [p,f,e] = fileparts(files.build{i});
    out = [release, p];
    if ~exist(out,'dir'); mkdir(out); end
    copyfile(in, [out,filesep,f,e], 'f');  
end

% Copy the bin
in = ['bin',filesep];
out = [release,'bin',filesep];
f = dir([in,'*.m']);
mkdir(out);
for i = 1:length(f);
    copyfile([in,f(i).name], out);
end