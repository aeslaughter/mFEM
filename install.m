function install(varargin)
%INSTALL Sets up the documentation for the mFEM library.
%
% Syntax
%   install
%   install('main')
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
    
% Define the type of install
type = 'all';
if nargin == 1 && ischar(varargin{1});
    type = varargin{1};
end

% Add the necessary paths and save 
addpath(cd);
addpath(fullfile(cd,'bin'));
addpath(fullfile(cd,'examples'));
savepath; 

% Set the pusblishing options
opt.format = 'html';
opt.showCode = false;
opt.evalCode = false;
opt.outputDir = ['doc',filesep,'html'];

% Publish main mFEM help page
if any(strcmpi(type,{'all','main'}));
    publish(fullfile(cd,'doc','main_page.m'), opt); % main mFEM page
end

% Publish the tutorial and examples
if any(strcmpi(type, {'all','examples'}));
    % List of examples page
    publish(fullfile(cd,'doc','list_of_examples.m'), opt); % main mFEM page

    % The tutorial
    opt.showCode = true;
    opt.evalCode = true;
    publish(fullfile(cd,'tutorial.m'),opt);

    % Examples
    opt.showCode = true;
    opt.evalCode = true;
    publish(fullfile(cd,'examples','example1a.m'), opt);
end

% Make the html folder available in the help browser
builddocsearchdb([cd,filesep,'doc',filesep,'html']);

% Set the mFEM root preference
setpref('MFEM_PREF','ROOT_DIR',cd);

