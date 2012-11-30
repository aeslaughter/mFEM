%INSTALL Sets up the documentation for the mFEM library.
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
    
% Add the necessary paths and save 
addpath(cd);
addpath('bin');
savepath; 

% Set the pusblishing options
opt.format = 'html';
opt.showCode = false;
opt.evalCode = false;
opt.outputDir = ['doc',filesep,'html'];

% Publish main mFEM help page
publish([cd,filesep,'doc',filesep,'main_page.m'], opt); % main mFEM page

% Publish the tutorial
opt.showCode = true;
opt.evalCode = true;
publish('tutorial', opt);

% Make the html folder available in the help browser
builddocsearchdb([cd,filesep,'doc',filesep,'html']);

