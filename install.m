%INSTALL Sets up the documentation for the mFEM library.

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

