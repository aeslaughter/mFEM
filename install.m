%INSTALL Sets up the documentation for the mFEM library.

% Add the necessary paths

% Set the pusblishing options
opt.format = 'html';
opt.showCode = false;
opt.evalCode = false;
opt.outputDir = ['doc',filesep,'html'];

% Publish the various pages
publish(['doc',filesep,'main_page'], opt); % main mFEM page

opt.showCode = true;
opt.evalCode = true;
publish('tutorial', opt);


% Make the html folder available in the help browser
% builddocsearchdb([cd,filesep,'doc',filesep,'html']);

