%% mFEM: An object-oriented MATLAB finite element library.
% 
%% Getting started with the mFEM Toolbox
% A 1D heat transfer boundary value problem is used to demonstrate the use
% of the mFEM library. To run the example, type the following on the 
% command line.
%
%   tutorial;
%
% The source code for this example is documented in a step-by-step manner:
% <tutorial.html *mFEM Tutorial Documentation*>.
% 
%% Examples
% The example problems fall into one of two types: manual and automatic
% assembly. The manual assembly examples demonstrate how to loop through
% the elements and assemble the stiffness matrix and force vector directly.
% The automatic assembly rely on the |System| class to perform the
% assembly. The are listed in the contents menu on the left:
%
% <list_of_examples.html *List of Examples*>.
% 
%% Documentation
% The documentation for all functions, classes, and members is accessible
% using MATLAB's existing |help| and |doc| functions.
%
% To see a linked list of all classes in the |+mFEM| library type the
% following, which will open in the help browser a list of all the classes
% in the library. Simply click the class to see the associated help.
%
%   doc mFEM
%
% All of the classes for the mFEM library are are located in the |+mFEM| 
% directory, which is a MATLAB package. The documention for any of the
% classes in this directory may be called using the |doc| command. For 
% example, the following opens the documentation for the FEmesh class.
%
%   doc mFEM.FEmesh
% 
% Within the |+mFEM| directory is another package that contains all of the
% element classes (which are subclasses of Element), to access the
% documentation for these subclasses, you must also include the element
% package.
%
%   doc mFEM.elements.Line2
%
