%% mFEM Tutorial: 1D Boundary Value Problem
% A simple 1D heat transfer problem that details the use of the mFEM
% library. 

%% Problem Statement
% Consider a bar with a uniformily distributed heat source (s) of 5 W/m.
% The bar has a uniform cross sectional area (A) of 0.1 m^2 and thermal
% conductivity (k) of 2 W/(Cm). The length of the bar is 4 m.
% The boundary conditions are T(0) = 0 C and q(4) = 5 W/m^2 as shown in the 
% homework assignment.

%% Function Call
% The function is executed from the command line simply by running the
% |tutorial.m| from the MATLAB command line, as follows.
%
%   tutorial;
%
function tutorial


%% Setting up the mFEM Library
% In order to utilize the mFEM library, the classes must be made available
% by using the IMPORT command.
import mFEM.*

%% Create a FE space
% The FEmesh class controls the geometric aspects of the FEM including
% minimal mesh generation and boundary indentification. The following
% commands create a simple 1D FEM mesh for the current problem.
%
% An instance of the |FEmesh| class is created. This
% class contains a method named |grid| wich is called and creates a 1D grid
% of 2-noded elements from 0 to 4. This grid is divided into two elements,
% for details on the |FEmesh| refer to the associated documentation using
% the MATLAB help browser (this is available for all classes in the +mFEM 
% directory).
%
%   doc mFEM.FEmesh
%
% Finally, after the mesh is created it must be initilized, this is done
% using the |init| method. As shown, the |init| method will print
% calculation times for creating the mesh and other required calculations.
mesh = FEmesh();
mesh.grid('Line2',0,4,2);
mesh.init();

%% Adding Boundary Identification
% The FEmesh class allows for the tagging of boundaries using keywords. The
% following commands tag the left and right boundaries with ids 1 and 2
% respectively. In 2D there are also 'top' and 'bottom' keywords available.
% It is also possible to mark any unmarked boundaries. For example, the
% following code would produce the same results, but the |add_boundary|
% method implemented without a keyword applies the id to any unmarked
% boundary, which for this case is only the right hand side.
%
%   mesh.add_boundary('left', 1);
%   mesh.add_boundary(2);
%
mesh.add_boundary('left', 1);   % T = 0 boundary (essential)    
mesh.add_boundary('right', 2);  % q = 5 boundary  

%% Define Finite Element Equations
% The finite element matrix and vector creation as well as their assembly
% is handled automatically using the System class, which requires the
% FEmesh object as input. It is possible to create multiple systems from
% the same FEmesh object.
sys = System(mesh);

%%
% Constants are added using the |add_constant| method, they may be entered
% as pairs as done here or with individual calls to the method or any
% combination.
sys.add_constant('k',2,'A',0.1,'b',5,'q_bar',5); 

%%
% It is also possible to use constants in the definition of
% other constants using strings. For example, the following code develops the
% constintutive matrix given the modulus of elasticity E and posion ratio
% v. Notice that the constants added may be non scalar.
%
%   sys.add_constant('E',3e7,'v',0.3);
%   sys.add_constant('D', 'E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]');
%
% Executing the above is equivelement the following code on the MATLAB
% command line.
%
%   E = 3e7;
%   v = 0.3;
%   D = eval('E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]');
%

%%
% Matrices and vectors are created using the |add_matrix| or the
% |add_vector| methods. For both, the syntax is the same. The first
% agrument is simply a name for later use. The second is a string that
% gives the element level calculations to perform. This expression should
% be given using standard MATLAB syntax that may be executed using the
% |feval| function (see |doc feval|).
%
% The element shape functions are denoted using |N| and the shape function
% derivatives using |B|. The integration and transformation of local to
% global coordinates is handled automatically (e.g., the Jacobian and Gauss
% integration weight functions do not need to be included). 
%
% Any constants there were previously defined may be included in the
% function string.
%
% The following adds the stiffness matrix to the system, in equation form
%
% $$K_e = \int_{\Omega_e} B^T k A B d\Omega.$
%
sys.add_matrix('K', 'B''*k*A*B');

%%
% In similar fashion, the force vector is defined in equation form as
%
% $$f_{s_e} = \int_{\Omega_e} N^T b d\Omega$
%
% and
%
% $$f_{q_e} = \int_{\Gamma_{q_e}} - N^T \bar{q} A d\Gamma.$
%
% Notice that an additional argument is given for the boundary term, this
% is the boundary id defined previously that this equation should be
% applied.
sys.add_vector('f_s', 'N''*b');
sys.add_vector('f_q', '-q_bar*A*N''', 2);

%% Assembly The Stiffness Matrix and Force Vector
% Assembly is handled autmatically by the System and it is performed using
% the |assemble| method.
K = sys.assemble('K');
f = sys.assemble('f_s') + sys.assemble('f_q');

%% Solve for the Temperature
% The solution is done manually, the first step is to extract the
% degree of freedom indices for the essential and non-essential boundaries.
% The FEmesh class |get_dof| method provides this functionality. First, the
% non-essential boundaries are extracted, where the 'ne' flag indicates not
% equal, so in this case the |get_dof| method returns the global degrees of
% freedom not associated with boundary id 1. 
non = mesh.get_dof(1,'ne'); % 2,3
ess = mesh.get_dof(1);      % 1

%%
% Solve for the temperatures first by initializing the solution vector and
% setting the values of the vector on the essential boundary to zero, as
% prescibed in the problem.
T = zeros(size(f));         % initialize the temperature vector
T(ess) = 0;                 % apply essential boundary condtions

%%
% The unknown temperature are then solved using the non-essential degrees
% of freedom and MATLAB built-in solver (see |doc \|).
T(non) = K(non,non)\f(non); % solve for T on the non-essential boundaries

%% Post Processing
% At this point the solution is complete, it took only 18 lines of code to
% solve the problem. However, the result is better understood if graphed.
% The FEmesh class as built in plotting capabilities.
%
% To create a plot of the results simply call the |plot| method with the
% temperature data as an input. The |plot| method is simply a wrapper to
% the |FEplot| function, see |doc FEplot| for details on this program.
mesh.plot(T);

%%
% The exact solution for this problem is known, this can easily be added to
% the plot created using standard MATLAB commands. Additionally, labels and
% legends may be added.
x = 0:0.1:4;
Tex = -12.5*x.^2 + 97.5*x;
hold on;
plot(x,Tex,'k-o','LineWidth',1);
legend({'FEM','Exact'},'location','best');
xlabel('x (m)','interpreter','tex');
ylabel('Temperature (\circC)','interpreter','tex');

