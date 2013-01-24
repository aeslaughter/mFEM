% Example 9.2 of Fish & Belytschko (2007)
function [stress, strain] = example4b(varargin)

% Gather options
opt.display = true;
opt = gatherUserOptions(opt, varargin{:});

% Import the mFEM library
import mFEM.*;
  
% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Space','vector');
mesh.addElement('Quad4',[0,1; 0,0; 2,0.5; 2,1]);
mesh.init();

% Label the boundaries
mesh.addBoundary(1, 'left');    % essential boundaries
mesh.addBoundary(2, 'top');     % distributed load (t = -20)

% Create system and add necessary components
sys = System(mesh);
sys.addConstant('E', 3e7, 'v', 0.3, 't_top', [0;-20]);
sys.addConstant('D', 'E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]');
sys.addMatrix('K', 'B''*D*B');
sys.addVector('f', 'N''*t_top');

% Assemble the matrix and vector
K = sys.assemble('K'); 
f = sys.assemble('f', 'Boundary', 2);

% Assemble and solve
solver = solvers.LinearSolver(sys);
solver.addEssential('Boundary',1,'value',0);
u = solver.solve();

% Display the results
if opt.display;
    u
end

% Compute the stress and strain at the Gauss points
% Loop through the elements
for e = 1:mesh.n_elements; % (include for illustration, but not needed)
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Collect the local values of T
    d(:,1) = u(elem.getDof());
    
    % Compute the stress and strain at the Gauss points
    for i = 1:length(elem.qp);
        strain(:,i) = elem.shapeDeriv(elem.qp{i})*d;
        stress(:,i) = sys.get('D')*strain(:,i);
    end
end    

% Display the stress and strain vectors
if opt.display;
    strain, stress
end

