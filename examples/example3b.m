% Example 8.2 of Fish & Belytschko (2007).
%
% Syntax:
%   example3
%
% Description:
%   example3 solves a simple single element heat conduction problem.
function varargout = example3b(varargin)

% Gather options
opt.debug = false;
opt = gatherUserOptions(opt, varargin{:});

% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh();
mesh.addElement('Quad4',[0,1; 0,0; 2,0.5; 2,1]);
mesh.init();

% Label the boundaries
mesh.addBoundary(1, 'top');     % q = 20 boundary
mesh.addBoundary(2, 'right');   % q = 0 boundary
mesh.addBoundary(3);       % essential boundaries (all others)

% Create the System
sys = System(mesh);
sys.addConstant('k', 5*eye(2), 'b', 6, 'q_top', 20);

% Create matrices
sys.addMatrix('K', 'B''*k*B');
sys.addVector('f', 'N''*b');
sys.addVector('f', 'N''*-q_top', 'Boundary', 1);

% Assemble and solve
solver = solvers.LinearSolver(sys);
solver.addEssential('boundary',3,'value',0);
T = solver.solve();

if ~opt.debug; 
    T
end

% Loop through the elements
for e = 1:mesh.n_elements; % (include for illustration, but not needed)
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Collect the local values of T
    d(:,1) = T(elem.getDof());
    
    % Compute the flux at the Gauss points
    for i = 1:length(elem.qp);
        q(:,i) = -sys.get('k')*elem.shapeDeriv(elem.qp{i})*d;
    end
end    

% Display the flux vectors
if ~opt.debug;
    q
else
    varargout = {T,q};
end