% Example 8.1 of Fish & Belytschko (2007).
%
% Syntax
%   example2b
%
% Description
%   example2b solves a simple two element heat conduction problem.
function varargout = example2b(varargin)

% Gather options
opt.debug = false;
opt = gatherUserOptions(opt, varargin{:});

% Create a FEmesh object, add the single element, and initialize it
mesh = mFEM.Mesh();
mesh.createNode([0,0; 2,0.5; 2,1; 0,1]);
mesh.createElement('Tri3',[1,2,4; 2,3,4]);
mesh.init();

% Label the boundaries
mesh.addBoundary(1, 'top');     % q = 20 boundary
mesh.addBoundary(2, 'right');   % q = 0 boundary
mesh.addBoundary(3);            % essential boundaries (all others)

% Create the System
sys = mFEM.System(mesh);
sys.addConstant('k', 5*eye(2), 's', 6, 'q_top', 20);

% Create matrices
sys.addMatrix('K', 'B''*k*B');
sys.addVector('f', 'N''*s');
sys.addVector('f', 'N''*-q_top', 'Tag', 1);

% Assemble and solve
solver = mFEM.solvers.LinearSolver(sys);
solver.addEssential('Tag',3,'Value',0);
T = solver.solve();

% Display the results
if ~opt.debug;
    T
end

% Compute the flux values for each element
% Loop through the elements
elements = mesh.getElements();
for e = 1:length(elements);
    
    % Extract the current element from the mesh object
    elem = elements(e);
    
    % Collect the local values of T
    d(:,1) = T(elem.getDof());
    
    % Compute the flux at the Gauss points
    q(:,e) = -sys.get('k')*elem.shapeDeriv([])*d;
end    

% Display the flux vectors
if ~opt.debug
    q
else
    varargout = {T,q};
end