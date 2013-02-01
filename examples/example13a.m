%EXAMPLE13a A truss example.
% Beer and Johnston, Ex. 6.2, p. 227

function varargout = example13a(varargin)

% Define user options
opt.debug = false;
opt = gatherUserOptions(opt, varargin{:});

% Import the mFEM library
import mFEM.*

% Create the truss structure
mesh = buildTruss();

% Label the various boundaries
mesh.addBoundary(1,{'x==8','y==10'},{'x==24','y==10'}); % 28 kip load
mesh.addBoundary(2,{'x==40','y==10'}); % 16 kip load
mesh.addBoundary(3,{'x==0','y==0'});   % pin connections
mesh.addBoundary(4,{'x==32','y==0'});  % roller connection

% Define the known parameters
P1 = -28e3;
P2 = 16e3;
E = 10^7;
A = 1;

% Initilize the stiffness matrix and force vector
K = Matrix(mesh);
f = zeros(mesh.n_dof,1);

% Loop through the elements and append global stiffness
for e = 1:mesh.n_elements;
    
    % The current element
    elem = mesh.element(e);
    
    % The element length
    L = elem.size();

    % Element stiffness matrix
    Ke = A*E/L*elem.stiffness();

    % Add to the global
    K.add(Ke,elem.getDof());
end

% Extract the essential boundary conditions
ess = mesh.getDof({'Boundary',3},{'Boundary',4,'Component','y'});

% Create the stiffness matrix
K = K.init();

% Apply the external forces
n1 = mesh.getDof('Boundary',1,'Component','y');
n2 = mesh.getDof('Boundary',2,'Component','x');
f(n1) = P1;
f(n2) = P2;

% Solve for the displacements
u(ess,1) = 0;
u(~ess,1) = K(~ess,~ess)\f(~ess);

% Solve/display for the reactions on essential boundaries
r = K*u - f;
r = r(ess);

% Display Results
if ~opt.debug;
    disp(['Bx = ', num2str(r(1)/1000), ' kips']);
    disp(['By = ', num2str(r(2)/1000), ' kips']);
    disp(['J = ', num2str(r(3)/1000), ' kips']);

    % Plot results
    mesh.plot([],'-new','NodeLabels',false,'ElementLabels',false);
    mesh.plot(u,'-deform','Component',2,'Scale',3,'colorbar','Magnitude of Displacement (in)');
else
    varargout = {u,r};
end
end

function mesh = buildTruss()
    %BUILDTRUSS creates the mesh for the example truss

    mesh = mFEM.FEmesh();

    % Positions of the nodes
    x = [0,8,16,24,32,40];
    y = 10;

    % Create the first 4 boxes with cross brace
    for i = 2:length(x)-1;
        mesh.addElement('Truss',[x(i-1),y; x(i-1),0]);
        mesh.addElement('Truss',[x(i-1),0; x(i),0]);
        mesh.addElement('Truss',[x(i),y; x(i-1),y]);
        mesh.addElement('Truss',[x(i),y; x(i-1),0]);   
    end

    % Add the final triangle
    mesh.addElement('Truss',[x(end-1),y; x(end-1),0]);
    mesh.addElement('Truss',[x(end-1),0; x(end),y]);
    mesh.addElement('Truss',[x(end-1),y; x(end),y]);
    mesh.init();
end