%EXAMPLE12a A truss example.

function example12a

% Import the mFEM library
import mFEM.*

% Create the truss structure
mesh = FEmesh();
mesh.addElement('Truss', [-1,1; 0,0]);
mesh.addElement('Truss', [0,1; 0,0]);
mesh.addElement('Truss', [1,1; 0,0]);
mesh.init();

% Label the various boundaries
mesh.addBoundary(1,{'x<=0','y>0'}); % pin connections
mesh.addBoundary(2,{'x>0','y>0'});  % roller connection
mesh.addBoundary(3,'bottom');       % applied load

% Define the known parameters
P = 10^3;
E = 10^7;
a = 10^-2;
A = [a,2*a,a];

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
    Ke = A(e)*E/L*elem.stiffness();

    % Add to the global
    K.add(Ke,elem.getDof());
end

% Extract the essential boundary conditions
ess = mesh.getDof({'Boundary',1},{'Boundary',2,'Component','y'});

% Create the stiffness matrix
K = K.init();

% Apply the external force
nat = mesh.getDof('Boundary',3,'Component','x');
f(nat) = P;

% Solve for the displacements
u(ess) = 0;
u(~ess) = K(~ess,~ess)\f(~ess);

% Plot results
mesh.plot(u,'-deform','colorbar','Magnitude of Displacement (m)');
