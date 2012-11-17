%EXAMPLE13b A truss example.
% Beer and Johnston, Ex. 6.2, p. 227

function example13b

% Cal subfunction to create the truss structure
mesh = build_truss(); 

% Label the various boundaries
mesh.add_boundary(1,{'x==8','y==10'},{'x==24','y==10'}); % 28 kip load
mesh.add_boundary(2,{'x==40','y==10'}); % 16 kip load
mesh.add_boundary(3,{'x==0','y==0'});   % pin connections
mesh.add_boundary(4,{'x==32','y==0'});  % roller connection

% Define the known parameters
P1 = -28e3;
P2 = 16e3;

% Define the system
sys = mFEM.System(mesh);
sys.add_constant('A', 1, 'E', 10^7);
sys.add_matrix('K', 'A*E/L*Ke');

% Initilize the force vector
f = zeros(mesh.n_dof,1);

% Extract the essential boundary conditions
ess = mesh.get_dof({'Boundary',3},{'Boundary',4,'Component','y'});

% Create the stiffness matrix
K = sys.assemble('K');

% Apply the external forces
n1 = mesh.get_dof('Boundary',1,'Component','y');
n2 = mesh.get_dof('Boundary',2,'Component','x');
f(n1) = P1;
f(n2) = P2;

% Solve for the displacements
u(ess,1) = 0;
u(~ess,1) = K(~ess,~ess)\f(~ess);

% Solve/display for the reactions on essential boundaries
r = K*u - f;
r(ess)

% Plot results
mesh.plot('-new','NodeLabels',false,'ElementLabels',false);
mesh.plot(u,'-deform','Component',2,'Scale',3,'colorbar','Magnitude of Displacement (in)');
end

function mesh = build_truss()
    %BUILD_TRUSS creates the mesh for the example truss

    mesh = mFEM.FEmesh('Element','Truss');

    % Positions of the nodes
    x = [0,8,16,24,32,40];
    y = 10;

    % Create the first 4 boxes with cross brace
    for i = 2:length(x)-1;
        mesh.add_element([x(i-1),y; x(i-1),0]);
        mesh.add_element([x(i-1),0; x(i),0]);
        mesh.add_element([x(i),y; x(i-1),y]);
        mesh.add_element([x(i),y; x(i-1),0]);   
    end

    % Add the final triangle
    mesh.add_element([x(end-1),y; x(end-1),0]);
    mesh.add_element([x(end-1),0; x(end),y]);
    mesh.add_element([x(end-1),y; x(end),y]);
    mesh.init();
end
