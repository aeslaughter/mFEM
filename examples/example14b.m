%EXAMPLE14b Beam element
% Fish and Belytschko (2009), Ex. 10.1, p. 262
%
% See Also EXAMPLE14a

function example14b

% Load the mFEM library
import mFEM.*

% Create the 2-element mesh
n = 100;
mesh = FEmesh('Element','Beam');
mesh.grid(0,12,n*3); % using ensures that nodes always are at nodes
mesh.init();

% Add flags for the boundary elements
mesh.add_boundary(1, 'left');    % fixed connection
mesh.add_boundary(2, 'right');   % prescribed force and moment

% Add flag for the distributed load and point loads
mesh.add_subdomain(10, 'x<8');    % distributed load, b
mesh.add_subdomain(11, 'x==4');   % point load, P1
mesh.add_subdomain(12, 'x==8');   % point load, P2

% Difine the paramters for the problem
b = -1;         % body force
P1 = -10;       % load at x = 4
P2 = 5;         % load at x = 8
EI = 10^4;      % modulus and moment of intertia
c = [-20;20];   % presecribed force and moment on boundary

% Initialize global stiffness matrix and force vector
K = Matrix(mesh);
f = zeros(mesh.n_dof,1);

% Loop through the elements
for e = 1:mesh.n_elements;
    
    % Current element
    elem = mesh.element(e);
    
    % Initilize local K and f
    Ke = zeros(elem.n_dof);    
    fe = zeros(elem.n_dof,1);    
    
    % Create functions for N and B
    B = @(xi) elem.shape_deriv(xi);
    N = @(xi) elem.shape(xi);
    
    % Extract quadrature points
    [qp,w] = elem.quad.rules();

    % Perform quadrature for stiffness matrix
    for i = 1:length(qp);
        Ke = Ke + w(i)*EI*B(qp(i))'*B(qp(i))*elem.detJ(qp(i));
    end

    % Account for the body force, only on subdomain
    if any(elem.subdomain == 10);
        for i = 1:length(qp);  
            fe = fe + w(i)*N(qp(i))'*b*elem.detJ(qp(i));
        end
    end
    
    % Add the local matrix and vector to the global
    dof = elem.get_dof();
    K.add_matrix(Ke, dof);
    f(dof) = f(dof) + fe;
end

% Create the stiffness matrix
K = K.init();

% Apply point loads
dof1 = mesh.get_dof('Subdomain', 11, 'Component', 1);
f(dof1) = f(dof1) + P1;

dof2 = mesh.get_dof('Subdomain', 12, 'Component', 1);
f(dof2) = f(dof2) + P2;

% Apply the natural boundary conditions
dof3 = mesh.get_dof('Boundary', 2);
f(dof3) = f(dof3) + c;

% Solve for the unknowns
ess = mesh.get_dof('Boundary' , 1);
u(ess) = 0;
u(~ess) = K(~ess,~ess)\f(~ess);

% Plot the results
h = subplot(4,1,1);
mesh.plot(u, 'Axes', h, 'Component', 1, '-deform', 'patch', {'EdgeColor','k','LineWidth',2});
xlabel('x (m)');
ylabel('Displacement (m)');

h = subplot(4,1,2);
mesh.plot(u, 'Axes', h, 'Component', 2, '-deform', 'patch', {'EdgeColor','k','LineWidth',2});
xlabel('x (m)');
ylabel('Rotation (rad.)');

% Compute the moment and shear
s = []; sx = []; m = []; mx = [];
for e = 1:mesh.n_elements;
    
    % Current element
    elem = mesh.element(e);
    dof = elem.get_dof();  
    
    % Create functions for N and B
    B = @(xi) elem.shape_deriv(xi);
    D = elem.dNdx3();

    % Shear
    s(end+1) = EI*D*u(dof)';
    s(end+1) = s(end);
    sx(end+1) = elem.nodes(1);
    sx(end+1) = elem.nodes(2);
    
    % Moments at quass points
    [qp,~] = elem.quad.rules();
    for i = 1:length(qp);
        m(end+1) = EI*B(qp(i))*u(dof)';
        mx(end+1) = elem.get_position(qp(i));
    end
end

% Plot the shear force
subplot(4,1,3);
plot(sx,s,'b-', 'LineWidth', 2);
xlabel('x (m)');
ylabel('Shear (N)');

% Plot the moment
subplot(4,1,4);
plot(mx,m,'b-', 'LineWidth', 2);
xlabel('x (m)');
ylabel('Moment (Nm)');
