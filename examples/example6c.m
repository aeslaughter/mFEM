% A transient heat transfer example, using Matrix class
%
% Syntax:
%   example6b
%   example6b('PropertyName', PropertyValue)
%
% Description:
%   example6 solves a simple transient heat conduction problem, with the
%   default settings.
%
%   example6a('PropertyName', PropertyValue) allows the user to customize
%   the behavior of this example using the property pairs listed below.
%
% Example6a Property Descriptions
%
% N
%   scalar
%   The number of elements in the x and y directions, the default is 32.
%
% Element
%   {'Quad4'} | 'Tri3' | 'Tri6'
%   Specifies the type of element for the mesh
%

function example9b(varargin)

% Set the default options and apply the user defined options
opt.n = 32;
opt.element = 'Quad4';
opt = gather_user_options(opt,varargin{:});

% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, add the single element, and initialize it
tic;
mesh = FEmesh('Element',opt.element);
mesh.grid(0,1,0,1,opt.n,opt.n);
mesh.init();

% Display time for mesh creation
disp(['Mesh generation time: ', num2str(toc), ' sec.']);

% Label the boundaries
mesh.add_boundary(1); % essential boundaries (all)

% Create Gauss objects for performing integration on the element and
% elements sides.
q_elem = Gauss('Order',2,'Type','quad');
[qp, W] = q_elem.rules();

% Problem specifics
D = 1 / (2*pi^2);           % thermal conductivity
theta = 0.5;                % numerical intergration parameter
dt = 0.1;                   % time-step

% Initialize storage
M = Matrix(mesh.n_dof);
K = Matrix(mesh.n_dof);

% Create mass and stiffness matrices by looping over elements
tic;
for e = 1:mesh.n_elements;

    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Define short-hand function handles for the element shape functions
    % and shape function derivatives
    B = @(xi,eta) elem.shape_deriv(xi, eta);
    N = @(xi,eta) elem.shape(xi,eta);

    % Initialize the local matrices and vector
    Me = zeros(elem.n_dof);     % mass matrix
    Ke = zeros(elem.n_dof);     % stiffness matrix
    
    % Loop over the quadrature points in the two dimensions to perform the
    % numeric integration
    for i = 1:length(qp);
        for j = 1:length(qp);
            Me = Me + W(i)*N(qp(i),qp(j))'*N(qp(i),qp(j))*elem.detJ(qp(i),qp(j));
            Ke = Ke + W(i)*B(qp(i),qp(j))'*D*B(qp(i),qp(j))*elem.detJ(qp(i),qp(j));
        end
    end
    
    % Add the local to global
    dof = elem.get_dof();
    M.add_matrix(Me, dof);
    K.add_matrix(Ke, dof);
end

% Initialize the sparse matrices (these method deletes the Matrix objects)
K = K.init();
M = M.init();

% Print assembly time
disp(['Matrix assembly time: ', num2str(toc), ' sec.']);

% Define dof indices for the essential dofs and non-essential dofs 
ess = mesh.get_dof('Boundary',1);      
non = ~ess;

% Initialize the temperatures
T_exact = @(x,y,t) exp(-t)*sin(pi*x).*sin(pi*y);
nodes = mesh.get_nodes();
T = T_exact(nodes(:,1),nodes(:,2),0);

% Collect the node positions for applying the essential boundary conditions
x = nodes(:,1);
y = nodes(:,2);

% Plot the initial condition
figure; hold on;
mesh.plot(T);
title('t = 0');
xlabel('x');
ylabel('y');
cbar = colorbar;
set(get(cbar,'YLabel'),'String','Temperature');

% Compute residual for non-essential boundaries, the mass matrix does not
% contribute because the dT/dt = 0 on the essential boundaries. This
% problem also does not have a force term.
R(:,1) = - K(non,ess)*T(ess);

% Use a general time integration scheme
K_hat = M(non,non) + theta*dt*K(non,non);
f_K   = M(non,non) - (1-theta)*dt*K(non,non);

% Perform 10 time-steps
for t = dt:dt:1;

    % Compute the force componenet using previous time step T
    f_hat = dt*R + f_K*T(non);

    % Solve for non-essential boundaries
    T(non) = K_hat\f_hat;
    
    % Set values for the essential boundary conditions the next time step 
    T(ess) = T_exact(x(ess), y(ess), t);

    % Plot the results
    pause(0.5);
    mesh.plot(T);
    title(['t = ', num2str(t)]);
end
