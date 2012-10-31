% A transient heat transfer example
%
% Syntax:
%   example6a
%   example6a('PropertyName', PropertyValue)
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
% 
%    The number of elements in the x and y directions, the default is 32.
%
% Element
%   {'Quad4'} | 'Tri3' | 'Tri6'
%
%   Specifies the type of element for the mesh
%
% Method
%   {'normal'} | 'alt'
%
%   Inidicates the type of sparse matrix assembly to utilize, the 'alt'
%   method is the index method that is faster for large matrices. For
%   example, for a 100 x 100 grid the normal assembly took 25.6 sec. and the
%   alternative method 23.8 sec.

function example6a(varargin)

% Import the mFEM library
import mFEM.*;

% Set the default options and apply the user defined options
opt.n = 32;
opt.element = 'Quad4';
opt.method = 'normal';
opt = gather_user_options(opt,varargin{:});

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh();
mesh.grid(opt.element,0,1,0,1,opt.n,opt.n);
mesh.init();

% Label the boundaries
mesh.add_boundary(1); % essential boundaries (all)

% Create Gauss objects for performing integration on the element and
% elements sides.
q_elem = Gauss(2,'quad');
[qp, W] = q_elem.rules();

% Problem specifics
D = 1 / (2*pi^2);           % thermal conductivity
theta = 0.5;                % numerical intergration parameter
dt = 0.1;                   % time-step

% Initialize storage

ticID = tmessage('Matrix assembly...');
if strcmpi(opt.method,'alt');
    I = NaN(mesh.n_elements * mesh.n_dim^2,1); % (guess)
    J = I;
    Mij = I;
    Kij = I;
else
    M = sparse(mesh.n_dof(), mesh.n_dof());
    K = sparse(mesh.n_dof(), mesh.n_dof());
end

% Create mass and stiffness matrices by looping over elements
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

    % Insert current values into global matrix using one of two methods
    if strcmpi(opt.method,'alt');
        % Get the local degrees of freedom for this element
        dof = (1:elem.n_dof)';

        % Compute indices for inserting into sparse matrix i,j,s vectors
        m = numel(Me);
        idx = m*(e-1)+1 : m*(e);

        % Build the i,j components for the sparse matrix creation
        i = repmat(dof, length(dof),1);
        j = sort(i);
        
        % Get the global degrees of freedom for this element
        dof = elem.get_dof();
        I(idx) = dof(i);
        J(idx) = dof(j);

        % Add the local mass and stiffness matrix to the sparse matrix values
         Mij(idx) = reshape(Me, numel(Me), 1);
         Kij(idx) = reshape(Ke, numel(Ke), 1);
    else
        % Add local mass, stiffness, and force to global (this method is slow)
        dof = elem.get_dof();
        M(dof,dof) = M(dof,dof) + Me;
        K(dof,dof) = K(dof,dof) + Ke;
    end
end

% If the alternative method of assembly is used, the sparse matrices must
% be created using the I,J,Mij,Kij vectors
if strcmpi(opt.method,'alt');
    % Assemble sparse matrices
    M = sparse(I,J,Mij);
    K = sparse(I,J,Kij);
end
tmessage(ticID);

% Define dof indices for the essential dofs and non-essential dofs
non = mesh.get_dof(1,'ne'); 
ess = mesh.get_dof(1);      

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
    pause(0.25);
    mesh.plot(T);
    title(['t = ', num2str(t)]);
end

% Clean up
delete(mesh);
