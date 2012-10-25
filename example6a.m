% A transient heat transfer example
%
% Syntax:
%   example6
%   example6(N)
%   example6(..,'alt');
%
% Description:
%   example6 solves a simple transient heat conduction problem, using the
%   index sparse assembly method
%
%   example6(N) specifies the number of division in the mesh, the default
%   is 32, which creats a 32x32 element mesh.
%
%   example6(..,'alt') uses the faster, but more complicated sparse matrix
%   creation method
%
% Sprase Matrix Creation:
% This example imploys a more efficient method for creating sparse
% matrices, see "doc sparse". Using the I,J method for the sparse matrix
% improved the matrix assembly time. As shown below, the increase was small
% for the 100 x 100 element mesh tested. Considering that the index method
% is more intuitive the time savings may not be worth it.
%
% Assembly time for I,J method: 22.02 s
% Assembly time for index method: 25.61 s

function example6a(varargin)

% Determine the number of elements to divide the mesh into
N = 32;
if nargin >= 1 && isnumeric(varargin{1});
    N = ceil(varargin{1});
end

% Specify the type of assembly to use
alt = false;
if nargin >= 1 && ischar(varargin{end}) && strcmpi(varargin{end},'alt');
    alt = true;
end

% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, add the single element, and initialize it
tic;
mesh = FEmesh('Quad4');
mesh.grid(0,1,0,1,N,N);

% Display time for mesh creation
disp(['Mesh generation time: ', num2str(toc), ' sec.']);

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
if alt;
    I = NaN(mesh.n_elements * mesh.n_dim^2,1); % (guess)
    J = I;
    Mij = I;
    Kij = I;
else
    M = sparse(mesh.n_dof(), mesh.n_dof());
    K = sparse(mesh.n_dof(), mesh.n_dof());
end

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
    
    % Insert current values into global matrix using one of two methods
    if alt;
        % Get the global degrees of freedom for this element
        dof = elem.get_dof();

        % Compute indices for inserting into sparse matrix i,j,s vectors
        m = numel(Me);
        idx = m*(e-1)+1 : m*(e);

        % Build the i,j components for the sparse matrix creation
        I(idx) = repmat(dof, length(dof), 1);
        J(idx) = sort(repmat(dof, length(dof), 1));

        % Add the local mass and stiffness matrix to the sparse matrix values
        Mij(idx) = reshape(Me, numel(Me), 1);
        Kij(idx) = reshape(Ke, numel(Ke), 1);

    else
        % Add local mass, stiffness, and force to global (this method is slow)
        dof = elem.get_dof();
        M(dof,dof) =  M(dof,dof) + Me;
        K(dof,dof) = K(dof,dof) + Ke;
    end
end

% If the alternative method of assembly is used, the sparse matrices must
% be created using the I,J,Mij,Kij vectors
if alt;
    % Assemble sparse matrices
    M = sparse(I,J,Mij);
    K = sparse(I,J,Kij);
end

% Print assembly time
disp(['Matrix assembly time: ', num2str(toc), ' sec.']);

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
    pause(0.5);
    mesh.plot(T);
    title(['t = ', num2str(t)]);
end

% Clean up
delete(mesh);
