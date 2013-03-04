% Example 8.1 from Bhatti, 2005 (p. 553)
%
% Syntax:
%   example5a
%
% Description:
%   example5a solves a simple transient heat conduction problem.
function varargout = example5a(varargin)
   
% Gather options
opt.debug = false;
opt = gatherUserOptions(opt, varargin{:});

% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh();
mesh.addElement('Tri3', 1/100*[0,0; 2,0; 2,4]);
mesh.addElement('Tri3', 1/100*[0,0; 2,4; 0,2]);
mesh.init();

% Label the boundaries
mesh.addBoundary(1,'left', 'right');  % insulated
mesh.addBoundary(2, 'bottom');        % convective
mesh.addBoundary(3);                  % essential

% Problem specifics
D = 3 * eye(2);         % thermal conductivity (W/(mC))
rho = 1600;             % density (kg/m^3)
c_p = 800;              % specific heat (J/kg)
h = 200;                % convective heat transfer coefficient (W/m^2)
T_s = 300;              % prescribed temperature on top (C)
T_inf = 50;             % ambient temperature (C)
T_0 = 50;               % initial temperature (C)
dt = 30;                % time step (sec)
theta = 0.5;            % time integration coefficient

% Initialize storage
M = sparse(mesh.n_dof(), mesh.n_dof());
K = sparse(mesh.n_dof(), mesh.n_dof());
f = zeros(mesh.n_dof(), 1);

% Create the mass and stiffness matrices and force vector by looping over
% elements.
for e = 1:mesh.n_elements;

    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Define short-hand function handles for the element shape functions
    % and shape function derivatives
    B = elem.shapeDeriv([]);
    N = @(xi) elem.shape(xi);

    % Initialize the local matrices and vector
    Me = zeros(elem.n_dof);     % mass matrix
    Ke = zeros(elem.n_dof);     % stiffness matrix
    fe = zeros(elem.n_dof,1);   % force vector

    % Loop over the quadrature points in the two dimensions to perform the
    % numeric integration
    for i = 1:length(elem.qp);
        Me = Me + elem.W(i)*rho*c_p*N(elem.qp{i})'*N(elem.qp{i})*elem.detJ(elem.qp{i});
        Ke = Ke + elem.W(i)*B'*D*B*elem.detJ(elem.qp{i});
    end

    % Loop through the sides of the element, if the side has the boundary
    % id of 2 (bottom), then add the convective flux term to stiffness
    % matrix and force vector using numeric integration via the quadrature 
    % points for element side.
    
    for s = 1:elem.n_sides;
        if any(elem.side(s).boundary_id == 2);
            side = elem.buildSide(s);
            N = @(xi) side.shape(xi);
            for i = 1:length(side.qp);
                d = elem.getDof('side',s,'-local'); % local dofs for side
                Ke(d,d) = Ke(d,d) + h *side.W(i)*N(side.qp{i})'*N(side.qp{i})*side.detJ(side.qp{i});
                fe(d) = fe(d) + h*T_inf*side.W(i)*N(side.qp{i})'*side.detJ(side.qp{i});              
            end
        end
    end      

    % Add local mass, stiffness, and force to global (this method is slow)
    dof = elem.getDof();
    M(dof,dof) = M(dof,dof) + Me;
    K(dof,dof) = K(dof,dof) + Ke;
    f(dof) = f(dof) + fe;
end

% Define dof indices for the essential dofs and non-essential dofs
ess = mesh.getDof('Boundary', 3); % 3,4
non = ~ess;

% Solve for the temperatures
T(:,1) = T_0 * ones(size(f));       % initialize the temperature vector
T(ess,1) = T_s;                     % apply essential boundaries

% Compute residual for non-essential boundaries, the mass matrix does not
% contribute because the dT/dt = 0 on the essential boundaries
R(:,1) = f(non) - K(non,ess)*T(ess,1);

% Use a general time integration scheme
K_hat = M(non,non) + theta*dt*K(non,non);
f_K   = M(non,non) - (1-theta)*dt*K(non,non);

% Perform 10 time-steps
for t = 1:10;

    % Compute the force componenet using previous time step T
    f_hat = dt*R + f_K*T(non,t);

    % Solve for non-essential boundaries
    T(non, t+1) = K_hat\f_hat;
    
    % Set values for the essential boundary conditions the next time step 
    T(ess, t+1) = T_s;
end

% Display the temperatures (in same order as p.556 of Bhatti, 2005)
if opt.debug; % debug, return M,K,f,T
    varargout = {M,K,f,T};
else
    ix = [1,4,2,3];
    T(ix,:)'
end
