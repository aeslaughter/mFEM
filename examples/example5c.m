% Example 8.1 from Bhatti, 2005 (p. 553)
%
% Syntax:
%   example5
%
% Description:
%   example5 solves a simple transient heat conduction problem, using a
%   mixed finite element mesh.
function example5c

warning('This function is not working correctly!');

% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh();
mesh.grid('Quad4',0,0.02,0,0.02,2,2);
mesh.add_element('Tri3', 1/100*[0,2; 1,2; 1,3]);
mesh.add_element('Tri3', 1/100*[1,3; 2,3; 2,4]);
mesh.grid('Quad4',0.01,0.02,0.02,0.03,1,1);
mesh.init();

% Label the boundaries
mesh.add_boundary(1, 'left','right');   % insulated (q = 0)
mesh.add_boundary(2, 'bottom');         % convective (q = h(T - Tinf))
mesh.add_boundary(3);                   % essential boundaries (all others)

% Create system
sys = System(mesh);

% Problem specifics
sys.add_constant('D', 3 * eye(2));      % thermal conductivity (W/(mC))
sys.add_constant('rho', 1600);          % density (kg/m^3)
sys.add_constant('c_p', 800);           % specific heat (J/kg)
sys.add_constant('h', 200);             % convective heat transfer coefficient (W/m^2)
sys.add_constant('T_s', 300);           % prescribed temperature on top (C)
sys.add_constant('T_inf', 50);          % ambient temperature (C)
sys.add_constant('T_0', 50);            % initial temperature (C)

% Add matrices
sys.add_matrix('M', 'rho*c_p*N''*N');
sys.add_matrix('K', 'B''*D*B');
sys.add_matrix('K_h', 'h*N''*N', 2);
sys.add_vector('f','h*T_inf*N''', 2);

% Assemble the matrices
M = sys.assemble('M');
K = sys.assemble('K'); + sys.assemble('K_h');
f = sys.assemble('f');

% Define dof indices for the essential dofs and non-essential dofs
ess = mesh.get_dof('Boundary',3); % 3,4
non = ~ess;

% Solve for the temperatures
T(:,1) = sys.get('T_0') * ones(size(f)); % initialize temperature vector
T(ess,1) = sys.get('T_s');               % apply essential boundaries

% Create the plot, at the intitial time step
mesh.plot(T);
title('t = 0');
xlabel('x');
ylabel('y');
cbar = colorbar;
set(get(cbar,'YLabel'),'String','Temperature');

% Compute residual for non-essential boundaries, the mass matrix does not
% contribute because the dT/dt = 0 on the essential boundaries
R(:,1) = f(non) - K(non,ess)*T(ess,1);

% Use a general time integration scheme
dt = 30;
theta = 0.5;
K_hat = M(non,non) + theta*dt*K(non,non);
f_K   = M(non,non) - (1-theta)*dt*K(non,non);

% Perform 10 time-steps
for t = 1:10;

    % Compute the force componenet using previous time step T
    f_hat = dt*R + f_K*T(non,t);

    % Solve for non-essential boundaries
    T(non, t+1) = K_hat\f_hat;
    
    % Set values for the essential boundary conditions the next time step 
    T(ess, t+1) = sys.get('T_s');
    
    % Plot the results
    pause(0.25);
    mesh.plot(T(:,t+1));
    title(['t = ', num2str(t)]);
end

% Display the temperatures (in same order as p.556 of Bhatti, 2005)
T'
