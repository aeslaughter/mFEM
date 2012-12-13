% Example 8.1 from Bhatti, 2005 (p. 553)
%
% Syntax:
%   example5b
%
% Description:
%   example5b solves a simple transient heat conduction problem.
function example5b
   
% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Element','Tri3');
mesh.add_element(1/100*[0,0; 2,0; 2,4]);
mesh.add_element(1/100*[0,0; 2,4; 0,2]);
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
sys.add_matrix('K', 'h*N''*N', 'Boundary', 2);
sys.add_vector('f','h*T_inf*N''', 'Boundary', 2);

% Create solver
solver = solvers.TransientLinearSolver(sys, 'dt', 30);

% Add essential boundary
solver.add_essential_boundary('id',3,'value','T_s');

% Initialize the temperatures
 T(:,1) = sys.get('T_0') * ones(mesh.n_dof,1); % initialize temperature vector
% T(ess,1) = sys.get('T_s');               % apply essential boundaries
solver.init(T(:,1));

% % Compute residual for non-essential boundaries, the mass matrix does not
% % contribute because the dT/dt = 0 on the essential boundaries
% R(:,1) = f(non) - K(non,ess)*T(ess,1);
% 
% % Use a general time integration scheme
% dt = 30;
% theta = 0.5;

% K_hat = M(non,non) + theta*dt*K(non,non);
% f_K   = M(non,non) - (1-theta)*dt*K(non,non);

% Perform 10 time-steps
for t = 1:10;

    % Compute the force componenet using previous time step T
%     f_hat = dt*R + f_K*T(non,t);
% 
%     % Solve for non-essential boundaries
%     T(non, t+1) = K_hat\f_hat;
    
    % Set values for the essential boundary conditions the next time step 
    T(:, t+1) = solver.solve();
end

% Display the temperatures (in same order as p.556 of Bhatti, 2005)
ix = [1,4,2,3];
T(ix,:)'
