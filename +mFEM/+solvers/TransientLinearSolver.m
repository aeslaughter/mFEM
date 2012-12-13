classdef TransientLinearSolver < mFEM.base.Solver
    %TRANSIENTLINEARSOLVER A basic transient linear solver.
    % This solver solve the basic Mdu/dt + Ku = f matrix equation, 
    % where u is the unknown. See the class constructor for details 
    % regarding intializing the solver correctly.
    %
    % See Also LINEARSOLVER SOLVER SYSTEM FEMESH EXAMPLE1C
    %
    %----------------------------------------------------------------------
    %  mFEM: An Object-Oriented MATLAB Finite Element Library
    %  Copyright (C) 2012 Andrew E Slaughter
    % 
    %  This program is free software: you can redistribute it and/or modify
    %  it under the terms of the GNU General Public License as published by
    %  the Free Software Foundation, either version 3 of the License, or
    %  (at your option) any later version.
    % 
    %  This program is distributed in the hope that it will be useful,
    %  but WITHOUT ANY WARRANTY; without even the implied warranty of
    %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %  GNU General Public License for more details.
    % 
    %  You should have received a copy of the GNU General Public License
    %  along with this program.  If not, see <http://www.gnu.org/licenses/>.
    %
    %  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
    %----------------------------------------------------------------------
   properties (Access = protected)
       opt = ...        % Solver options
           struct('mass', 'M', 'stiffness', 'K', 'force', 'f', ...
           'theta', 0.5, 'dt', [], 'assemblemass', false, ...
           'assemblestiffness', false, 'assembleforce', false, ...
           'assembleall', false);

       initialized = false;     % flag of initlization state
       K;                       % stiffness matrix
       M;                       % mass matrix
       f;                       % force vector
       f_old;                   % force vector from previous time step
       u_old;                   % solution from previous time step
   end
   
   methods
       function obj = TransientLinearSolver(input, varargin)  
           %TRANSIENTLINEARSOLVER An automatic linear solver for Ku = f
           %
           % Syntax
           %    TransientLinearSolver(system)
           %    TransientLinearSolver(mesh)
           %    TransientLinearSolver(..., 'PropertyName', PropertyValue, ...
           %
           % Description
           %    TransientLinearSolver(system) creates the solver using an existing
           %    System class, this will exploit the system for automatic
           %    assembly of the required stiffness matrix and force vector.
           %
           %    TransientLinearSolver(mesh) creates the 
           %    solver using an existing FEmesh class, this requires the 
           %    assembled mass and stiffness matrices and force vector to 
           %    be inputed. This may be done using the property pairing as 
           %    shown described below or by using the SET method.
           %
           %    TransientLinearSolver(..., 'PropertyName', PropertyValue, ...)
           %    same as above but allows the user to alter the behavior
           %    using the property pairs liste below. Any of the properties
           %    may be changed after instantation using the SET method.
           %
           % TRANSIENTLINEARSOLVER Property Descriptions
           %    Note: All of the following properties may be
           %    defined/re-defined using the SET method.
           %
           %    stiffness
           %        char | matrix
           %        When a character it should be the name of the stiffness
           %        matrix that is desired to be used when solving the
           %        equation, Mdu/dt + Ku = f. Using a numeric array 
           %        explicitly gives the assembled matrix. The default is 
           %        the character 'K'.
           %
           %    force
           %        char | vector
           %        When a character it should be the name of the force
           %        vector that is desired to be used when solving the
           %        equation, Mdu/dt + Ku = f. Using a numeric vector 
           %        explicitly gives the assembled vector. The default is 
           %        the character 'f'.          
           %
           %    mass
           %        char | matrix
           %        When a character it should be the name of the mass
           %        matrix that is desired to be used when solving the
           %        equation, Mdu/dt + Ku = f. Using a numeric array 
           %        explicitly gives the assembled matrix. The default is 
           %        the character 'M'.   
           %
           %    dt
           %        scalar
           %        Defines the time step to use. This may be changed at
           %        any point using the SET method.
           %            theta = 0 results in forward difference integration
           %            theta = 1/2 (default) is the Crank-Nicolson scheme
           %            theta = 2/3 is the Galerkin method
           %            theta = 1 is the backward difference scheme
           %
           %    theta
           %        scalar
           %        Defines numerical integration coefficient. This value
           %        must be between 0 and 1. 
           %
           %    AssembleStiffness
           %        true | {false}
           %        Toggles the automatic assembly of the stiffness matrix,
           %        set this to false if the matrix does not change with
           %        time.
           %
           %    AssembleMass
           %        true | {false}
           %        Toggles the automatic assembly of the mass matrix,
           %        set this to false if the matrix does not change with
           %        time.
           %
           %    AssembleForce
           %        true | {false}
           %        Toggles the automatic assembly of the force vector,
           %        set this to false if the matrix does not change with
           %        time.
           %
           %    AssembleAll
           %        true | {false}
           %        Toggles the automatic assembly for all three of the
           %        above components.

           % Call the base class constructor
           obj@mFEM.base.Solver(input)
           
           % Collect the inputs
           obj.opt = gather_user_options(obj.opt, varargin{:});
       end

       function init(obj, x)
           %INIT Initilizes the TransientLinearSolver
           %
           % Syntax
           %    init(u)
           %
           % Description
           %    init(u) initlizes the transient solution to the vector u
           %    supplied by the user
           
           % Test size
           if length(x) ~= obj.mesh.n_dof;
               error('TransientLinearSolver:init', 'Expected a vector of length %d for the solution, but recieved on of length %d.', obj.mesh.n_dof, length(x));
           end
           
           % Store the old solution
           obj.u_old = x; 
           
           % Apply the boundary constraints
           obj.u_old = obj.apply_constraints(obj.u_old);

           % Get the initial force vector
           obj.f_old = obj.get_component('force', 'vector');
          
           % Set the initlized flag
           obj.initialized = true;
       end
       
       function u = solve(obj)
            %SOLVE Solve the transient system, Mdu/dt + Ku = f.
            %
            % Syntax
            %    u = solve()
            %
            % Description
            %    u = solve() returns the solution to the transient
            %    system of equations, Mdu/dt + Ku = f.

            % Produce error if not initilized
            if ~obj.initialized;
                error('TransientLinearSolver:solve','The solver must be initlized using the init method.');
            end
            
            % Produce error if numerical integration value is out of range
            if obj.opt.theta < 0 || obj.opt.theta > 1;
                error('TransientLinearSolver:solve','The value of theta must be between 0 and 1.');
            end
            
            % Extract/assemble the stiffness matrix
            if isempty(obj.K) || obj.opt.assemblestiffness || obj.opt.assembleall;
                obj.K = obj.get_component('stiffness', 'matrix');
            end

            % Extract/assemble the force vector
            if isempty(obj.f) || obj.opt.assembleforce || obj.opt.assembleall;
                obj.f = obj.get_component('force', 'vector');
            end

            % Extract/assemble the mass matrix
            if isempty(obj.M) || obj.opt.assemblemass || obj.opt.assembleall;
               obj.M = obj.get_component('mass', 'matrix');
            end

            % Apply boundary conditions to old and new solution
            [u,ess] = obj.apply_constraints();

            % Numerical constants
            theta = obj.opt.theta; % numerical intergration parameter
            dt = obj.opt.dt;       % time-step

            % Use a general time integration scheme
            K_hat = obj.M + theta*dt*obj.K;
            f_hat = dt*(theta*obj.f + (1-theta)*obj.f_old) + ...
                (obj.M - (1-theta)*dt*obj.K)*obj.u_old;

            % Solve for the unknowns
            u(~ess) = K_hat(~ess,~ess)\(f_hat(~ess) - K_hat(ess,~ess)'*u(ess));

            % Update stored values
            obj.u_old = u;
            obj.f_old = obj.f;
       end
   end
end