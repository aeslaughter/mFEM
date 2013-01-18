classdef LinearSolver < mFEM.solvers.base.Solver
    %LINEARSOLVER A basic linear solver.
    % This solver solve the basic Ku = f matrix equation, where u is the
    % unknown. See the class constructor for details regarding initlizeing
    % the solver correctly.
    %
    % See Also SOLVER SYSTEM FEMESH EXAMPLE1C
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
           struct('stiffness', 'K', 'force', 'f');
   end
   
   methods
       function obj = LinearSolver(input, varargin)  
           %LINEARSOLVER An automatic linear solver for Ku = f
           %
           % Syntax
           %    LinearSolver(system)
           %    LinearSolver(mesh)
           %    LinearSolver(..., 'PropertyName', PropertyValue, ...
           %
           % Description
           %    LinearSolver(system) creates the solver using an existing
           %    System class, this will exploit the system for automatic
           %    assembly of the required stiffness matrix and force vector.
           %
           %    LinearSolver(mesh) creates the 
           %    solver using an existing FEmesh class, this requires the 
           %    assembled mass and stiffness matrices and force vector to 
           %    be inputed. This may be done using the property pairing as 
           %    shown described below or by using the SET method.
           %
           %    LinearSolver(..., 'PropertyName', PropertyValue, ...)
           %    same as above but allows the user to alter the behavior
           %    using the property pairs liste below. Any of the properties
           %    may be changed after instantation using the SET method.
           %
           % LINEARSOLVER Property Descriptions
           %    stiffness
           %        char | matrix
           %        When a character it should be the name of the stiffness
           %        matrix that is desired to be used when solving the
           %        equation, Ku = f. Using a numeric array 
           %        explicitly gives the assembled matrix. The default is 
           %        the character 'K'.
           %
           %    force
           %        char | vector
           %        When a character it should be the name of the force
           %        vector that is desired to be used when solving the
           %        equation, Ku = f. Using a numeric vector 
           %        explicitly gives the assembled vector. The default is 
           %        the character 'f'.      
           
           % Call the base class constructor
           obj@mFEM.base.Solver(input)

           % Collect the inputs
           obj.opt = gather_user_options(obj.opt, varargin{:});
       end

       function u = solve(obj)
           %SOLVE Solve the linear system, Ku = f.
           %
           % Syntax
           %    u = solve()
           %
           % Description
           %    u = solve() returns the solution to the linear system of
           %    equations, Ku = f.

           % Extract/assemble the stiffness matrix
           K = obj.get_component('stiffness', 'matrix');
          
           % Extract/assemble the force vector
           f = obj.get_component('force', 'vector');
                
           % Apply essential boundary constraings to the solution
           [u,ess] = obj.apply_constraints();
   
           % Solve the equations
           u(~ess) = K(~ess,~ess)\(f(~ess) - K(ess,~ess)'*u(ess));
       end
   end
end