classdef LinearSolver < mFEM.base.Solver
   properties 
       opt = ...        % Solver options
           struct('stiffness', 'K', 'force', 'f');

       %solution;

   end
   
   methods
       function obj = LinearSolver(input, varargin)  
           
           obj@mFEM.base.Solver(input)

           obj.opt = gather_user_options(obj.opt, varargin{:});
           
       end

       function u = solve(obj)
          
           % Stiffness matrix
           K = obj.get_component('stiffness', 'matrix');
           f = obj.get_component('force', 'vector');
                
           u = zeros(size(f));
           
           dof = zeros(length(f), length(obj.essential),'uint32');
           for i = 1:length(obj.essential);
               dof(:,i) = obj.mesh.get_dof('Boundary', obj.essential(i).id);
               u(logical(dof(:,i))) = obj.essential(i).value;
           end
           
           ess = any(dof,2);
           u(~ess) = K(~ess,~ess)\(f(~ess) - K(ess,~ess)'*u(ess));

       end
   end
   

   
end