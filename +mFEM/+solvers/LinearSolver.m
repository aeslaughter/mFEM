classdef LinearSolver < mFEM.base.Solver
   properties 
       matrix = 'K';
       vector = 'f';
       %solution;

   end
   
   methods
       function obj = LinearSolver(sys)  
           obj@mFEM.base.Solver(sys)
       end

       function u = solve(obj)
          
           K = obj.system.assemble(obj.matrix);
           f = obj.system.assemble(obj.vector);
           u = zeros(size(f));
           
           dof = zeros(length(f), length(obj.essential),'uint32');
           for i = 1:length(obj.essential);
               dof(:,i) = obj.system.mesh.get_dof('Boundary', obj.essential(i).id);
               u(logical(dof(:,i))) = obj.essential(i).value;
           end
           
           ess = any(dof,2);
           u(~ess) = K(~ess,~ess)\(f(~ess) - K(ess,~ess)'*u(ess));

       end
   end
end