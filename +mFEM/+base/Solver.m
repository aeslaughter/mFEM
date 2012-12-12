classdef Solver < handle
    %SOLVER A base class for defining solvers.
    
    properties %(GetAccess={?mFEM.System})
       system;
       essential = struct('id',{},'value',{});
    end
    
    methods (Abstract, Access = public)
        solve(obj);
    end
    
    
    methods (Access = public)
        function obj = Solver(sys)
            if ~isa(sys,'mFEM.System');
                error('LinearSolver:LinearSolver','Input error, expected a mFEM.System as input but recieved a %s',class(sys));
            end
            obj.system = sys;
        end
        
        function essential_boundary(obj, id, value)

            idx = length(obj.essential);
            obj.essential(idx+1).id = id;
            obj.essential(idx+1).value = value;

        end
    end
    
end

