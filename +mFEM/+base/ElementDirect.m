classdef ElementDirect < mFEM.base.ElementCore
    %ELEMENTDIRECT A base class for direct assembly type elements
    %
    
    methods (Abstract)
        stiffness(obj, varargin);
        force(obj, varargin);
    end
    
    methods  
        % Define the Truss constructor
        function obj = ElementDirect(id, nodes, varargin)
           % Class constructor; calls base class constructor
           
           % Call the base class constructor
           obj = obj@mFEM.base.ElementCore(id, nodes, varargin{:}); 
        end 
    end
end