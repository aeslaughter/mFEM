classdef System < mFEM.registry.ConstantRegistry
    %SYSTEM A class for automatic assembly of finite element equations.
    % This class allows the specification of the finite element equations
    % for matrices and vectors as strings, the assembly is handled
    % automatically.
    %
    % Example
    %   mesh = FEmesh()
    %   mesh.grid('Line2',0,1,10)
    %   mesh.init()
    %
    %   sys = System(mesh)
    %   sys.add_constant('k',10);
    %   sys.add_matrix('K','B''*k*B');
    %   K = sys.assemble('K');
    %   full(K)
    %
    % See Also FEMESH
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
    
    properties(Access = public);
        time = 0;                   
    end
    
    properties (SetAccess = private, GetAccess = public)
        mesh = mFEM.FEmesh.empty;       % mesh object
        options = ...                   % struct of default user options
            struct('time', false);
    end
    
%     properties(Access = private)
%         m_reg;
%         c_reg;
%         f_reg;
%     end

    methods (Access = public)
        function obj = System(mesh, varargin)
            %SYSTEM Class constructor.    
            %
            % Syntax
            %   sys = System(mesh)
            %   sys = System(mesh, 'PropertyName', PropertyValue, ...)
            %
            % Description
            %   sys = System(mesh) creates a System object based on the
            %   supplied FEmesh object.
            %
            %   sys = System(mesh, 'PropertyName', PropertyValue) same as 
            %   above but allows the options to be changes
            %
            % System Properties Descriptions
            %   Time
            %       true | {false}
            %       A toggle for displaying the matrix and vector 
            %       assembly time.
            %   
            % See Also FEMESH
            
            obj = obj@mFEM.registry.ConstantRegistry();
            
            % Parse the user-defined options
            obj.options = gather_user_options(obj.options, varargin{:});
            
            % Store the mesh object
            obj.mesh = mesh;
            
%             obj.m_reg = mFEM.registry.MatrixKernelRegistry(mesh);
%             obj.c_reg = mFEM.registry.ConstantKernelRegistry();
%             obj.f_reg = mFEM.registry.FunctionKernelRegistry();
        end
%           
%         function add_constant(obj, varargin)
%             %ADD_CONSTANT Adds constant(s) variables to the system.
%             %
%             % Syntax
%             %   add_constant('ConstantName', ConstantValue, ...)
%             %
%             % Description
%             %   add_constant('ConstantName', ConstantValue, ...) adds
%             %   a single or multiple constants to the System, the
%             %   'ConstantName' may be any string identifier and the value
%             %   may be a numeric (scalar or matrix) or a string that may be
%             %   evaluated using MATLAB's eval function.
%             %
%             % Examples
%             %   sys.add_constant('k', 10, 'r', 2);
%             %   sys.add_cosntant('D','k^2');
%             
%             obj.c_reg.add(varargin{:});
%         end
%         
%         function add_function(obj, varargin)
%             %ADD_CONSTANT Adds constant(s) variables to the system.
%             %
%             % Syntax
%             %   add_function('FunctionName', FunctionHandle, ...)
%             %
%             % Description
%             %   add_function('FunctionName', FunctionHandle) adds
%             %   a functions to the System, the 'FunctionName' may be 
%             %   any string identifier and the FunctionHandle may be a
%             %   valid MATLAB function handle or character string, see
%             %   details below.
%             %
%             % Description of Function Handle Input
%             %   If the FunctionHandle property is a MATLAB function handle
%             %   it must be in the following form:
%             %
%             %       FunctionHandle = @(elem, x, t) ...
%             %       
%             %   elem = An Element class object.
%             %   x = The current position, where x(1) is the x-direction, 
%             %   x(2) is the y-direction, and x(3) is the z-direction
%             %   t = The current time, which is automatically passed to the
%             %   function from the t variable of the System class itself.
%             %
%             % Description of Text Input
%             %   It is also possible to specify a function as a text string
%             %   such that the following is valid:
%             %       fhandle = str2func(['@(elem,x,t)', <Text Input Here>]);
%             %
%             % Examples
%             %   sys.add_function('k',@(elem,x,t) x(1)^2 + x(2)^2);
%             %   sys.add_fucntion('p','x(1)^3');
%             
%             obj.f_reg.add(varargin{:});
%         end
% 
%         function add_matrix(obj, varargin)  
%             %ADD_MATRIX Create a sparse finite element matrix for assembly
%             %
%             % Syntax
%             %   add_matrix('MatrixName', MatrixEqn, ...)
%             %   add_matrix('MatrixName', MatrixEqn, ... 'PropertyName', PropertyValue, ...)
%             %
%             % Description
%             %   add_matrix('MatrixName', MatrixEqn) adds a matrix or
%             %   matrices the associated equations for finite element
%             %   assembly. The MatrixEqn is a string that is valid MATLAB
%             %   (i.e., it works with the eval function) that gives the 
%             %   integrands of the finite element equations at the element
%             %   level. The variable N is predefined as the element shape 
%             %   functions vector and B is pre-defined as grad(N). Any
%             %   constant added with the ADD_CONSTANT method may be used.
%             %
%             %   add_matrix('MatrixName', MatrixEqn, 'PropertyName',
%             %       PropertyValue) same as above, but the equation is 
%             %   limited according to the supplied options. For example, the
%             %   'Boundary' property limits the equation application to the 
%             %   boundaries with the same id (see FEMESH.ADD_BOUNDARY)
%             %
%             % ADD_MATRIX Property Descriptions
%             %   Boundary
%             %       numeric
%             %       Limits the application of the supplied equation to the 
%             %       boundaries with the id, see FEMESH.ADD_BOUNDARY
%             %
%             %   Subdomain
%             %       numeric
%             %       Limits the application of the supplied equation to the 
%             %       elements on the subdomain, see FEMESH.ADD_SUBDOMAIN
%             %
%             %   Component
%             %       numeric
%             %       Limits the application of the supplied equation to the
%             %       specified component for vector unknowns.
%             %
%             % Examples
%             %   sys.add_matrix('M','N''*N');
%             %   sys.add_matrix('K','B''*B','Boundary',1);
%   
%             obj.add_matrix_private('matrix', varargin{:});
%         end
% 
%         function add_vector(obj, varargin)  
%             %ADD_VECTOR Create a finite element vector for assembly
%             %
%             % Syntax
%             %   add_vector('VectorName', VectorEqn, ...)
%             %   add_vector('VectorName', NumericVector, ...)
%             %   add_vector(..., 'PropertyName', PropertyValue)
%             %
%             % Description
%             %   add_vector('VectorName', VectorEqn) adds
%             %   a vector and the associated equation for finite element
%             %   assembly. The VectorEqn is a string that is valid MATLAB
%             %   (i.e., it works with the eval function) that gives the 
%             %   integrands of the finite element equations at the element
%             %   level. The variable N is predefined as the element shape 
%             %   functions vector and B is pre-defined as grad(N). Any
%             %   constant added with the ADD_CONSTANT method may be used.
%             %
%             %   add_vector('VectorName', NumericVector) directly insert a
%             %   numeric vector.
%             %
%             %   add_vector('VectorName', VectorEqn, 'PropertyName',
%             %       PropertyValue) same as above, but the equation is 
%             %   limited according to the supplied options. For example, the
%             %   'Boundary' property limits the equation application to the 
%             %   boundaries with the same id (see FEMESH.ADD_BOUNDARY)
%             %
%             % ADD_VECTOR Property Descriptions
%             %   Boundary
%             %       numeric
%             %       Limits the application of the supplied equation to the 
%             %       boundaries with the id, see FEMESH.ADD_BOUNDARY
%             %
%             %   Subdomain
%             %       numeric
%             %       Limits the application of the supplied equation to the 
%             %       elements on the subdomain, see FEMESH.ADD_SUBDOMAIN
%             %
%             %   OverWrite
%             %       true | {false}
%             %       If true the vector is overwritten if the same name is
%             %       encounter, otherwise the new vector is added to any
%             %       existing.
%             %
%             %   Dof
%             %       numeric vector
%             %       Limits the application of the numeric input to specific
%             %       dofs.
%             %
%             % Example
%             %   sys.add_vector('f','N''*b'); % b is a constant
% 
%             obj.add_matrix_private('vector', varargin{:});
%         end
%         
%         function output = point_value(obj, name, elem, x)
%             %POINT_VALUE extract the value for a vector at a point x
%             
%             [type,idx] = obj.locate(name);
%             if ~strcmp(type,'vector');
%                 error('System:point_value', 'This method only works with vectors');
%             end
%             
%             % Degrees of freedom
%             dof = elem.get_dof();
% 
%             % Extract local vector
%             u = obj.vec(idx).vector.get_local(dof);
% 
%             % Shape functions at point x
%             N = elem.shape(x{:});
% 
%             % Output the value
%             output = N*u;
%         end
%         
%         function output = point_gradient(obj, name, elem, x)
%             %POINT_GRADIENT extract the gradient of a vector at a point x
%             
%             [type,idx] = obj.locate(name);
%             if ~strcmp(type,'vector');
%                 error('System:point_value', 'This method only works with vectors');
%             end
%             
%             % Degrees of freedom
%             dof = elem.get_dof();
% 
%             % Extract local vector
%             u = obj.vec(idx).vector.get_local(dof);
% 
%             % Shape functions at point x
%             B = elem.shape_deriv(x{:});
% 
%             % Output the value
%             output = B*u;
%         end
%        
%         function X = get(obj, name)
%             %GET Returns the value for the specified name.
%             %
%             % Syntax
%             %   X = get(name)
%             %
%             % Description
%             %   X = get(name) returns the value of constant, vector, or 
%             %   matrix given the string name. In the case of of matrices
%             %   the unassemble Matrix class object is returned.
%             
%             % Determine the type and location in the structure
%             [type, idx] = obj.locate(name);
%             
%             % Throw and error if the name was not found
%             if isempty(idx);
%                 error('System:locate', 'The entity with name %s was not found.', name);
%             end   
%             
%             % Extract the value based on the type
%             switch type;
%                 case 'matrix';   X = obj.mat(idx).matrix;
%                 case 'vector';   X = obj.vec(idx).vector;
%                 case 'constant'; X = obj.const(idx).value;
%                 case 'function'; X = obj.func(idx).handle;
%             end
%         end
%         
%         function TF = exists(obj, name)
%            %EXISTS Returns true if name exists as a type in System
%            %
%            % Syntax
%            %    TF = exists(name)
%            %    [TF, type] = exists(name)
%            %
%            % Description
%            %    TF = exists(name) returns a true value if a constant,
%            %    matrix, vector, or function exists in the System with the
%            %    name given.
%            %    [TF, type] = exists(name) same as above but also returns
%            %    the type ('mat', 'vec', 'const', or 'func').
%             
%            kern = obj.get_kernel(name);
%            TF = false;
%            if ~isempty(kern); TF = true; end
%            
%         end
%         
%         function X = assemble(obj, name, varargin)
%             %ASSEMBLE Assembles matrix or vector given by name
%             %
%             % Syntax
%             %   X = assemble(name)
%             %
%             % Description
%             %   X = assemble(name) assembles the finite element matrix or
%             %   vector. In the case of a matrix the sparse matrix is
%             %   returned, for vectors a tradional MATLAB array is given.
% 
%             %kern = obj.find_kernel(name);
%             X = obj.m_reg.assemble(name, varargin{:});
%         end
%         
%         function kern = get_kernel(obj, name)
% 
%            reg = {obj.m_reg, obj.c_reg, obj.f_reg};
% 
%            kern = [];
%            for i = 1:length(reg);
%                kern = reg{i}.locate(name);
%                if ~isempty(kern); return; end
%            end
% 
%         end  
%        end 
%     
%     methods (Hidden = true, Access = private)
% 
%         function add_matrix_private(obj, type, varargin)
%             
%             opt.boundary = [];
%             opt.subdomain = [];
%             opt.component = [];
%             [opt, unknown] = gather_user_options(opt, varargin{:});
% 
%             % Location of last input flag
%             n = length(unknown) - 1;
%             
%             % Loop through each name and store in the const property
%             for i = 1:2:n;
%                 name = unknown{i};
%                 kern = unknown{i+1};
%                 
%                 if isa(kern, 'mFEM.kernels.base.MatrixKernel');
%                     obj.m_reg.add(name, kern);
%                 elseif ischar(kern);
%                     kern = mFEM.kernels.AutoKernel(obj.mesh, name, kern,...
%                         'Constants', obj.c_reg,...
%                         'Functions', obj.f_reg,...
%                         'boundary', opt.boundary,...
%                         'subdomain', opt.subdomain,...
%                         'component', opt.component,...
%                         'type', type);
%                     obj.m_reg.add(name, kern);
%                 else
%                     error('System:add_matrix_private','Input of type %s is not allowed, you must supply a valid character string or a MatrixKernel object', class(kern));
%                 end
%             end  
%         end
    end
end