classdef ElementCore < handle & matlab.mixin.Heterogeneous
    %ELEMENTCORE class for defining core behavior common to all elements.
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

    % Abstract Properties (must be redefined in subclass)
    properties (Abstract, Constant, GetAccess = public) 
        n_sides;      % no. of sides
        side_dof;     % array of local node ids for each side
        side_type;    % defines the type of element that defines the sides
        quad;         % Instance of Gauss quadrature class to use
    end
        
    % Public properties (read only)
    properties (SetAccess = protected, GetAccess = public)
        id = uint32([]);            % element id
        n_nodes = uint32([]);       % no. of nodes 
        n_dim = uint32([]);         % no. of spatial dimensions
        n_dof = uint32([]);         % no. of global degrees of freedom
        n_dof_node = uint32(1);     % no. of dofs per node (scalar = 1)  
        nodes = [];                 % global coordinates (no. nodes, no. dim)
        node_plot_order;            % node plotting order (only needed if nodes are out of order)
        W = [];                     % Gauss quadrature weights
        qp = [];                    % Gauss quadrature points        
        opt = ...                   % Options structure
          struct('space', 'scalar');
      
    end
    
    % Public properties (read only; except FEmesh and Element)
    properties (GetAccess = public, SetAccess = {?mFEM.FEmesh, ?mFEM.Element})
        on_boundary;                % flag if element is on a boundary
        boundary_id = uint32([]);   % list of all boundary ids for element
        subdomain = uint32([]);     % list of all subdomain ids for element
        side;                       % side info, see constructor        
        local_n_dim;                % local dimensions
        direct = false;             % a flag for using direct assembly, see Truss
    end
    
    % Protected properties
     properties (Hidden = true, Access = {?mFEM.FEmesh, ?mFEM.Element}) 
         neighbors;       % storage of nieghbor elements (see FEmesh.find_neighbors)
         global_dof = []; % Global dof for nodes of element
                          % (these are not the true dofs (except in scalar 
                          % space) as such the user should always access the
                          % dofs for an element with the get_dof() function.) 
     end    
   
    % Public Methods
    % (These methods are accessible by the user to create the element and
    % access the shape functions and other necessary parameters)
    methods (Access = public)
        function obj = ElementCore(id, nodes, varargin)
            %ELEMENTCORE Class constructor.
            %
            % This is an abstract class, it must be inherited by a subclass
            % to operate, see Line2.m for example. The following syntax and
            % descriptions apply to all subclasses unless noted otherwise
            % in the documentation for the specific element.
            %
            % Syntax
            %   Element(id, nodes)
            %   Element(id, nodes, 'PropertyName', PropertyValue, ...)
            %
            % Description
            %   Element(id, nodes) creates an element given, where id is a
            %   unique identification number for this element and nodes is 
            %   a matrix of node coordinates (global) that should be 
            %   arranged as column matrix (no. nodes x no. dims).
            %
            %   Element(id, nodes, 'PropertyName', PropertyValue, ...) 
            %   allows the user to customize the behavior of the element, 
            %   the available properties are listed below.
            %
            % Element Property Descriptions
            %   space
            %       {'scalar'} | 'vector'  | integer
            %       Allows the type of FEM space to be set: scalar sets the 
            %       number of dofs per node to 1, vector  sets it to the 
            %       no. of space dimension, and  specifing a number sets it
            %       to that value.

            % Insert required values into object properties
            obj.id = id;
            obj.nodes = nodes;
            [obj.n_nodes, obj.n_dim] = size(nodes);
            obj.local_n_dim = obj.n_dim;

            % Change dofs per node
            if nargin == 4 && strcmpi(varargin{2},'vector');
                obj.n_dof_node = obj.n_dim;  
                obj.opt.space = varargin{2};
            elseif nargin == 4 && isnumeric(varargin{2});
            	obj.n_dof_node = varargin{2};
                obj.opt.space = varargin{2};
            end   
            
            % Determine the total number of global dofs
            obj.n_dof = obj.n_nodes * obj.n_dof_node;
            
            % Intialize the side data structure
            obj.side = struct('on_boundary', true, ...
                'boundary_id', cell(obj.n_sides,1),...
                'neighbor',[], 'neighbor_side', []);
            
            % Initialize neighbor array
            obj.neighbors = feval([class(obj),'.empty']);
            
            % Get the quadrature points and weights (N/A to Points)
            if ~isempty(obj.quad);
                [obj.qp,obj.W] = obj.quad.rules('-cell');
            end
        end
               
        function size(obj)
            %SIZE Returns the length, area, or volume of the element
            %
            % Syntax
            %   L = size()
            %
            % Description 
            %   L = size() returns the size of the element, the meaning of
            %   the size depends on the element, 1D, 2D, or 3D, which
            %   should return the length, area, or volume, respectively.
            %
            % This function must be defined within the subclass, e.g., see
            % Line2. If it is not defined you will recieve an error if an
            % attempt use is made. This was not made an abstract method,
            % because in many cases this property is not needed. However,
            % the System class has a variable reserved L for it, so it
            % should be available for all Elements, even if it is an error.
            %
            % See Also
            %   LINE2
            
            % The default behavior is an error message
            error('Element:size','The size method is not defined for the %s element type',class(obj));
        end
        
        function stiffness(obj)
            %STIFFNESS Returns the complete stiffness matrix.
            %
            % Syntax
            %   Ke = stiffness()
            %
            % Description
            %   Ke = stiffness() retuns the complete stiffness matrix, this
            %   is used for special elements such as the Truss and Beam
            %   element that have a defined stiffness matrix that may be
            %   assembled directly. Use the '-direct' flag when calling the
            %   ADD_MATRIX command in the SYSTEM class.
            %
            % See Also
            %   BEAM TRUSS
            
            % The default behavior is an error message
            error('Element:stiffness','The stiffness method is not defined for the %s element type',class(obj));
        end
        
        function L = hMax(obj)
            %HMAX Return maximum length between ANY two points
            %
            % Syntax
            %   L = max();
            %
            % Description
            %   L = max() returns the maximum distance between any two
            %   points on the element.
            L = max(obj.distance());
        end
        
        function L = hMin(obj)
            %HMIN Return inimum length between ANY two points
            %
            % Syntax
            %   L = min();
            %
            % Description
            %   L = min() returns the minimum distance between any two
            %   points on the element.
            L = min(obj.distance());
        end
                  
        function varargout = getPosition(obj, varargin)
            %GETPOSITION Returns the real coordinates given xi, eta, ...
            %
            % Syntax
            %   x = ggetPosition(xi)
            %   [x,y] = getPosition(xi,eta)
            %   [x,y,z] = getPosition(xi,eta,zeta)
            %   [...] = getPosition(..., 'PropertyName', PropertyValue)
            %   xyz = getPosition(...)
            %
            % Description
            %   [...] = getPosition(...) returns the position in global
            %   system given a value(s) for the local position, the number 
            %   of outputs varies according the number of spatial 
            %   dimensions.
            %
            %   [...] = getPosition(..., 'PropertyName', PropertyValue)
            %   allow the limiting of which shape functions are used for
            %   the mapping, see EXAMPLE14b for an example.
            %
            %   xyz = getPosition(...) same as above but it returns the
            %   positions as a single array.
            %
            % GETPOSITION Property Descriptions
            %   index
            %       logical array | numeric array
            %       An array that limits what shape functions are used for
            %       the mapping, this is useful elements that have mutiple
            %       dofs per node. The value of index should be defined
            %       such that the following is valid
            %           N = obj.shape(xi,...)
            %           N = N(opt.index)
            %       
            %       Using this property overrides the default behavior,
            %       which is to use opt.index = 1:n_dof_node:end when a the
            %       no. of nodes is different from no. of shape functions,
            %       as is the case for the Beam element. For the Beam
            %       element there are 2 dofs per node and the 1 and 3 value
            %       for the shape functions may be used for maping the
            %       displacement, this behavior is the default. If you have
            %       an element that behaves otherwise, use the index
            %       property.
            
            if strcmpi(obj.opt.space,'vector');
                warning('Element:get_position','This is not tested for vector space and likely does not work correctly');
            end
            
            % Set the default options
            options.index = [];
            
            % Parse inputs
            for i = 1:length(varargin);
                TF = ischar(varargin{i});
                if TF;
                    options = gatherUserOptions(options,varargin{i:end});
                    varargin = varargin(1:i-1);
                    break;  
                end  
            end

            % The number of spatial dimensions available
            n = obj.n_dim;

            % Loop through the dimensions and return the desired position
            xyz = zeros(1,n);
            for i = 1:n
                
               % Evaluate shape functions
               N = obj.shape(varargin{:},'-scalar');

               % The index property was set, so limit the shape functions
               if ~isempty(options.index);
                   N = N(options.index);
                   
               % The no. of shape functions is not equal to no. of nodes,
               % so assume that the 1:n_dof_node:end shape functions are used to
               % map the position, this is the case for Beam elements.
               elseif length(N) > length(varargin);
                   N  = N(1:obj.n_dof_node:end);
               end

               % Map the position
               xyz(i) = N*obj.nodes(:,i);
            end
   
            % Reduce to single array if only a single output is given
            if nargout > 1;
                varargout = num2cell(xyz,1);
            else
                varargout{1} = xyz;
            end
        end
        
        function n = getNormal(obj, varargin)
            %GETNORMAL Returns the normal vector at a location
            %
            % Syntax
            %   n = getNormal(xi)
            %   n = getNormal(xi,eta)
            %
            % Description
            %   n = getNormal(...) returns the normal vector for a given
            %   local location. The number of inputs depends on the number 
            %   of space dimensions of the element. 
            %  
            % This function is still under development and is not fully
            % tested, expect changes.
            
            % The number of spatial dimensions available
            n = obj.n_dim; 
            
            % 2D: Side is defined by a line
            if n == 2;
                % Compute the tangent at the point
                n = obj.localGradBasis(varargin{:}) * obj.nodes;
                
                % Re-arrange tangent to give the normal (outward from
                % element face is positive)
                n = [n(2), -n(1)]/norm(n);

            % 3D: Side is defined by a plane    
            elseif n == 3;
                error('Element:get_normal', 'Not yet supported');
                
            % Only defined for 1D and 2D sides
            else
                error('Element:get_normal', 'Not defined for an element with %d spatial coordinates', n);
            end 
        end
                
        function side = buildSide(obj, id)
            %BUILDSIDE Build an element for the side
            %
            % Syntax
            %   side = buildSide(id)
            % 
            % Description
            %   side = buildSide(id) creates an element for the specified
            %   side, the type of element is specified in the side_type
            %   property.

            % Extract the nodes for the side
            dof = obj.side_dof(id,:);
            node = obj.nodes(dof,:);

            % Create the side element
            side = feval(['mFEM.elements.',obj.side_type], NaN, node, ...
                'Space', obj.opt.space);
            
            % Set the global dofs
            side.global_dof = obj.global_dof(dof);
        end
               
        function dof = getDof(obj, varargin)
            %GETDOF The global degrees of freedom, account for type of space
            %
            % Syntax
            %   dof = getDof()
            %   dof = getDof('PropertyName', PropertyValue,...)
            % 
            % Description           
            %   dof = getDof() returns all of the global dofs for element
            %
            %   dof = getDof('PropertyName', PropertyValue) returns the
            %   dofs for the element subject to the restrictions specified
            %   by the properties, see the descriptions below for details.  
            %
            % GET_DOF Property Descriptions
            %   Side
            %       integer
            %       Indicates that only the degrees of freedom for the
            %       specific side should be returned.
            %
            %   Local
            %       true | {false}
            %       Toggles the type of degrees of freedom to return, which
            %       is generally import for sides. For example, the
            %       following returns the local degrees of freedom for side
            %       number 1 of an element.
            %           dof = get_dof('Side',1,'local',true) or
            %           dof = get_dof('Side',1,'-local')
            %
            %   Component
            %       scalar | 'x' | 'y' | 'z'
            %       Returns the dof associated with the vector space or
            %       multipe dof per node elements like the Beam element.
            
            % Set default options and gather the options
            options.side = [];
            options.local = false;
            options.component = [];
            options = gatherUserOptions(options, varargin{:});
            
            % Extract ALL of the degrees of freedom
            if isempty(options.side)
                if options.local
                    dof = 1:obj.n_dof;
                else
                    dof = obj.global_dof;
                end
            
            % Extract dofs for the specified side    
            elseif isnumeric(options.side) && options.side <= obj.n_sides;
                s = options.side;
                if options.local;
                    dof = obj.side_dof(s,:);
                else
                    idx = obj.side_dof(s,:);
                    dof = obj.global_dof(idx);
                end
            else
                error('Element:get_dof','Input error.');
            end
            
            % Extract boundaries for certain component of vector space
            if ~isempty(options.component);
                
                % Convert string to numeric
                if ischar(options.component);
                    switch options.component;
                        case 'x'; ix = 1;
                        case 'y'; ix = 2;
                        case 'z'; ix = 3;
                    end  
                else
                    ix = options.component;
                end

                % Test the dimensionality
                if ix > obj.n_dof_node;
                    error('Element:get_dof','Desired vector component does not exist.');
                end
                
                % Extract the dofs
                ix = ix:obj.n_dof_node:length(dof);
                dof = dof(ix);
            end   
            
            % Non-scalar FE space
            if obj.n_dof_node > 1;
                dof = obj.transformDof(dof);
            end
        end  
    end
    
    methods (Hidden = true, Access = {?mFEM.FEmesh})        
        function D = transformDof(obj, d)
            %TRANSFROMDOF Converts the dofs for vector element space
            %
            % Syntax
            %   D = transformDof(d)
            %   
            % Description
            %   D = transformDof(d) converts the scalar degrees of freedom
            %   for to vector based degrees of freedom. For example,
            %   inputing d = [1,3] returns D = [1,2,5,6].
            
            D = transformDof(d,obj.n_dof_node); % call func in bin
        end
        
        function d = distance(obj, varargin)
            %DISTANCE Compute distances between all points
            %
            % Syntax
            %   d = distance()
            %   d = distance(node)
            %
            % Description
            %   d = distance() computes the distance between all nodes, so
            %   is returns N! distances, where N is the number of nodes
            %   minus one.
            %
            %   d = distance(node) distance of all nodes from the node
            %   given.
           
            if nargin == 1;
                d = zeros(factorial(obj.n_nodes-1),1);
                k = 0;
                for i = 1:obj.n_nodes;
                    for j = i+1:obj.n_nodes;
                        k = k + 1;
                        d(k) = norm(obj.nodes(i,:) - obj.nodes(j,:));
                    end
                end
                
            elseif nargin == 2
                node = varargin{1};
                d = zeros(obj.n_nodes,1);
                for i = 1:obj.n_nodes;
                    d(i) = norm(obj.nodes(i,:) - node);
                end 
            end
        end
    end
end