classdef Solver < mFEM.base.handle_hide
    %SOLVER A base class for defining solvers.
    % This class serves as the backbone of the various built-in solvers
    % located in the +solvers name space. This is an abstract class that
    % requires two things be defined, the opt property which should be a
    % strcuture (see LinearSolver) and the solve method.
    %
    % See Also SYSTEM FEMESH
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
    properties (Abstract, Access = protected)
        opt;        % abstract options property
    end
    
    methods (Abstract)
        solve(obj); % abstract solver method
    end   
    
    properties %(SetAccess = private, GetAccess={?mFEM.System})
       system;  % storage for the System class instance (optional)
       mesh;    % storage for the FEmesh class instance (required)
    end
    
    properties (Access = protected)
        essential = ... % Data structure containing essential boundary info
            struct('id',{},'value',{});
    end

    methods (Access = public)
        function obj = Solver(input)
            %SOLVER Base class constructor
            %
            % Syntax
            %   Solver(system)
            %   Solver(mesh)
            %
            % Description
            %   Solver(system) constructs a Solver instanance from an
            %   existing System class instance, this will allow for
            %   the solver to automatically assemble the matrix and vectors
            %   required to solve the problem.
            %
            %   Solver(mesh) constructs a Solver instances from an existing
            %   FEmesh, this requires the user to supply a the assemble
            %   matrices and vectors required to solve the problem.
            
            if ~isa(input,'mFEM.System') && ~isa(input, 'mFEM.FEmesh');
                error('LinearSolver:LinearSolver','Input error, expected a mFEM.System or mFEM.FEmesh as input but recieved a %s',class(input));
            elseif isa(input,'mFEM.System')
                obj.system = input;
                obj.mesh = input.mesh;
            elseif isa(input, 'mFEM.FEmesh')
                obj.system = mFEM.System.empty();
                obj.mesh = input;
            end
        end
        
        function add_essential_boundary(obj, varargin)
            %ADD_ESSENTIAL_BOUNDARY label a boundary as essential
            %
            % Syntax
            %   add_essential_boundary('PropertyName', PropertyValue,...);
            %   add_essential_boundary(C1,C2,...); 
            %
            % Description
            %   add_essential_boundary('PropertyName', PropertyValue,...)
            %   adds based on the properties an essential boundary
            %   condition that will be applied to the solution
            %
            %   add_essential_boundary(C1,C2,...) same as above where each
            %   input is a seperate cell array of property pairings, this
            %   simply allows multiple boundary conditions to be assigned
            %   together.
            %
            % ADD_ESSENTIAL_BOUNDARY Property Descriptions
            %   id
            %   scalar | row vector
            %   A scalar or row vector of boundary ids (see FEMESH) that are
            %   identify the essential boundary conditions.
            %
            %   value
            %   scalar | array
            %   A scalar or vector of values that are assigned to the
            %   supplid boundary ids. If a scalar all of the dofs for the
            %   given boundary ids are assigned the value, if an array it
            %   must be have the same number of columns as the ids.
            %
            %   Component
            %       scalar | 'x' | 'y' | 'z'
            %       Returns the dof associated with the vector space or
            %       multipe dof per node elements like the Beam element.
            %
            
            % Cell input case
            if nargin > 0 && iscell(varargin{1});

                % Loop through each input, they should all be cells
                for i = 1:length(varargin);
                    
                    % Give an error if the input is not a cell
                    if ~iscell(varargin{i});
                        error('Solver:add_essential_boundary', 'Expected a cell, but recieved a %s', class(varargin{i}));
                    end
                    
                    % Add to the dof
                    obj.add_essential_boundary_private(varargin{i}{:});  
                end
                
            % Traditional input case    
            else
                obj.add_essential_boundary_private(varargin{:});  
            end
        end
        
        function set(obj, varargin)
           %SET Sets/updates the solver options
           %
           % Syntax
           %    set('PropertyName', PropertyValue, ...)
           %
           % Description
           %    set('PropertyName', PropertyValue, ...) allows user to 
           %    change the values of the options for the Solver class.
           
           obj.opt = gather_user_options(obj.opt, varargin{:});
        end
    end
    
    methods (Hidden = true, Access = protected)
       function x = get_component(obj, name, type)
           %GET_COMPONENT returns the matrix of vector ready for solving
           %
           % Syntax
           %    x = get_component(name, type)
           %
           % Description
           %    x = get_component(name, type) returns the matrix or vector
           %    associated with the name (i.e., obj.opt.(name)), in the
           %    case when a System is passed to the constructor of the
           %    Solver class and the actual text name of the matrix or
           %    vector is given the System.assemble method is called to
           %    build the component. Otherwise whatever is stored in the
           %    obj.opt.(name) structure is returned.
          
           % System case, call the assemble function        
           if ~isempty(obj.system) && ischar(obj.opt.(name));
               % Test if the component exists in the system and is the
               % correct type, if both tests pass call the assemble method
               [TF, sys_type] = obj.system.exists(obj.opt.(name));
               if TF && strcmp(sys_type, type);
                    x = obj.system.assemble(obj.opt.(name)); 
               else
                    error('Solver:get_component', 'An error occured when attempting to assemble %s, either this variable does not exist in the system or it is the wrong type.', obj.opt.(name));
               end
               
           % Generic case, the user supplied the actual matrix or vector    
           else
                x = obj.opt.(name);
           end   
       end
         
       function add_essential_boundary_private(obj, varargin)
           
            % Gather the input
            options.id = [];
            options.value = [];
            options.component = [];
            options = gather_user_options(options, varargin{:});
            
            % Test that id and value are given
            if isempty(options.id) || isempty(options.value);
                error('Solver:add_essential_boundary_private','Both the id and value properties must be set.');
            end
            
            % Append storage data structure
            idx = length(obj.essential);
            obj.essential(idx+1).id = options.id;
            obj.essential(idx+1).value = options.value;  
            obj.essential(idx+1).component = options.component;                 
       end
       
       function [u,ess] = solution_init(obj)
           
           % Initlize the solution
           u = zeros(obj.mesh.n_dof,1);
           
           % Apply the essential boundary condions
           dof = zeros(obj.mesh.n_dof, length(obj.essential),'uint32');
           for i = 1:length(obj.essential);
               dof(:,i) = obj.mesh.get_dof('Boundary', obj.essential(i).id, 'Component', obj.essential(i).component);
               u(logical(dof(:,i))) = obj.essential(i).value;
           end
           
           % Build the essential dofs
           ess = any(dof,2);
       end
   end
end

