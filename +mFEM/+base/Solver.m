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
            %
            % Description
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
            %   Examples
            %       add_essential_boundary('id',1,'value',0);
            %       add_essential_boundary('id',[1,2],'value',0);
            %       add_essential_boundary('id',[1,2],'value',[0,10]);
            %       add_essential_boundary('id',1,'value',[10,20,30]);
            %
            %   The last example is a special case, the length of the value
            %   array must be the same as the length of the number of dofs
            %   such that the following is valid.
            %       dof = mesh.get_dof('Boundary',1);
            %       x(dof) = value; % x is the solution vector
            %   The id property must be a scalar in this case.
            %

            % Gather the input
            opt.id = [];
            opt.value = [];
            opt = gather_user_options(opt, varargin{:});
            
            % Test that id and value are given
            if isempty(opt.id) || isempty(opt.value);
                error('Solver:add_essential_boundary','Both the id and value properties must be set.');
            end
            
            % Scalar input
            if isscalar(opt.id);
                
                % Append storage data structure
                idx = length(obj.essential);
                obj.essential(idx+1).id = opt.id;
                obj.essential(idx+1).value = opt.value;
                
            % Vector input
            else
                
                % Expand values to same size as id if given as a sclar
                if isscalar(opt.value);
                    opt.value = repmat(opt.value, length(opt.id), 1);
                end
                
                % Loop through each id and assign append storage structure
                for i = 1:length(opt.id);
                    idx = length(obj.essential);
                    obj.essential(idx+1).id = opt.id(i);
                    obj.essential(idx+1).value = opt.value(i);
                end
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
    
    methods (Access = protected)
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
   end
end

