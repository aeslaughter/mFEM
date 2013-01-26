classdef Solver < handle
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
        options;        % abstract options property
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
            struct('boundary',{},'value',{},'component',{});
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
        
        function addEssential(obj, varargin)
            %ADDESSENTIAL label a boundary as essential
            %
            % Syntax
            %   addEssential('PropertyName', PropertyValue,...);
            %   addEssential(C1,C2,...); 
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
            %   boundary
            %       scalar | row vector
            %       A scalar or row vector of boundary ids (see FEMESH) that are
            %       identify the essential boundary conditions.
            %
            %   value
            %       scalar | array | char
            %       A scalar or vector of values that are assigned to the
            %       supplid boundary ids. If a scalar all of the dofs for the
            %       given boundary ids are assigned the value, if an array it
            %       must be have the same number of columns as the ids. It is
            %       also possible to use a char, where the value is extracted
            %       from the System class using the get method.
            %
            %   Component
            %       scalar | 'x' | 'y' | 'z'
            %       Returns the dof associated with the vector space or
            %       multipe dof per node elements like the Beam element.
            %
            %   Clear
            %       true | {false}
            %       Removes any existing essential boundaries, useful for
            %       changing the conditions with time. If clear is used in
            %       the C1, C2 style of input it will clear all boundaries
            %       including those defined previously in the same call.
            %
            
            % Cell input case
            if nargin > 0 && iscell(varargin{1});

                % Loop through each input, they should all be cells
                for i = 1:length(varargin);
                    
                    % Give an error if the input is not a cell
                    if ~iscell(varargin{i});
                        error('Solver:addEssentialBoundary', 'Expected a cell, but recieved a %s', class(varargin{i}));
                    end
                    
                    % Add to the dof
                    obj.addEssentialPrivate(varargin{i}{:});  
                end
                
            % Traditional input case    
            else
                obj.addEssentialPrivate(varargin{:});  
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
           
           obj.opt = gatherUserOptions(obj.opt, varargin{:});
        end
    end
    
    methods (Hidden = true, Access = protected)          
        function x = getComponent(obj, name)
           %GETCOMPONENT returns the matrix or vector ready for solving
           %
           % Syntax
           %    x = getComponent(name)
           %
           % Description
           %    x = getComponent(name) returns the matrix or vector
           %    associated with the name (i.e., obj.opt.(name)), in the
           %    case when a System is passed to the constructor of the
           %    Solver class and the actual text name of the matrix or
           %    vector is given the System.assemble method is called to
           %    build the component. Otherwise whatever is stored in the
           %    obj.opt.(name) structure is returned.
           
            % System case, call the assemble function        
            if ~isempty(obj.system) && ischar(obj.options.(name));
                % Test if the component exists in the system and is the
                % correct type, if both tests pass call the assemble method
                TF = obj.system.exists(obj.options.(name));
                
                if TF;% && strcmp(sys_type);
                    x = obj.system.assemble(obj.options.(name),'-zero'); 
                else
                    error('Solver:getComponent', 'The %s was not found in the system when attempting to assemble.', obj.opt.(name));
                end

            % Generic case, the user supplied the actual matrix or vector    
            else
                x = obj.options.(name);                    
            end               
        end
         
        function addEssentialPrivate(obj, varargin)
            %ADDESSENTIALPRIVATE label a boundary as essential
            %
            % Syntax
            %   addEssentialPrivate('PropertyName', PropertyValue,...);
            %
            % Description
            %   add_essential_boundary('PropertyName', PropertyValue,...)
            %   adds based on the properties an essential boundary
            %   condition that will be applied to the solution
            %
            % ADDESSENTIALPRIVATE Property Descriptions
            %   see ADDESSENTIAL

            % Gather the input
            opt.boundary = [];
            opt.value = [];
            opt.component = [];
            opt.clear = false;
            opt = gatherUserOptions(opt, varargin{:});

            if opt.clear
                obj.essential = struct('id',{},'value',{},'component',{});
            end

            % Test that id and value are given
            if isempty(opt.boundary) || isempty(opt.value);
                error('Solver:addEssentialPrivate','Both the id and value properties must be set.');
            end

            % Append storage data structure
            idx = length(obj.essential);
            obj.essential(idx+1).boundary = opt.boundary;
            obj.essential(idx+1).component = opt.component;   

            % Apply the value
            if ischar(opt.value) && ~isempty(obj.system);
                obj.essential(idx+1).value = obj.system.get(opt.value);

            elseif isnumeric(opt.value) || isa(opt.value,'function_handle');
                obj.essential(idx+1).value = opt.value;
                
           else
                error('Solver:addEssentialPrivate', 'Error extacting value from the System');
            end

        end

        function [u,ess] = applyConstraints(obj, varargin)
            %APPLYCONSTRAINTS Cretes solution and applies essential boundaries
            %
            % Syntax
            %   [u,ess] = applyConstraints(obj)
            %   [u,ess] = applyConstraints(obj,u)
            %
            % Description
            %   [u,ess] = applyConstraints(obj) return the solution (u) with
            %   essential boundaries applied and logcal vector of the
            %   essential boundaries (ess) as extracted from mesh.get_dof
            %
            %   [u,ess] = applyConstraints(obj,u) sames as above but uses an
            %   existing solution vector

            % Initlize the solution
            if nargin == 1;
                u = zeros(obj.mesh.n_dof,1);
            else
                u = varargin{1};
            end

            % Apply the essential boundary condions
            dof = false(obj.mesh.n_dof, length(obj.essential));
            for i = 1:length(obj.essential);
               dof(:,i) = logical(obj.mesh.getDof('Boundary', obj.essential(i).boundary, 'Component', obj.essential(i).component));
               
               if isa(obj.essential(i).value, 'function_handle');
                   x = obj.mesh.get_nodes(); 
                   u(dof(:,i)) = obj.essential(i).value(x(dof(:,i),:));
               else
                   u(dof(:,i)) = obj.essential(i).value;
               end
            end

            % Build the essential dofs
            ess = any(dof,2);
       end
   end
end

