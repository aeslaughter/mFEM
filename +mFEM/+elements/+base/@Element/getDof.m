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
    % GETDOF Property Descriptions
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
    %           dof = getDof('Side',1,'local',true) or
    %           dof = getDof('Side',1,'-local')
    %
    %   Component
    %       scalar | 'x' | 'y' | 'z'
    %       Returns the dof associated with the vector space or
    %       multipe dof per node elements like the Beam element.

    % Set default options and gather the options
    options.side = [];
    options.local = false;
    options.component = [];
    opt = gatherUserOptions(options, varargin{:});
    
    glbl = [obj.nodes.dof]
    local = 1:length(glbl)
    
    
    
    
    
    
    
    