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
    opt.side = [];
    opt.local = false;
    opt.component = [];
    opt = gatherUserOptions(opt, varargin{:});
    
    % Local only operates on single element
    if opt.local && length(obj) > 1;
        error('Element:getDof:InvalidUseOfLocal','The local property only functions on a single element');
    elseif opt.local
        glbl = obj.dof; % all dofs, global
        x = obj.nodes.getDof('Component',opt.component); % desired dofs
        ix = ismember(glbl,x);     % location in glbl of desired dofs
        lcl = 1:length(obj.dof);   % local dofs
        dof = lcl(ix);             % desired local dofs
    end
    
    
    
    
    
    
    
    