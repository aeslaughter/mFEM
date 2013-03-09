function addBoundary(obj, id, varargin)
    %ADDBOUNDARY Labels elements and sides with numeric tags.
    %
    % Syntax
    %   addBoundary(id, Limit1,...)
    %   addBoundary(id)
    %
    % Description
    %   addBoundary(id) labels all unidentified elements and sides
    %       with the specified id, which must be an interger
    %       value greater than zero.
    %
    %   addBoundary(id, Limit1,...) allows 
    %       for customization to what boundaries are tagged, see 
    %       the descriptions below. It is possible to
    %       supply multiple limits For example,
    %           add_boundary(1,'left','right','x==1')
    %       will loop through each value and apply the boundary
    %       id to each.
    %
    % Descriptions of Limits
    %   Locations:
    %       'left' | 'right' | 'top' | 'bottom' | 'front' | 'back'
    %       Adds a boundary id to a side identified by a string
    %       that describes its location. Note, in 1D use 'left' or
    %       'right', in 2D 'top' and 'bottom' are added, and in 3D
    %       all options are available.
    %
    %   Functions:
    %       string | cell array of strings
    %       It is possible to tag boundaries using test functions.
    %       For example, 'x==1' will tag boundaries that have a
    %       node location at 1. Mutliple conditions should be
    %       expressed in a cell array. For example, {'x==1','y<2'}
    %       will id elements with nodes at x = 1 AND y < 2.
    %
    % Examples
    %   The following labels the right hand side of the mesh and
    %   the node point (1,1) with id 1.
    %       add_boundary(1,'right',{'x==1','y==1'});
    %
    %   The following labels the right hand side of the mesh and
    %   the elements with nodes on x == 1 OR y == 1 to 1.
    %       add_boundary(1,'right','x==1','y==1'});
    % 
    % See Also
    %   ADDSUBDOMAIN

    % Check if system is initialized
    if obj.initialized;
        warning('Mesh:addBoundary:MeshInitialized',...
            'The Mesh was previously initailized, for these tags to be applied the init() method must be called again.');
    end

    % Special case, tag all untagged
    if nargin == 2;
        obj.idEmptyBoundary(id);
        
    % Loop through the supplied options    
    else
        obj.addTag(id, 'boundary', varargin{:});
    end
end