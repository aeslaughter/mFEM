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

    error('Not done yet');
    
    % Extract the nodes for the side
    dof = obj.side_dof(id,:);
    node = obj.nodes(dof,:);

    % Create the side element
    side = feval(['mFEM.elements.',obj.side_type], NaN, node, ...
        'Space', obj.opt.space);

    % Set the global dofs
    side.global_dof = obj.global_dof(dof);
end