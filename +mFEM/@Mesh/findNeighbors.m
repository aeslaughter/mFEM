function findNeighbors(obj)
    % Locates elements that share a side
    %
    % Syntax
    %   findNeighbors()   
    %
    % Description
    %   findNeighbors() locates neighboring elements and
    %   neighboring sides for the finite element mesh, it sets the
    %   various neighbor and side related properties for the
    %   elements in the mesh.

    % Display wait message
    if obj.options.time;
        ticID = tMessage('Locating neighbor elements...');
    end

    elements = obj.elements;
    spmd
        elements.findNeighbors();
    end

    % Complete message
    if obj.options.time;
        tMessage(ticID);
    end;
end 