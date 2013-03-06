function findNeighbors(elements,time_flag)
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
    if time_flag;
        ticID = tMessage('Locating neighbor elements...');
    end
    
    % Use the Element class function
    spmd
        elements.findNeighbors();
    end

    % Complete message
    if time_flag;
        tMessage(ticID);
    end;
end 