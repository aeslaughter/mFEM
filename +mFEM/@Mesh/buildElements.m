function [elements,nodes] = buildElements(obj,nodes) 
    %BUILDELEMENTS (protected) Create elements from element map and nodes
    %
    % Syntax
    %   [elements,nodes] = buildElements(nodes)
    %
    % Description
    %
    %
    
    % Get variables for use
    elem_map = obj.elem_map;
    node_map = obj.node_map;
    
    % Get the element type (todo) this will allow for mixed meshes
    elem_type = obj.elem_type;
    type = gather(elem_type(1));type = type{1};
    
    % Check that enough elements are present for parallel
    if matlabpool('size') > obj.n_elements;
        error('Mesh:buildElements:ParallelFail','The number of elements (%d) must exceed the number of processors (%d)',obj.n_elements,matlabpool('size'));
    end
    
    % Serial case
    if matlabpool('size') == 0;
        spmd
            n_elem = size(elem_map,1); % no. of elements
            elements(n_elem) = feval(['mFEM.elements.',type]);
            elements.init(1:n_elem, nodes(gather(elem_map)),'-addParents');
        end
        return % done
    end

    % Parallel case
    spmd
        % Extract local part of the element map 
        e_id  = globalIndices(elem_map, 1); % global element no.

        % Extract all the nodes that are needed on this processor
        no = obj.getOffLabNodes(elem_map, node_map, nodes); 

        % Create the elements
        n = size(no,1);
%         elements(n,1) = feval(['mFEM.elements.',type]);
%         elements.init(e_id,no,'-addParents');
    end
    elements = {};
end