function createElement(obj,type,nodes)
    for i = 1:size(nodes,1);
        id = size(obj.elem_map,1) + 1;
        obj.elem_map(id,:) = nodes(i,:);
        obj.elem_type{id} = type;
    end
end