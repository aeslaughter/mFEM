function createElement(obj,type,nodes)
    id = size(obj.elem_map,1) + 1;
    obj.elem_map(id,:) = nodes;
    obj.elem_type{id} = type;
end