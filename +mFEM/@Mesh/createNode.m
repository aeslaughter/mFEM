function createNode(obj,x)
    id = size(obj.node_map,1) + 1;
    obj.node_map(id,:) = x;
end