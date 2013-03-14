function createNode(obj,x)
    for i = 1:size(x,1);
        id = size(obj.node_map,1) + 1;
        obj.node_map(id,:) = x(i,:);
    end
end