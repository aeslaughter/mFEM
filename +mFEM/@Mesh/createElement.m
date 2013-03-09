function elem = createElement(obj,type,nodes)

    if isnumeric(nodes);
        nodes = obj.nodes(nodes);
    end

    id = length(obj.elements) + 1;
    elem = feval(['mFEM.elements.',type],id,nodes);
    obj.elements{id} = elem;
end