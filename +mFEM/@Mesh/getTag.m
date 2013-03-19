function tag = getTag(obj,name)

    if ~iscell(name);
        name = {name};
    end
    
    n = length(name);
    tag = struct('name',{},'type',{});
    k = 1;
    for i = 1:n

        if isnumeric(name{i}); 
            name{i} = num2str(name{i});
        end

        ix = strcmp(name{i},{obj.tag.name});

        if ~isempty(ix);
            tag(k) = obj.tag(ix);
            k = k + 1;
        end
    end
    
end