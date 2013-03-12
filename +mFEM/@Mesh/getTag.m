function tag = getTag(obj,name)

    if isnumeric(name); 
        name = num2str(name);
    end
    
    ix = strcmp(name,{obj.tag.name});
    
    if ~isempty(ix);
        tag = obj.tag(ix);
    else
        tag = [];
    end
    
end