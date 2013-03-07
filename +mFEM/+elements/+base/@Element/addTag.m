function addTag(obj,tag,type)

% Update the elements
    for i = 1:length(obj);
%         % Go to next iteration if element is not on boundary
%        if strcmpi(type,'boundary');
%             if ~obj(i).on_boundary; continue; end
%        end

       % Loop through sides, mark side if dofs match
       found = false;
       for s = 1:obj(i).n_sides;
           if strcmpi(type,'boundary');
                if ~obj(i).sides(s).on_boundary; continue; end
           end
                tf = all(obj(i).nodes(obj(i).side_ids(s,:)).hasTag(tag));
%                 tf = all([obj(i).nodes(obj(i).side_ids(s,:)).on_boundary]);
                if tf;
                    obj(i).sides(s).tag{end+1} = tag;
                    found = true;
                end
       end

       if found;
           obj(i).tag{end+1} = tag;
           found = false;
       end
    end 
end