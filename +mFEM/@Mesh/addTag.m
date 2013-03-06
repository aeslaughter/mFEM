function addTag(obj, tag, func)
    %ADD_TAG Adds boundary and subdomain tags based on function or string
    %
    % Syntax
    %    add_tag(tag, FuncString)
    %    add_tag(tag, FuncCell)
    %
    % Description
    %    addTag(tag, FuncString) adds a tag based on the string expression
    %    in FuncString.
    %
    %    addTag(tag, FuncCell) same as above, but adds a tag
    %    based on all of the string expressions in FuncCell.
    %
    % See Also ADDBOUNDARY ADDSUBDOMAIN
    
    node_map = obj.node_map;
    nodes = obj.nodes;
    
    % Convert numeric tags to strings
    if isnumeric(tag); tag = num2str(tag); end
    
    % Parse the input
    [fcn,col,value] = parseInput(func, node_map);
    
    % Evaluate the functions on each lab

    elements = obj.elements;
    spmd
        local = getLocalPart(node_map);
        idx = zeros(length(local),length(fcn));
        
        for i = 1:length(fcn);
            idx(:,i) = feval(fcn{i},local(:,col(i)),value(i));
        end
        idx = all(idx,2); 
        nodes(idx).addTag(tag);
        
        % Update the elements
        for e = 1:length(elements);
            
           % The current element
           elem = elements(e);
           
           % Go to next iteration if element is not on boundary
           if ~elem.on_boundary; continue; end
           
           % Loop through sides, mark side if dofs match
           found = false;
           for s = 1:elem.n_sides;
               if ~elem.sides(s).on_boundary; continue; end
               
                    tf = all(elem.nodes(elem.side_ids(s,:)).hasTag(tag));
                    if tf;
                        elem.sides(s).tag{end+1} = tag;
                        found = true;
                    end
           end

           if found;
               elem.tag{end+1} = tag;
               found = false;
           end
           
        end 
    end
end

function [fcn,col,value] = parseInput(str,node_map)

    C = {'left','right','top','bottom','front','back'};
    
    % Convert character input into a cell
    if ischar(str)
       str = {str};
    end

    % Loop throug each function specified and apply the id
    value = zeros(length(str),1);
    col = value;
    fcn = cell(size(col));
    for i = 1:length(str);

        % Extract the current function
        s = strtrim(str{i});
        
        % Account for boundary flag input
        if any(strcmpi(s,C));
            [fcn{i},col(i),value(i)] = parseFlagInput(s,node_map);
        else
            [fcn{i},col(i),value(i)] = parseFuncInput(s);
        end
    end
end

function [fcn,col,value] = parseFlagInput(str,node_map)

    switch lower(str);
        case 'left';   fcn = 'le'; col = 1; value = min(node_map(:,1)) + eps;
        case 'right';  fcn = 'ge'; col = 1; value = max(node_map(:,1)) - eps;
        case 'bottom'; fcn = 'le'; col = 2; value = min(node_map(:,2)) + eps;
        case 'top';    fcn = 'ge'; col = 2; value = max(node_map(:,2)) - eps;
        case 'front';  fcn = 'le'; col = 3; value = min(node_map(:,3)) + eps;
        case 'back';   fcn = 'ge'; col = 3; value = max(node_map(:,3)) - eps;
    end
    value = gather(value);
end

function [fcn,col,value] =  parseFuncInput(str)

    A = {'<','>','==','<=','>=','~='};
    B = {'lt','gt','eq','le','ge','ne'};
    
    % Determine the column to operate on, also test that the
    % mesh is the proper dimension for the desired test
    if strcmpi(str(1),'x'); 
       col = 1;
    elseif strcmpi(str(1),'y');
       col = 2;
    elseif strcmpi(str(1),'z');
       col = 3;
    else
       error('FEmesh:addTag:parseFuncInput:InvalidInput',...
            'The input string ''%s'' is invalid.', s);
    end

    idx = regexp(str,'[0-9]'); 
    value = str2double(str(idx));
    idx = strcmp(str(2:idx-1),A);
    fcn = B{idx};
end
    