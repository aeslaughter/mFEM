function addTag(obj, tag, type, varargin)
    %ADDTAG Adds boundary and subdomain tags based on function or string
    %
    % Syntax
    %    addTag(tag, FuncString, type)
    %    addTag(tag, FuncCell, type)
    %
    % Description
    %    addTag(tag, FuncString, type) adds a tag based on the string 
    %    expression in FuncString. The type should be either 'boundary' or
    %    'subdomain', the former of restricts the application of the tag to
    %    elements and nodes on boundaries.
    %
    %    addTag(tag, FuncCell) same as above, but adds a tag based on all 
    %    of the string expressions in FuncCell.
    % 
    % See Also ADDBOUNDARY ADDSUBDOMAIN

    node_map = obj.node_map;
    nodes = obj.nodes;
    
    % Convert numeric tags to strings
    if isnumeric(tag); tag = num2str(tag); end
    
    % Parse the input
    [fcn,col,value,use_all] = parseInput(node_map, varargin{:});
    
    % Evaluate the functions on each lab

    elements = obj.elements;
    spmd
        local = getLocalPart(node_map);
        idx = zeros(length(local),length(fcn));
        
        for i = 1:length(fcn);
            idx(:,i) = feval(fcn{i},local(:,col(i)),value(i));
        end
        if use_all;
            idx = all(idx,2); 
        else
            idx = any(idx,2); 
        end
        nodes(idx).addTag(tag);
        
        elements.addTag(tag,type);
    end
end

function [fcn,col,value,use_all] = parseInput(node_map, varargin)

    flags = {'left','right','top','bottom','front','back'};

    
    if nargin == 2;
        str = varargin{1};
        use_all = true;
    else
        str = varargin;
        use_all = false;
    end


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
        if any(strcmpi(s,flags));
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
            'The input string ''%s'' is invalid.', str);
    end

    idx = regexp(str,'[0-9]'); 
    value = str2double(str(idx));
    idx = strcmp(str(2:idx-1),A);
    fcn = B{idx};
end
    