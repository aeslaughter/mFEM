function addTag(obj, tag, type, varargin)
    %ADDTAG (protected) Adds tags based on string based flags and functions
    %   Acts to add tags to nodes, elements, and sides. This is a private
    %   function, the user should refer to the addBoundary or addSubdomain
    %   functions for adding tags.
    %
    %   ADDTAG only interacts with the node_tag_map and the elem_tag_map.
    %   It also adds the tags to the nodes on each process, the tagging of
    %   the elements takes place when the update method is called.
    %
    % Syntax
    %    addTag(tag, type, FuncString,...)
    %    addTag(tag, type, FuncCell)
    %
    % Description
    %    addTag(tag, type, FuncString,...) adds tag based on a single or
    %    multiple string inputs. The tag is added if ANY of the criteria
    %    are met. The type should be either 'boundary' or 'subdomain', the 
    %    former of restricts the application of the tag to elements and 
    %    nodes on boundaries. The later makes no restrictions on the
    %    location.
    %
    %    addTag(tag, FuncCell) similar to above but each string is placed
    %    inside a single cell array. In this case the tag is added if ALL
    %    of the criteria are met.
    % 
    % See Also addBoundary addSubdomain
    %
    %----------------------------------------------------------------------
    %  mFEM: A Parallel, Object-Oriented MATLAB Finite Element Library
    %  Copyright (C) 2013 Andrew E Slaughter
    % 
    %  This program is free software: you can redistribute it and/or modify
    %  it under the terms of the GNU General Public License as published by
    %  the Free Software Foundation, either version 3 of the License, or
    %  (at your option) any later version.
    % 
    %  This program is distributed in the hope that it will be useful,
    %  but WITHOUT ANY WARRANTY; without even the implied warranty of
    %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %  GNU General Public License for more details.
    % 
    %  You should have received a copy of the GNU General Public License
    %  along with this program. If not, see <http://www.gnu.org/licenses/>.
    %
    %  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
    %----------------------------------------------------------------------
    
    % Extract the variables for use
    nodes = obj.nodes;                  
    node_map = obj.node_map;            
    node_tag_map = obj.node_tag_map;
    
    % Convert numeric tags to strings
    if isnumeric(tag); 
        tag = num2str(tag); 
    end
    
    % Account for single text input, make it a cell
    if length(varargin) == 1 && ischar(varargin{1});
        varargin{1} = varargin(1);
    end
    
    % Build tag structure, append to the complete list of tags
    tag = struct('name',tag,'type',type);
    obj.tag(end+1) = tag;
    cnt = length(obj.tag); 

    % Parse the input
    [fcn,col,value,use_all] = parseAddTagInput(node_map, varargin{:});
    
    % Evaluate the functions on each lab
    spmd
        % Get the local node_map and global indices 
        local = getLocalPart(node_map);
        n_idx = zeros(length(local),length(fcn));
        
        % Evaluate the node based functions
        for i = 1:length(fcn);
            n_idx(:,i) = feval(fcn{i},local(:,col(i)),value(i));
        end
        
        % Account for the two types of input AND (all) or OR (any)
        if use_all;
            n_idx = all(n_idx,2); 
        else
            n_idx = any(n_idx,2); 
        end

        % Limit to ids to those on the boundary
        if strcmpi(type, 'boundary');
           n_idx(:,2) =  [nodes.on_boundary];
           n_idx = all(n_idx,2);
        end
 
        % Apply the nodal tags
        nodes(n_idx).addTag(tag);
        
        % Update the node tag map
        gi = globalIndices(node_tag_map,1);
        node_tag_map(gi,cnt) = n_idx;

        % Update the node tag map
        elem = [nodes(n_idx).parents];      % elements affected
        if ~isempty(elem);              
            [~,ix] = unique([elem.id]);     % ids of elements
            elem = elem(ix);                % element classes
            if ~isempty(elem);
                % Limit elements altered to this lab
                e_idx = [elem.lab]==labindex;
%                 elem(e_idx).addTag(tag,type);

                % Update the element map
                gi = [elem(e_idx).id];
                elem_tag_map(gi,cnt) = e_idx;
            end
        else
            elem_tag_map = [];
        end
    end
    
    % Return codistributed maps  to the object
    obj.node_tag_map = node_tag_map;
    obj.elem_tag_map = elem_tag_map;
    obj.nodes = nodes;
end

function [fcn,col,value,use_all] = parseAddTagInput(node_map, varargin)
    %PARSEADDTAGINPUT Seperates input for the two possible inputs
    %   The addTag method accepts two types of inputs: position tags (e.g.,
    %   'right') and boolean tests (e.g., 'x>1'). This subfunction parses
    %   those two options and calls functions that builds the output
    %   variable needed by the addTag method.
    
    % Single cell as input case
    if nargin == 2; 
        str = varargin{1};
        use_all = true;
        
    % Case when a set of strings is given    
    else % string input
        str = varargin;
        use_all = false;
    end

    % Check that all input are strings
    if ~iscellstr(str);
        error('Mesh:addTag:InvalidInput','Input must be either a single cell or a list of strings.');
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
        if any(strcmpi(s,{'left','right','top','bottom','front','back'}));
            [fcn{i},col(i),value(i)] = parseAddTagFlagInput(s,node_map);
        else
            [fcn{i},col(i),value(i)] = parseAddTagFuncInput(s);
        end
    end
end

function [fcn,col,value] = parseAddTagFlagInput(str,node_map)
    %PARSEADDTAGFLAGINPUT Parses function, column, and value for flag input
    %   One possible type of input for addTag are keywords such as 'right'
    %   or 'left', this function gives the proper information for the
    %   addTag method to perform the necessary analysis.

    switch lower(str);
        case 'left';   fcn = 'le'; col = 1; value = min(node_map(:,1));% + eps;
        case 'right';  fcn = 'ge'; col = 1; value = max(node_map(:,1));% - eps;
        case 'bottom'; fcn = 'le'; col = 2; value = min(node_map(:,2));% + eps;
        case 'top';    fcn = 'ge'; col = 2; value = max(node_map(:,2));% - eps;
        case 'front';  fcn = 'le'; col = 3; value = min(node_map(:,3));% + eps;
        case 'back';   fcn = 'ge'; col = 3; value = max(node_map(:,3));% - eps;
    end
    value = gather(value);
end

function [fcn,col,value] =  parseAddTagFuncInput(str)
    %PARSEADDTAGFUNCINPUT Parses function, column, and value for functions
    %   One possible type of input for addTag are functions such as 'x<1'
    %   or 'y==x', this function gives the proper information for the
    %   addTag method to perform the necessary analysis.    
    
    % Cell strings for boolean operators accepted
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

    % Assign value and fcn output
    idx = regexp(str,'[><=~]'); 
    value = str2double(str(idx(end)+1:end));
    fcn = B{strcmp(str(idx),A)};
end
    