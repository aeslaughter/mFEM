function out = gatherComposite(obj,varargin)
    %GATHERCOMPOSITE (protected) Returns selection of nodes or elements
    %
    % Syntax
    %   no = getComposite()
    %   no = getComposite(id)
    %   no = getComposite(...,'PropertyName',PropertyValue,...)
    %
    % Description
    %   no = getComposite() returns the node or element objects for the 
    %   entire mesh.
    %
    %   no = getComposite(id) returns the objects corresponding to the
    %   global ids supplied in id.
    %
    %   no = getComposite(...,'PropertyName',PropertyValue,...) returns the 
    %   objects as above but limits the selection further based on the
    %   property parings supplied, which are detailed below.
    %
    % GETNODES Property Descriptions
    %   Name (Required)
    %       'elements' | 'nodes'
    %       Triggers which type of data to return.
    %
    %   Tag
    %       scalar | char
    %       Limits the objects returned to only those with the
    %       supplied tag, the term tag is applied to both those added with
    %       the addBoundary and addSubdomain commands.
    %
    %   Lab
    %       scalar | vector
    %       Limits the objects returned to the given processors in
    %       parallel applications
    %
    %   ID
    %       scalar | vector
    %       Same as supplied id discussed in description above, this is
    %       simply an alternative method for supplying the ids. The values
    %       given in this option will overwrite those supplied directly.
    %

    % Gather the input and extract the desired objects
    [comp,id,codist,all_ids,opt] = parseGatherCompositeInput(varargin{:});

    % Extract in serial case    
    if matlabpool('size') == 0;
        elem = obj.elements{1};
        if ~all_ids;
            elem = elem(id);
        end
    else
            
        if ~all_ids;
            spmd
                % IDs of the elements stored on this lab
                local_id = globalIndices(map,1);

                % local ids for the desired elements
                idx = id >= min(local_id) & id <= max(local_id);

                % Build a Composite
                comp_idx = comp(idx);
            end
        else
            comp_idx = comp;
        end
    
        if isempty(lab)
            out = comp_idx{1};
            for i = 2:length(comp_idx);
                out = [out; comp_idx{i}];
            end   
        else
            out = comp_idx{lab};
        end
    end
end

function [comp,id,codist,all_ids,opt] = parseGatherCompositeInput(varargin)
    %PARSEGATHERCOMPOSITEINPUT Prepares input for use by main function
    
    % No input given by user, get all nodes
    if nargin == 0;
        opt.id = 1:obj.n_nodes;
        properties = {};
    
    % Ids without properties given by user    
    elseif nargin == 1 && isnumeric(varargin{1});
        opt.id = varargin{1};
        properties = {};
        
    % Ids with properties given by the user    
    elseif nargin >= 2 && isnumeric(varargin{1});
        opt.id = varargin{1};
        properties = varargin(2:end);
        
    % Only properties supplied    
    else
        opt.id = 1:obj.n_nodes;
        properties = varargin;
    end
    
    % Gather the properties
    opt.tag = {};
    opt.lab = [];
    opt.name;
    opt = gatherUserOptions(opt,properties{:});
    
    % Input the ids property to the output
    id = opt.id;
    
    % Extract the data to work with
    switch lower(opt.name)
        case {'elem','element','elements'};
            codist = obj.elem_map_codist;
            comp = obj.elements;
            n = sise(obj.elem_map,1);
        case {'nodes','node'};
            codist = obj.elem_map_codist;
            comp = obj.elements;
            n = sise(obj.elem_map,1);
        otherwise
            error('Mesh:gatherComposite:InvalidName','The name input must be either ''elements'' or ''nodes''.');
    end
    
    % Determine if all of the elements are requested
    all_ids = false;
    if islogical(id) && sum(id) == n;
        all_ids = true;
    elseif length(id) == n;
        all_ids = true;
    end
end
