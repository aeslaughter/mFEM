function grid(obj, type, varargin)

    % Test for valid type
    d = dir(fullfile(getpref('MFEM_PREF','ROOT_DIR'),'+mFEM','+elements','*.m'));
    c = struct2cell(d);
    if ~any(strcmp([type,'.m'],c(1,:)));
        msg = 'The %s element was not found in the +mFEM/+elements directory.';
        error('Mesh:grid:InvalidType', msg, type);
    end

    % Display wait message
    if obj.options.time;
        ticID = tMessage('Generating Grid...');
    end

    % Create the node map using the element
    node_map = feval(['mFEM.elements.',type,'.buildNodeMap'],varargin{:});
    [obj.node_map, obj.node_map_codist] = createCodistributed(node_map);
    
    % Create the element map using the element
    elem_map = feval(['mFEM.elements.',type,'.buildElementMap'], obj.node_map, varargin{:});
    [obj.elem_map, obj.elem_map_codist] = createCodistributed(elem_map);

    % Create a map of the type
    obj.elem_type = distributed(repmat({type},size(elem_map,1),1));
    
    % Complete message
    if obj.options.time;
        tMessage(ticID);
    end

    % Setup the Mesh object
    obj.setup();
end

function [in,codist] = createCodistributed(in)
    spmd 
        if ~iscodistributed(in);
            codist = codistributor1d(1,codistributor1d.unsetPartition,size(in));
            in = codistributed(in, codist);
        else
            codist = getCodistributor(in);
        end
    end
end