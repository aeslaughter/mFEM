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

    obj.node_map = feval(['mFEM.elements.',type,'.buildNodeMap'],varargin{:});

    obj.elem_map = feval(['mFEM.elements.',type,'.buildElementMap'], obj.node_map, varargin{:});

    obj.nodes = obj.buildNodes(obj.node_map);

    obj.elements = obj.buildElements(type, obj.elem_map, obj.node_map, obj.nodes);
    
    % Complete message
    if obj.options.time;
        tMessage(ticID);
    end
    
    obj.init();