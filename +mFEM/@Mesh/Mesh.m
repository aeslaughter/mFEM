classdef Mesh < handle
    %MESH Summary of this class goes here
    %   Detailed explanation goes here
    %
    % mention dimension independance, is possible to mix elements of
    % different dimensions (not tested)
    
    properties (GetAccess = public, SetAccess = protected)
        elements = Composite();
        nodes = Composite(); 
        n_nodes = uint32([]);                
        n_elements = uint32([]);            
        n_dof = uint32([]);            
        options = struct('time', false, 'space', 'scalar');
    end
    
    properties (Hidden, Access = protected)
        node_map = Composite();
        elem_map = Composite();
        elem_type = Composite();
        initialized = false;
        node_map_codist;
        elem_map_codist;
        node_tag_map = Composite();
        elem_tag_map = Composite();
        tag = {};
    end
    
    methods
        el = createElement(obj,type,nodes);
        grid(obj,varargin);
        addBoundary(obj,id,varargin);
        addSubdomain(obj,id,varargin);
        el = getElements(obj,varargin);
        no = getNodes(obj,varargin);
        dof = getDof(obj,varargin);
        plot(obj,varargin);

        function obj = Mesh(varargin)
            obj.options = gatherUserOptions(obj.options, varargin{:});  
        end
    end
    
    methods (Hidden, Access = protected)
        init(obj);
        addTag(obj, id, type, varargin)
        idEmptyBoundary(obj,id);
        out = gatherComposite(obj,name,id,tag,lab);
        nodes = buildNodes(obj);
        [elements, nodes] = buildElements(obj,nodes);
    end
end

