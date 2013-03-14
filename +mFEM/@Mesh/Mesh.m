classdef Mesh < handle
    %MESH Summary of this class goes here
    %   Detailed explanation goes here
    %
    % mention dimension independance, is possible to mix elements of
    % different dimensions and n dofs (not tested)
    
    properties (GetAccess = public, SetAccess = protected)
        elements = Composite();
        nodes = Composite(); 
        n_nodes = uint32([]);                
        n_elements = uint32([]);            
        n_dof = uint32([]);            
        options = struct('time', false, 'space', 'scalar');
    end
    
    properties (Hidden, Access = protected)
        node_map = [];
        elem_map = uint32([]);
        elem_type = {};
        initialized = false;
        node_map_codist;
        elem_map_codist;
        node_tag_map = Composite();
        elem_tag_map = Composite();
        tag = struct('name',{},'type',{});
    end
    
    methods
        createNode(obj,x);
        createElement(obj,type,nodes);       
        init(obj);
        grid(obj,varargin);
        addBoundary(obj,id,varargin);
        addSubdomain(obj,id,varargin);
        update(obj);
        el = getElements(obj,varargin);
        no = getNodes(obj,varargin);
        dof = getDof(obj,varargin);
        tag = getTag(obj,name);
        plot(obj,varargin);

        function obj = Mesh(varargin)
            obj.options = gatherUserOptions(obj.options, varargin{:});  
        end
    end
    
    methods (Hidden, Access = protected)
        addTag(obj, id, type, varargin)
        out = gatherComposite(obj,varargin);
        nodes = buildNodes(obj);
        [elements, nodes] = buildElements(obj,nodes);
    end
    
    methods (Hidden, Static, Access = protected)
        no = getOffLabNodes(e_map, node_map, nodes)
    end
end