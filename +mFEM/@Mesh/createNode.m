function node = createNode(obj,x,varargin)
    
    opt.space = obj.options.space;
    opt = gatherUserOptions(opt,varargin{:});

    id = length(obj.nodes) + 1;
    node = mFEM.elements.base.Node(id,x,opt.space);     
    obj.nodes{id} = node; 
end