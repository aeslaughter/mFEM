function init(obj,id,nodes)

    % The no. of elements
    n = length(obj);

    % Set/find the known constants
    lab = labindex; 
    [n_nodes,n_dim] = size(nodes);
    n_sides = size(obj(1).side_ids,1);
%     s = struct('neighbor',[],...
%                          'neighbor_side',[],...
%                          'on_boundary',...
%                          num2cell(true(n_sides,1)),...
%                          'tag',[]);

    n_nodes = num2cell(repmat(n_nodes,n,1));
    n_sides = num2cell(repmat(n_sides,n,1));
    n_dim = num2cell(repmat(n_dim,n,1));
    lab = num2cell(repmat(lab,n,1));
%     s = num2cell(repmat(s,n,1));
    id = num2cell(id);
%     no = mat2cell(nodes,ones(n,1),n_dim)

    [obj.n_nodes] = n_nodes{:};
    [obj.n_sides] = n_sides{:};
    [obj.n_dim] = n_dim{:};
    [obj.id] = id{:};
%     [obj.sides] = s{:};
    [obj.lab] = lab{:};
%     [obj.nodes] = no{:};


    for i = 1:length(obj);
%            obj(i).id = id(i);
    obj(i).nodes = nodes(i,:);
%            obj(i).lab = lab;
    obj(i).nodes.addParent(obj(i));
    obj(i).sides = struct('neighbor',[],...
                         'neighbor_side',[],...
                         'on_boundary',...
                         num2cell(true(obj(i).n_sides,1)),...
                         'tag',[]);
    end
end