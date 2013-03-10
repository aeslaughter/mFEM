function init(obj,id,nodes)
   lab = labindex; 
   for i = 1:length(obj);
       obj(i).id = id(i);
       obj(i).nodes = nodes(i,:);
       obj(i).lab = lab;
       obj(i).nodes.addParent(obj(i));
       obj(i).sides = struct('neighbor',[],...
                             'neighbor_side',[],...
                             'on_boundary',...
                             num2cell(true(obj(i).n_sides,1)),...
                             'tag',[]);
   end
end