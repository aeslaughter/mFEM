function work
% 
% node_map = mFEM.elements.Hex8.buildNodeMap(0,1,0,1,0,1,2,2,2);
% elem_map = mFEM.elements.Hex8.buildElementMap(0,1,0,1,0,1,2,2,2);
mesh = mFEM.Mesh();
mesh.grid('Quad4',0,1,0,1,4,4);

elem = mesh.getElements('-parallel');

spmd
    x = [elem.nodes.getCoord()]
    [elem.side_ids]
    
    
    
    
end


