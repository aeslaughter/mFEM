function test_Element
%
%
%     7-----8-----9
%     |  3  |  4  |  
%     4-----5-----6
%     |  1  |  2  |  
%     1-----2-----3


node_map = [0,0; 1,0; 2,0; 0,1; 1,1; 2,1; 0,2; 1,2; 2,2];
node(9) = mFEM.elements.base.Node();
node.init(1:9, node_map);

elems(4) = mFEM.elements.Quad4();
elems.init(1:4,node([1,2,5,4; 2,3,6,5; 4,5,8,7; 5,6,9,8]));

elems.findNeighbors();

T = mFEM.Test();
% T.compare(elems(1).nodes(3),elems(2).nodes(4), 'Elements sharing nodes');
% T.compare(elems(1).nodes(3),node(5), 'Elements nodes match input nodes');
% T.compare(elems(1).nodes(4).coord, [0,1,0], 'Element coordinates correct');
% T.compare(node(2).coord, [1,0,0], 'Input nodes coordinates correct');
% T.compare(elems(1).nodes(3).parents, elems, 'Parent test');
% T.compare(elems(1).sides(2).neighbor.id,2,'Find neighbor element 1, side 2'); 
% T.compare(elems(1).sides(3).neighbor.id,3,'Find neighbor element 1, side 3'); 
% T.compare(elems(2).sides(3).neighbor.id,4,'Find neighbor element 2, side 3'); 
% T.compare(elems(2).sides(4).neighbor.id,1,'Find neighbor element 2, side 4'); 
% T.compare(elems(3).sides(1).neighbor.id,1,'Find neighbor element 3, side 1'); 
% T.compare(elems(3).sides(2).neighbor.id,4,'Find neighbor element 3, side 2');
% T.compare(elems(4).sides(1).neighbor.id,2,'Find neighbor element 4, side 1'); 
% T.compare(elems(4).sides(4).neighbor.id,3,'Find neighbor element 4, side 4');


delete(elems);
delete(node);
clear nodes elems;
close all;