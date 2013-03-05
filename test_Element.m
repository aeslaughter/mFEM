function test_Element


node_map = [0,0; 1,0; 2,0; 0,1; 1,1; 2,1];
node(6) = mFEM.elements.base.Node();
node.init(1:6, node_map);

elems(2) = mFEM.elements.Quad4();
elems.init(1:2,node([1,2,5,4; 2,3,6,5]));

T = mFEM.Test();
T.compare(elems(1).nodes(3),elems(2).nodes(4), 'Elements sharing nodes');
T.compare(elems(1).nodes(3),node(5), 'Elements nodes match input nodes');
T.compare(elems(1).nodes(4).coord, [0,1,0], 'Element coordinates correct');
T.compare(node(2).coord, [1,0,0], 'Input nodes coordinates correct');
T.compare(elems(1).nodes(3).parents, elems, 'Parent test');
