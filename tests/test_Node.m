function test_Node
    %TEST_NODE Tests the mFEM.elements.base.Node class

n = mFEM.elements.base.Node(1,[1,2,3]);
T = mFEM.Test();

[x,y,z] = n.get();
T.compare([x,y,z],[1,2,3],'individual x,y,z output');

coord = n.get();
T.compare(coord,[1,2,3],'numeric array output');

x = n.get(1);
T.compare(x,1,'numeric x output');

x = n.get([1,2]);
T.compare(x,[1,2],'numeric x,y output');

y = n.get('y');
T.compare(y,2,'character y output');



