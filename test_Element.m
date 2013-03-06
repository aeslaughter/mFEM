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
    T.compare(elems(1).nodes(3),elems(2).nodes(4), 'Elements sharing nodes');
    T.compare(elems(1).nodes(3),node(5), 'Elements nodes match input nodes');
    T.compare(elems(1).nodes(4).coord, [0;1;0], 'Element coordinates correct');
    T.compare(node(2).coord, [1;0;0], 'Input nodes coordinates correct');
    T.compare(elems(1).nodes(3).parents, elems, 'Parent test');

    % Test the neighbors were located
    testNeighbors(T,elems);
    testBoundarySides(T,elems);
    testNodeBoundary(T,node);

    delete(elems);
    delete(node);
    clear nodes elems classes;
end

function testNeighbors(T,elems)
    x = [1,2,2; 1,3,3; 2,3,4; 2,4,1; 3,1,1; 3,2,4; 4,1,2; 4,4,3];
    for i = 1:length(x);
        e = x(i,1); s = x(i,2); id = x(i,3);
        msg = sprintf('Finding neighbor of element %d, side %d',e,s);
        T.compare(elems(e).sides(s).neighbor.id,id,msg); 
    end
end

function testBoundarySides(T,elems)
    x = [1,4; 1,2; 3,4; 2,3];
    for i = 1:length(x);
        for s = 1:2;
            msg = sprintf('Setting boundary flag of element %d, side %d',i,x(i,s));
            T.compare(elems(i).sides(x(i,s)).on_boundary,true,msg); 
        end
    end
end

function testNodeBoundary(T,node)
    tf = true(9,1);
    tf(5) = false;
    
    for i = 1:9;
        msg = sprintf('Setting Node %d boundary flag',i);
        T.compare(node(i).on_boundary,tf(i),msg);
    end
end