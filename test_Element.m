function test_Element
    %
    %     12---13----14---15
    %     |  7  |  8  |  9 |
    %     8-----9----10---11
    %     |  4  |  5  |  6 | 
    %     4-----5-----6----7
    %     |  1  |  2  |  3 | 
    %     1-----2-----3----4


    mesh = mFEM.Mesh();
    mesh.grid('Quad4',0,3,0,3,3,3);
    elem = mesh.getElements();
    node = mesh.getNodes();
    
    T = mFEM.Test();
    T.compare(elem(1).nodes(3),elem(2).nodes(4), 'Elements sharing nodes');
    T.compare(elem(1).nodes(3),node(5), 'Elements nodes match input nodes');
    T.compare(elem(1).nodes(4).coord, [0;1;0], 'Element coordinates correct');
    T.compare(node(2).coord, [1;0;0], 'Input nodes coordinates correct');
    T.compare(elem(1).nodes(3).parents, elem, 'Parent test');

    % Test the neighbors were located
%     testNeighbors(T,elems);
%     testBoundarySides(T,elems);
%     testNodeBoundary(T,node);

    delete(elem);
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