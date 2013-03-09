function T = test_Node(varargin)
    %TEST_NODE Tests all methods of the node class.
    %
    %     4-----5-----6
    %     |  1  |  2  |
    %     1-----2-----3

    % Create the test class
    T = mFEM.Test('name','Node',varargin{:});   
    
    % Create the node objects (use empty construtor method)
    try
        n = 6;   
        coord = rand(n,2);
        node(n) = mFEM.elements.base.Node();
    catch err
        T.caught(err);
        return;
    end

    % Intialize the nodes
    T.handle = node; % use test class to perform methods
    T.eval('init',{1:n,coord,3}); % initialize
    T.compare(node(4).getCoord(),coord(4,:)','Extracted coordinates');
    
    % Test dof set/get
    T.eval('setDof',{1}); % set 3 dofs per node
    T.handle = node(6);
    dof = T.eval('getDof',{},'Nout',1);
    T.compare(dof,[16,17,18],'Extracted completed dofs from node 6');
    dof = T.eval('getDof',{'Component',[2,3]},'Nout',1);
    T.compare(dof,[17,18],'Extracted 2,3 dof component from node 6');
    dof = T.eval('getDof',{'Component','x'},'Nout',1);
    T.compare(dof,16,'Extracted x-component dof from node 6');
    
    % Set the boundary id for nodes 3 & 6
    T.handle = node([3,6]);
    T.eval('setBoundaryFlag',{});
    T.compare(node(3).on_boundary,true,'Node 3 on boundary');
    T.compare(node(6).on_boundary,true,'Node 6 on boundary');
    
    % Add/get parent to nodes 5 and 2
    T.handle = node([2,5]);
    T.eval('addParent',{2});
    p = T.eval('getParents',{},'Nout',1);
    T.compare(p,2,'Parents extracted from nodes');

    % Test the tag
    T.eval('addTag',{'A'});
    T.compare(node(5).tag{1},'A','Tag added correctly');
    
    % Complete the test
    delete(T);



