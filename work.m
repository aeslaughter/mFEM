function work

n = 100000;

tic;
node = mFEM.elements.base.Node.empty(n,0);
for i = 1:n;
    node(i) = mFEM.elements.base.Node(i,rand(2,1));
    node(i).addParent(i);
end
c = num2cell(node);
toc;
    
tic;
node(n) = mFEM.elements.base.Node();
for i = 1:length(node)
    node(i).init(i,rand(2,1));
end
node(1:n).addParent(1:n)
c = num2cell(node);
toc;

tic;
node(n) = mFEM.elements.base.Node();
node(1:n).init(1:n, rand(n,2));
c = num2cell(node);
toc;

tic;
node = mFEM.elements.base.Node(1:n, rand(n,2));
c = num2cell(node);
toc;