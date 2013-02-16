function test_Vector
    %TEST_VECTOR Tests the Vector class (serial and parallel)
    
% Start a parallel session, if it does not exist
if matlabpool('size') == 0;
    matlabpool 2;
end

% Create a instance of the mFEM.Test class
T = mFEM.Test();

% Create a constant vector
v = mFEM.Vector(11,53);

% Check that out-of-range error is working
try
    v.get(12)
catch err
   T.compare('Vector:get:OutOfRange', err.identifier, 'Out of range error');
end

% Check extraction of specific values of the vector
T.compare(53,v.get(6),'Parallel: Constant initilized, get scalar');
T.compare([53;53], v.get(6:7), 'Parallel: Constant initilized, get vector (1 lab)');
T.compare([53;53;53], v.get([1,7,8]), 'Parallel: Constant initilized, get vector (2 lab)');

% Create a Codistributed array
X = distributed(1:11);
spmd
    local = getLocalPart(X);
    codist = getCodistributor(X); 
    vec = codistributed.build(2*local, codist);
end

% Create a vector from the existing codistributed array
v = mFEM.Vector(vec);

% Check extraction of specific values of the vector
T.compare(6,v.get(3), 'Parallel: Codistributed initialized, get scalar');
T.compare([6;8],v.get(3:4), 'Parallel: Codistributed initialized, get vector (1 lab)');
T.compare([4;12;22],v.get([2,6,11]), 'Parallel: Codistributed initialized, get vector (2 lab)');

% Check that the add function works
v.add(10,1);
T.compare(12,v.get(1), 'Parallel: add scalar');

v.add(100,[2,6,11]);
T.compare(100+[4;12;22],v.get([2,6,11]), 'Parallel: add vector (2 lab)');

