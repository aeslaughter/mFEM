function test_Hex8

% Call the Test class for this file
T = mFEM.Test(mfilename('fullfile'));

% Create a one-element mesh
mesh = mFEM.FEmesh();
mesh.grid(-1,1,-1,1,-1,1,1,1,1,'Element','Hex8');
mesh.init();

% System
sys = mFEM.System(mesh);
sys.addConstant('E',1,'v',0.3);
sys.addConstant('D','E/((1-2*v)*(1+v))*[1-v,v,v,0,0,0; v,1-v,v,0,0,0; v,v,1-v,0,0,0; 0,0,0,(1-2*v)/2,0,0; 0,0,0,0,(1-2*v)/2,0; 0,0,0,0,0,(1-2*v)/2]');
sys.addMatrix('K', 'B''*D*B');

elem = mesh.element(1);
N = elem.shape([0,0.25,0])

B = elem.shapeDeriv([0,0.25,0])



