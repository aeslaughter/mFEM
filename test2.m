function test2
clear;
import mFEM.*;

mesh = FEmesh('Quad4','scalar','DG');
mesh.grid(0,1,0,1,100,100);

% E = mesh.element(1:2)


