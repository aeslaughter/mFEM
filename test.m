function test

clear;
import mFEM.elements.*

% elem = Line2(1,[0,0;1,1],'space','vector');
% 
% x = elem.get_position(1)

elem = Beam(1,[0;1],'space','vector');

x = elem.get_position(0,'index',[1,3])
x = elem.get_position(0)