% Example 5.1 of Fish & Belytschko (2007).
function example1

% Clear all variables, including classes
clear;

% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, 2 node linear mesh from 0 to 4
mesh = FEmesh('Linear2');
mesh.grid(0,4,2); % the mesh is initialized automatically with grid function


