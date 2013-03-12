%% Example 1c: Reproduces Example 5.1 of Fish & Belytschko (2007)
% This example uses kernels directly.
% 
% Syntax
%   example1c
%   example1c(n)
%
% Description
%   example1c runs example exactly as done in textbook, with two elements
%   example1c(n) runs the example with n number of elements
%
% See also EXAMPLE1A EXAMPLE1B
%
%--------------------------------------------------------------------------
%  mFEM: An Object-Oriented MATLAB Finite Element Library
%  Copyright (C) 2012 Andrew E Slaughter
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
%--------------------------------------------------------------------------

%% Initilization
% Functin declaration
function T = example1c(varargin)

% Parse optional input
nel = 2;
if nargin == 1 && isnumeric(varargin{1}); 
    nel = varargin{1};
end

%% Create a FEmesh Object 
% Build mesh of 2-node linear elemnts from 0 to 4
mesh = mFEM.Mesh();
mesh.grid('Line2',0,4,nel); 
mesh.init();

%% Label The Boundaries
mesh.addBoundary(1,'left');    % T = 0 boundary (essential)    
mesh.addBoundary(2,'right');   % q = 20 boundary   

%% Define Constants
k = 2;          % thermal conductivity 
A = 0.1;        % cross sectional area
b = 5;          % heat source (defined over entire domain)
q_bar = 5;      % right boundary prescribed heat flux
T_bar = 0;      % known temperatures

%% Define Kernels
diffusion = mFEM.kernels.Diffusion(mesh, 'D', k*A);
source = mFEM.kernels.Source(mesh, 'b', b);
% flux = Source(mesh, 'b', -q_bar*A, 'boundary', 2);
flux = mFEM.kernels.Source(mesh, 'b', -q_bar*A);

%% Assemble Stiffness Matrix and Force Vector
K = diffusion.assemble();
% f = source.assemble() + flux.assemble();
f = source.assemble() + flux.assemble('Tag', 2);

%% Define Variables for Essential and Non-essential Degrees-of-freedom
ess = mesh.getDof('Tag', 1);  % 1
non = ~ess;                         % 2,3

%% Solve for the Temperatures
T = zeros(size(f));         % initialize the temperature vector
T(ess) = T_bar;             % apply essential boundary condtions
T(non) = K(non,non)\f(non); % solve for T on the non-essential boundaries

%% Compute the Temperature Gradients
% Loop through the elements
elements = mesh.getElements();
for e = 1:length(elements);
    
    % Extract the current element from the mesh object
    elem = elements(e);
    
    % Collect the local values of T
    d(:,1) = T(elem.getDof());
    
    % Compute the temperature gradient at the gauss point, store the value
    % twice for each element for creating graph, TGx is the node locations
    % used for plotting
    TG(1:2,e) = elem.shapeDeriv(elem.qp(1))*d;
    TGx(1:2,e) = elem.nodes.getCoord();
end    

%% Generate Figure for T and TG Solutions

% Create Exact Solutions
x0 = 0:0.1:4;
Tex = -12.5*x0.^2 + 97.5*x0;
TGex = -25*x0 + 97.5;

% Initilize the figure
figure('Color','w','Name','Example 1 Results');

% Create Temperature Plot
h = subplot(2,1,1);
mesh.plot(T,'-ShowNodes'); hold on;
plot(h,x0,Tex,'k-','LineWidth',1);
legend({'FEM','Exact'},'location','best');
xlabel('x (m)','interpreter','tex');
ylabel('Temperature (\circC)','interpreter','tex');

% Create TG Plot
h = subplot(2,1,2);
plot(h,x0,TGex,'k-',TGx,TG,'b-o','LineWidth',1);
legend({'Exact','FEM'},'location','best');
xlabel('x (m)','interpreter','tex');
ylabel('Temp. Gradient (\circC/m)','interpreter','tex');
