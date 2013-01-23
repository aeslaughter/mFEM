function varargout = tmessage(varargin)
    % Function for displaying time messages.
    %
    % Syntax
    %   ticID = tmessage(X)
    %   ticID = tmessage(formatSpec,A1,...,An)
    %   tmessage(ticID)
    %
    % Description
    %   ticID = tmessage(X) displays the string X and starts recoreding the
    %   exectution time (see doc tic)
    %
    %   ticID = tmessage(formatSpec,A1,...,An) same as above put allows for
    %   string to accept additional inputs, use format specified by the MATLAB
    %   sprintf function.
    %
    %   tmessage(ticID) completes the time message.
    %
    %----------------------------------------------------------------------
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
    %----------------------------------------------------------------------

% Start a message
if nargin >= 1 && ischar(varargin{1});
    str = sprintf(varargin{:});
    disp(str);
    varargout{1} = tic;
    
% Stop a message
else
    disp(['   Done in ', num2str(toc(varargin{1})),' sec.']);
end
    
    