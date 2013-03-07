function elem = getElements(obj, varargin)
    %GETELEMENT returns the desired elements
    %
    % Syntax
    %   elem = getElement(id)
    %   elem = getElement(id,lab)
    %
    % Description
    %   elem = getElement(id) returns an array of elements with the ids
    %   given the numeric vector id.
    %
    %   elem = getElement(id,lab) operates the same as above but limits the
    %   search to the specified processor given in lab.
    %
    %----------------------------------------------------------------------
    %  mFEM: An Object-Oriented MATLAB Finite Element Library
    %  Copyright (C) 2013 Andrew E Slaughter
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

    % Parse input
    if nargin == 1;
        id = 1:obj.n_elements;
    elseif nargin == 2;
        id = varargin{1};
        lab = [];
    else
        id = varargin{1};
        lab = varargin{2}; 
    end
    
    elem = obj.gatherComposite('elem',id,lab); 
end




