function elem = getElements(obj, varargin)
    %GETELEMENTS returns the desired elements
    %   Each element is a class, this allows for the extraction of these
    %   objects using the global id and/or a combination properties. If
    %   operating in parallel this will perform all the necessary
    %   communication and gather the elements.
    %
    % Syntax
    %   elem = getElements()
    %   elem = getElements(id)
    %   elem = getElements(id,'PropertyName',PropertyValue,...)
    %
    % Description
    %   elem = getElement(id) returns an array of elements with the ids
    %   given the numeric vector id.
    %
    %   elem = getElement(id,lab) operates the same as above but limits the
    %   search to the specified processor given in lab.
    %
    % GETELEMENTS Property Descriptions
    %   Tag
    %       scalar | char
    %       Limits the objects returned to only those with the
    %       supplied tag, the term tag is applied to both those added with
    %       the addBoundary and addSubdomain commands.
    %
    %   Lab
    %       scalar | vector
    %       Limits the objects returned to the given processors in
    %       parallel applications, the objects are automatically gathered
    %       to the calling lab.
    %
    %   Gather
    %       {false} | true
    %       If true the objects stored in parallel are gathered to the
    %       calling lab.
    %
    % See Also addBoundary addSubdomain getNodes
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
    %---------------------------------------------------------------------

    elem = obj.gatherComposite(varargin{:},'name','elements');
end