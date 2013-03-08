function nodes = getNodes(obj,varargin)
    %GETNODES Returns node objects
    %   Each node is a class, this allows for the extraction of these node
    %   objects using the global id and/or a combination properties. If
    %   operating in parallel this will perform all the necessary
    %   communication and gather the nodes.
    %
    % Syntax
    %   no = getNodes()
    %   no = getNodes(id)
    %   no = getNodes(...,'PropertyName',PropertyValue,...)
    %
    % Description
    %   no = getNodes() returns the node objects for the entire mesh.
    %
    %   no = getNodes(id) returns the nodes corresponding to the global ids
    %   supplied in id.
    %
    %   no = getNodes(...,'PropertyName',PropertyValue,...) returns the 
    %   node objects as above but limits the selection further based on the
    %   property parings supplied, which are detailed below.
    %
    % GETNODES Property Descriptions
    %   Tag
    %       scalar | char
    %       Limits the node objects returned to only those with the
    %       supplied tag, the term tag is applied to both those added with
    %       the addBoundary and addSubdomain commands.
    %
    %   Lab
    %       scalar | vector
    %       Limits the node objects returned to the given processors in
    %       parallel applications
    %
    % See Also addBoundary addSubdomain getElements
    %
    %----------------------------------------------------------------------
    %  mFEM: A Parallel, Object-Oriented MATLAB Finite Element Library
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
    %  along with this program. If not, see <http://www.gnu.org/licenses/>.
    %
    %  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
    %----------------------------------------------------------------------

    % Collect the nodes
    nodes = obj.gatherComposite(varargin{:},'name','nodes');
end