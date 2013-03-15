function cmp = parseComponentInput(cmp)
    %PARSECOMPONENTINPUT (protected) preps getDof component input for use
    %   The GETDOF method allows for certain components of the degrees-of-
    %   freedom, this function parse the input for use by GETDOF.
    %
    % Syntax
    %   cmp = parseComponentInput(cmp)
    %
    % Description
    %   cmp = parseComponentInput(cmp) parse the cmp input for preperation
    %   for use by GETDOF, see the help for getDof the acceptable forms for
    %   the cmp input.
    %
    % See Also getDof
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
   
    % Define the error message, which are used in various locations
    msgid = 'Node:getDof:InvalidComponent';
    errmsg = 'The supplied component must be number or ''x'', ''y'', or ''z''.';

    % Create a cell array for single char string input
    if ischar(cmp); 
        cmp = {cmp};
    end
    
    % Convert string input into numerical
    if iscell(cmp);
        for i = 1:length(cmp);
            if ischar(cmp{i});
                switch lower(cmp{i});
                    case 'x'; cmp{i} = 1;
                    case 'y'; cmp{i} = 2;
                    case 'z'; cmp{i} = 3;
                    otherwise
                        error(msgid,errmsg);
                end
            elseif ~isnumeric(cmp{i});
                error(msgid,errmsg);
            end
        end
        
        % Convert to numeric array
        cmp = cell2mat(cmp);
    
    % If not a numeric input case, throw an error
    elseif ~isnumeric(cmp);
        error(msgid,errmsg);
    end
end
    