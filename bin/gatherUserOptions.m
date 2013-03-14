function [opt, unknowns] =  gatherUserOptions(opt, varargin)
    % Collects property parings input into a function
    %
    % Syntax
    %   opt = gatherUserOptions(opt,'PropertyName',PropertyValue,...);
    %   opt = gatherUserOptions(...,-PropertyName,...);
    %   opt = gatherUserOptions(...,'GatherUserOptions',{'PropName', PropValue,...});
    %   opt = gatherUserOptions(opt,optUser);
    %   [opt, unknown] = gatherUserOptions(...);
    %
    % Description
    %   opt = gatherUserOptions(opt,'PropertyName',<PropertyValue>,...)
    %   compares the property names with the default values stored in the
    %   data structure and applies the property value if the field matches
    %   the property name. The opt input variable is a data structure 
    %   containing the default values, 
    %       e.g. opt = struct('Prop1',true,'Prop2','install','Coeffient',1);
    %
    %   opt = gatherUserOptions(...,-PropertyName)
    %   adds an alternative method for flipping boolean values without
    %   specifing the value. For example if the options structure is:
    %       opt.flag = true;
    %   The the following two commands are equivalent.
    %       opt = gatherUserOptions(opt,'-flag');
    %       opt = gatherUserOptions(opt,'flag',false);
    %
    %   opt = gatherUserOptions(..., 'GatherUserOptions',{'PropName', PropValue,...}) options
    %   may be passed to this function itself by entering the pairs in a
    %   cell array that occupys the last input. 
    %
    %   opt = gatherUserOptions(opt,optUser) in this case the user suplies a
    %   data structure similar to that of opt, which is compared according
    %   to the fieldnames.
    %
    %   [opt, unknown] = gatherUserOptions(...) also returns a cell array
    %   of property pairs for the unknown options, i.e. the options not
    %   included in the opt input structure.
    %
    % gatherUserOptions Property Descriptions
    %   DisableWarn
    %       true | {false}
    %       Toggles the display of the unknown option warning.
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

% Automatically disable warnings if unknowns are output
options.disablewarn = false;
if nargout == 2;
    options.disablewarn = true;
end

% Apply options for this function
for i = 1:length(varargin)
    if ischar(varargin{i}) && strcmpi(varargin{i},'GatherUserOptions');
        input = varargin{i+1};
        options = gatherUserOptions(options, input{:});
        varargin = varargin(1:i-1);
        break;
    end
end

% Intilize the data
q = varargin;           % User supplied input
k = 1;                  % Intilize the counter
u = 1;                  % unknown counter
unknowns = {};          % unknown output cell array
N = length(q);          % Number of property inputs

% Convert data input as a structure to a cell array
if length(q) == 1 && isstruct(q{1});
    S = q{1}; q = {};
    fn = fieldnames(S); % Fieldnames from input structure
    c = struct2cell(S); % Corresponding setings desired
    n = length(fn);     % The number if user inputs 
    q(1:2:2*n-1) = fn;  % Inerts the 'PropName'
    q(2:2:2*n) = c;     % Inserts the PropValue
    N = length(q);      % Number of inputs
end

% Compare input with the defaults
while k <= N
    % Seperate the name
    itm = q{k};
    
    % Search for flag style input
    idx = strfind(itm,'-');
    
    % Case when a precedding '-' is detected
    if numel(idx) > 0 && idx(1) == 1;

       itm = lower(itm(2:end)); % remove the leading '-'
       k = k + 1;        % increment the counter
       
       % Change the value, if it existing in the options structure
       if isfield(opt, itm) 
           value = ~opt.(itm); 
            
       % Produce a warning and move on, if the item is not recognized
       else
           unknowns{u} = itm; 
           unknowns{u+1} = true;
           u = u + 2;
             if ~options.disablewarn;
                mes = ['The option, ', itm,', was not recognized and is being ignored.'];
                warning('gatherUserOptions:UnknownProperty', mes);
            end
            continue;
       end
       
    % Standard case, the next item is the value
    else
       value = q{k+1}; k = k + 2;
    end
 
    % When the property matches, update the structure
    if isfield(opt, lower(itm));   
        opt.(lower(itm)) = value;

    % Produce a warning if the property is not recongnized    
    else
        unknowns{u} = itm; 
        unknowns{u+1} = value;
        u = u + 2;
        if ~options.disablewarn;
            mes = ['The option, ',itm,', was not recognized and is being ignored.'];
            warning('gatherUserOptions:UnknownProperty', mes);
        end
    end
end
