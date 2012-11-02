function opt =  gather_user_options(opt,varargin)
% GATHER_USER_OPTIONS collects property parings input into a function
%__________________________________________________________________________
% SYNTAX:
%   opt = gatheruseroptions(opt,'PropertyName',<PropertyValue>,...);
%   opt = gatheruseroptions(...,-PropertyName);
%   opt = gatheruseroptions(opt,optUser);
%
% DESCRIPTION:
%   opt = gatheruseroptions(opt,'PropertyName',<PropertyValue>,...)
%       compares the property names with the default values stored in the
%       data structure and applies the property value if the field matches
%       the property name.
%
%   opt = gatheruseroptions(...,-PropertyName)
%       adds an alternative method for flipping boolean values without
%       specifing the value. For example if the options structure is:
%           opt.flag = true;
%       The the following two commands are equivalent
%           opt = gather_user_options(opt,'-flag');
%           opt = gather_user_options(opt,'flag',false);
%
%   opt = gatheruseroptions(opt,optUser) in this case the user suplies a
%       data structure similar to that of opt, which is compared according
%       to the fieldnames
%
% INPUT:
%   opt = a data structure containing the default values, e.g. opt =
%       struct('Prop1',true,'Prop2','install','Coeffient',1);.
%__________________________________________________________________________

% INTILIZE THE DATA
    q = varargin;           % User supplied input
    k = 1;                  % Intilize the counter
    N = nargin - 1;         % Number of property inputs

% CONVERT DATA INPUT AS A STRUCTURE TO A CELL ARRAY
if length(q) == 1 && isstruct(q{1});
    S = q{1}; q = {};
    fn = fieldnames(S); % Fieldnames from input structure
    c = struct2cell(S); % Corresponding setings desired
    n = length(fn);     % The number if user inputs 
    q(1:2:2*n-1) = fn;  % Inerts the 'PropName'
    q(2:2:2*n) = c;     % Inserts the PropValue
    N = length(q);      % Number of inputs
end

% COMPARE INPUT WITH THE DEFAULTS
while k < N
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
           disp('changing the value')
            value = ~opt.(itm); 
            
       % Produce a warning and move on, if the item is not recognized
       else
            mes = ['The option, ',itm,', was not recoignized.'];
            warning(mes);
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
        mes = ['The option, ',itm,', was not recoignized.'];
        warning(mes);
    end
end


