function opt =  gatheruseroptions(opt,varargin)
% GATHERUSEROPTIONS collects property parings input into a function
%__________________________________________________________________________
% SYNTAX:
%   opt = gatheruseroptions(opt,'PropertyName',<PropertyValue>,...);
%   opt = gatheruseroptions(opt,optUser);
%
% DESCRIPTION:
%   opt = gatheruseroptions(opt,'PropertyName',<PropertyValue>,...)
%       compares the property names with the default values stored in the
%       data structure and applies the property value if the field matches
%       the property name.
%   opt = gatheruseroptions(opt,optUser) in this case the user suplies a
%       data structure similar to that of opt, which is compared according
%       to the fieldnames
%
% INPUT:
%   opt = a data structure containing the default values, e.g. opt =
%       struct('Prop1',true,'Prop2','install','Coeffient',1);.
%__________________________________________________________________________

% 1 - INTILIZE THE DATA
    q = varargin;           % User supplied input
    list = fieldnames(opt); % Fieldnames of default values
    k = 1;                  % Intilize the counter
    N = nargin - 1;         % Number of inputs

% 2 - CONVERT DATA INPUT AS A STRUCTURE TO A CELL ARRAY
if length(q) == 1 && isstruct(q{1});
    S = q{1}; q = {};
    fn = fieldnames(S); % Fieldnames from input structure
    c = struct2cell(S); % Corresponding setings desired
    n = length(fn);     % The number if user inputs 
    q(1:2:2*n-1) = fn;  % Inerts the 'PropName'
    q(2:2:2*n) = c;     % Inserts the PropValue
    N = length(q);      % Number of inputs
end

% 3 - COMPARE INPUT WITH THE DEFAULTS
while k < N
    % Seperate the name from the value; increment counter
    itm = q{k}; value = q{k+1}; k = k + 2;
    
    % When the property matches, update the structure
    if strmatch(lower(itm),list,'exact');   
        opt.(lower(itm)) = value;
        
    % Produce a warning if the property is not recongnized    
    else 
        mes = ['The option, ',itm,', was not recoignized.'];
        disp(mes);
    end
end