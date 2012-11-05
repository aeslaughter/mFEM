function opt =  gather_user_options(opt,varargin)
    % Collects property parings input into a function
    %
    % Syntax
    %   opt = gather_user_options(opt,'PropertyName',<PropertyValue>,...);
    %   opt = gather_user_options(...,-PropertyName,...);
    %   opt = gather_user_options(opt,optUser);
    %
    % Description
    %   opt = gather_user_options(opt,'PropertyName',<PropertyValue>,...)
    %   compares the property names with the default values stored in the
    %   data structure and applies the property value if the field matches
    %   the property name. The opt input variable is described below.
    %
    %   opt = gather_user_options(...,-PropertyName)
    %   adds an alternative method for flipping boolean values without
    %   specifing the value. For example if the options structure is:
    %       opt.flag = true;
    %   The the following two commands are equivalent.
    %       opt = gather_user_options(opt,'-flag');
    %       opt = gather_user_options(opt,'flag',false);
    %
    %   opt = gather_user_options(opt,optUser) in this case the user suplies a
    %   data structure similar to that of opt, which is compared according
    %   to the fieldnames.
    %
    % Input
    %   opt = a data structure containing the default values, e.g. opt =
    %       struct('Prop1',true,'Prop2','install','Coeffient',1);. 
    %
    %----------------------------------------------------------------------
    % Copyright 2012 Andrew E. Slaughter
    % This software is for educational purposes only and may not be used
    % without written permession.
    %----------------------------------------------------------------------

% Intilize the data
    q = varargin;           % User supplied input
    k = 1;                  % Intilize the counter
    N = nargin - 1;         % Number of property inputs

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
