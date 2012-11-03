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
    % Copyright 2012 Andrew E. Slaughter
    % This software is for educational purposes only and may not be used
    % without written permession.
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
    
    