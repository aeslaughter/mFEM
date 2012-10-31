%TMESSAGE Function for displaying time messages
%
% Syntax
%   ticID = tmessage(X);
%   ticID = tmessage(formatSpec,A1,...,An)
%   tmessage(ticID)

function varargout = tmessage(varargin)

% Start a message
if nargin >= 1 && ischar(varargin{1});
    str = sprintf(varargin{:});
    disp(str);
    varargout{1} = tic;
    
% Stop a message
else
    disp(['   Done in ', num2str(toc(varargin{1})),' sec.']);
end
    
    