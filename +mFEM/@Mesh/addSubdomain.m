function addSubdomain(obj, id, varargin)
    %ADDSUBDOMAIN Labels elements with numeric tags.
    %
    % Syntax
    %   addSudomain(id, Limit1,...)
    %
    % Description
    %   addSudomain(id, Limit1,...) allows user to seperate the
    %   domain into subdomains based on a variety of spatial
    %   criteria functions, see the description below.
    %
    % Descriptions of Limits
    %   Functions:
    %       string | cell array of strings
    %       It is possible to tag domains using test functions.
    %       For example, 'x==1' will tag elements that have a
    %       node location at 1. Mutliple conditions should be
    %       expressed in a cell array. For example, {'x==1','y<2'}
    %       will id elements with nodes at x = 1 AND y < 2.

    % Check if system is initialized
    if ~obj.initialized;
        error('FEmesh:addSubdomain:NonInitializedMesh',...
            'The FEmesh object must be initialized');
    end

    % Apply the subdomain flags
    obj.addTag(id,'subdomain',varargin{:});
end