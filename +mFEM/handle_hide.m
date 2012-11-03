classdef handle_hide < handle
    %HANDLE_HIDE A wrapper class causing methods of handle to be hidden.
    %
    % Syntax
    %   classdef MyHandleClass < handle_hide
    %
    % Description
    %   classdef MyHandleClass < handle_hide makes MyHandleClass a subclass
    %   of the handle class, but hides all the methods, thus when the
    %   documentation is created the handle based methods to not appear.
    %
    % see also HANDLE ELEMENT FEMESH SYSTEM
    %
    %----------------------------------------------------------------------
    % Copyright 2012 Andrew E. Slaughter
    % This software is for educational purposes only and may not be used
    % without written permession.
    %----------------------------------------------------------------------
    methods(Sealed)  
      function TF = eq(varargin)
         TF = eq@handle(varargin{:});
      end
    end
    methods(Hidden)
      function lh = addlistener(varargin)
         lh = addlistener@handle(varargin{:});
      end
      function notify(varargin)
         notify@handle(varargin{:});
      end
      function delete(varargin)
         delete@handle(varargin{:});
      end
      function Hmatch = findobj(varargin)
         Hmatch = findobj@handle(varargin{:});
      end
      function p = findprop(varargin)
         p = findprop@handle(varargin{:});
      end
      function TF = ne(varargin)
         TF = ne@handle(varargin{:});
      end
      function TF = lt(varargin)
         TF = lt@handle(varargin{:});
      end
      function TF = le(varargin)
         TF = le@handle(varargin{:});
      end
      function TF = gt(varargin)
         TF = gt@handle(varargin{:});
      end
      function TF = ge(varargin)
         TF = ge@handle(varargin{:});
      end
   end
end