%EIG_GCOH Performs eigendecomposition for a cross-spectrum -- unimplemented
%Optional file header info (to give more details about the function than one line)
%Optional file header info (to give more details about the function than one line)
%
% Usage:  
%    [output1,output2] = eig_gcoh(input1)
%    [output1,output2] = eig_gcoh(input1,input2)
%
% Inputs:
%    input1 - Description
%    input2 - Description
%
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Example: 
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Other requirements: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%
% Copyright Apr-2020, David Zhou, dwzhou@mit.edu
% Last revision 06-Apr-2020
%------------------------------------------------
 
function [outputs,varargout] = eig_gcoh(inputs,varargin)

%---------------------------------------
% PROCESS INPUTS
%---------------------------------------
 
%VARARGIN
switch length(varargin)
    case 1
        newvar1=varargin{1};
    case 2
        newvar1=varargin{1};
        newvar2=varargin{2};
   otherwise
        newvar1=[];
        newvar2=[];
end

%VARARGOUT
switch length(varargout)
    case 1
        execute_procedure1=1;
    case 2
        execute_procedure1=1;
        execute_procedure2=1;
end

%INPUT ERRORS
if inputerror1
    error('eig_gcoh:inputError','Error:  Message goes here.');
end

%GLOBAL VARIABLES
globalvar1=some_value;
globalvar2=some_other_value;

%OTHER INPUT PROCESSING
    %insert here


%---------------------------------------
% MAIN BODY OF FUNCTION
%---------------------------------------

try 
    %PROCEDURE 1 NAME
    %code goes here
catch err
        % issue a warning
    warning('eig_gcoh:warningType',...
        ['Warning: The following went wrong in Procedure 1 Name --\n' err.message ]);
    %error handling code here

        % or issue an error & exit function:
    error('eig_gcoh:errorType',...
        ['Error: The following went wrong in Procedure 1 Name --\n' err.message ]);
end

try
    %PROCEDURE 2 NAME
    %code goes here
catch err
        % issue a warning
    warning('eig_gcoh:warningType',...
        ['Warning: The following went wrong in Procedure 2 Name --\n' err.message ]);
    %error handling code here

        % or issue an error & exit function:
    error('eig_gcoh:errorType',...
        ['Error: The following went wrong in Procedure 2 Name --\n' err.message ]);
end


%---------------------------------------
% SUBFUNCTIONS
%---------------------------------------


%---------------------------------------
% END CODE
%---------------------------------------
end
