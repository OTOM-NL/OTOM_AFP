function varargout = ADS_variables_window(varargin)
% ADS_VARIABLES_WINDOW MATLAB code for ADS_variables_window.fig
%      ADS_VARIABLES_WINDOW, by itself, creates a new ADS_VARIABLES_WINDOW or raises the existing
%      singleton*.
%
%      H = ADS_VARIABLES_WINDOW returns the handle to a new ADS_VARIABLES_WINDOW or the handle to
%      the existing singleton*.
%
%      ADS_VARIABLES_WINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADS_VARIABLES_WINDOW.M with the given input arguments.
%
%      ADS_VARIABLES_WINDOW('Property','Value',...) creates a new ADS_VARIABLES_WINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ADS_variables_window_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ADS_variables_window_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ADS_variables_window

% Last Modified by GUIDE v2.5 30-Aug-2017 10:41:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ADS_variables_window_OpeningFcn, ...
                   'gui_OutputFcn',  @ADS_variables_window_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ADS_variables_window is made visible.
function ADS_variables_window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ADS_variables_window (see VARARGIN)

% Choose default command line output for ADS_variables_window
handles.output = hObject;
javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ADS_variables_window wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ADS_variables_window_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
