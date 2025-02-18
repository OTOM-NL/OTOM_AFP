function varargout = Inlinemonitoring_control_window(varargin)
%INLINEMONITORING_CONTROL_WINDOW MATLAB code file for Inlinemonitoring_control_window.fig
%      INLINEMONITORING_CONTROL_WINDOW, by itself, creates a new INLINEMONITORING_CONTROL_WINDOW or raises the existing
%      singleton*.
%
%      H = INLINEMONITORING_CONTROL_WINDOW returns the handle to a new INLINEMONITORING_CONTROL_WINDOW or the handle to
%      the existing singleton*.
%
%      INLINEMONITORING_CONTROL_WINDOW('Property','Value',...) creates a new INLINEMONITORING_CONTROL_WINDOW using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to Inlinemonitoring_control_window_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      INLINEMONITORING_CONTROL_WINDOW('CALLBACK') and INLINEMONITORING_CONTROL_WINDOW('CALLBACK',hObject,...) call the
%      local function named CALLBACK in INLINEMONITORING_CONTROL_WINDOW.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Inlinemonitoring_control_window

% Last Modified by GUIDE v2.5 17-Aug-2018 14:41:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Inlinemonitoring_control_window_OpeningFcn, ...
                   'gui_OutputFcn',  @Inlinemonitoring_control_window_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before Inlinemonitoring_control_window is made visible.
function Inlinemonitoring_control_window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for Inlinemonitoring_control_window
handles.output = hObject;
javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Inlinemonitoring_control_window wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Inlinemonitoring_control_window_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_Time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Time as text
%        str2double(get(hObject,'String')) returns contents of edit_Time as a double


% --- Executes during object creation, after setting all properties.
function edit_Time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_inc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_inc as text
%        str2double(get(hObject,'String')) returns contents of edit_inc as a double


% --- Executes during object creation, after setting all properties.
function edit_inc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_Inline_code.
function checkbox_Inline_code_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Inline_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Inline_code



function edit_init_Temp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_init_Temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_init_Temp as text
%        str2double(get(hObject,'String')) returns contents of edit_init_Temp as a double


% --- Executes during object creation, after setting all properties.
function edit_init_Temp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_init_Temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_Live_output.
function checkbox_Live_output_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Live_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Live_output


% --- Executes on button press in checkbox_Video.
function checkbox_Video_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Video
