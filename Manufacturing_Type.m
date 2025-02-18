function varargout = Manufacturing_Type(varargin)


% Attention for GUI update
% Save shode GUI ro yek bar Run begir Save kon, badesh mishe estefade
% kard!!



% MANUFACTURING_TYPE MATLAB code for Manufacturing_Type.fig
%      MANUFACTURING_TYPE, by itself, creates a new MANUFACTURING_TYPE or raises the existing
%      singleton*.
%
%      H = MANUFACTURING_TYPE returns the handle to a new MANUFACTURING_TYPE or the handle to
%      the existing singleton*.
%
%      MANUFACTURING_TYPE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANUFACTURING_TYPE.M with the given input arguments.
%
%      MANUFACTURING_TYPE('Property','Value',...) creates a new MANUFACTURING_TYPE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Manufacturing_Type_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Manufacturing_Type_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Manufacturing_Type

% Last Modified by GUIDE v2.5 23-Aug-2019 11:54:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Manufacturing_Type_OpeningFcn, ...
                   'gui_OutputFcn',  @Manufacturing_Type_OutputFcn, ...
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


% --- Executes just before Manufacturing_Type is made visible.
function Manufacturing_Type_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Manufacturing_Type (see VARARGIN)

% Choose default command line output for Manufacturing_Type


C = imread('Winding_direction.png');
image(C);
axis off;

% set(gcf,'color','w');
javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));


handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Manufacturing_Type wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Manufacturing_Type_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_delayTime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_delayTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_delayTime as text
%        str2double(get(hObject,'String')) returns contents of edit_delayTime as a double


% --- Executes during object creation, after setting all properties.
function edit_delayTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_delayTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pitchAngle_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pitchAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pitchAngle as text
%        str2double(get(hObject,'String')) returns contents of edit_pitchAngle as a double


% --- Executes during object creation, after setting all properties.
function edit_pitchAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pitchAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Tape Global Temperature

[file,path] = uigetfile('.\Tape_Temp_Global.txt','Global Temperature Tape');
% fileID3 = fopen('Comp_Parmeter.txt','w');
% fileID1 = fopen(strcat(path,file),'w');
set(handles.edit_globalTape,'String',strcat(path,file));


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Substrate Global Temperature

[file,path] = uigetfile('.\Sub3D_Temp_Global.txt','Global Temperature Substrate');
% fileID3 = fopen('Comp_Parmeter.txt','w');
% fileID1 = fopen(strcat(path,file),'w');
set(handles.edit_globalSub,'String',strcat(path,file));



function edit_LocationNo_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LocationNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_LocationNo as text
%        str2double(get(hObject,'String')) returns contents of edit_LocationNo as a double


% --- Executes during object creation, after setting all properties.
function edit_LocationNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LocationNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_Wdir.
function radiobutton_Wdir_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_Wdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_Wdir


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Tape initial Temperature
[file,path] = uigetfile('.\T_Final_Tape.txt','Initial Temperature Tape');
% fileID3 = fopen('Comp_Parmeter.txt','w');
% fileID1 = fopen(strcat(path,file),'w');
set(handles.edit_initial_Tape,'String',strcat(path,file));



% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Substrate initial Temperature

[file,path] = uigetfile('.\T_Final_Sub_3D.txt','Initial Temperature Substrate');
% fileID3 = fopen('Comp_Parmeter.txt','w');
% fileID1 = fopen(strcat(path,file),'w');
set(handles.edit_initial_sub,'String',strcat(path,file));



function edit_initial_Tape_Callback(hObject, eventdata, handles)
% hObject    handle to edit_initial_Tape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_initial_Tape as text
%        str2double(get(hObject,'String')) returns contents of edit_initial_Tape as a double


% --- Executes during object creation, after setting all properties.
function edit_initial_Tape_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_initial_Tape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_initial_sub_Callback(hObject, eventdata, handles)
% hObject    handle to edit_initial_sub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_initial_sub as text
%        str2double(get(hObject,'String')) returns contents of edit_initial_sub as a double


% --- Executes during object creation, after setting all properties.
function edit_initial_sub_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_initial_sub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_globalTape_Callback(hObject, eventdata, handles)
% hObject    handle to edit_globalTape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_globalTape as text
%        str2double(get(hObject,'String')) returns contents of edit_globalTape as a double


% --- Executes during object creation, after setting all properties.
function edit_globalTape_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_globalTape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_globalSub_Callback(hObject, eventdata, handles)
% hObject    handle to edit_globalSub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_globalSub as text
%        str2double(get(hObject,'String')) returns contents of edit_globalSub as a double


% --- Executes during object creation, after setting all properties.
function edit_globalSub_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_globalSub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function text4_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to text4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function text3_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to text3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function text5_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to text5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function text2_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_Cooling.
function pushbutton_Cooling_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Cooling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% get(handles.popupmenu_Wdir,'Value'); %1 or 2
% Usecooling=get(handles.radiobutton_use_cooling,'Value');


  fileID_input_reader = fopen('.\Inputs_mode.txt');
        [Inputs_Vel_power] = textscan(fileID_input_reader,'%f %f','delimiter',',');
        % should be modified later for general cases
        Vel_Tank=cell2mat(Inputs_Vel_power(1,1));%/1000;  % should
  
    V_ave=mean(Vel_Tank);


 
    
Tape_global=get(handles.edit_globalTape,'String');
Sub_global=get(handles.edit_globalSub,'String');
Analysislog_global=get(handles.edit_analysis_log,'String');
  
Inner_gas_rho_cp=str2num(get(handles.edit_air_in_pr,'String'));


          T_nip_Sub = dlmread(Sub_global) ;
     T_nip_Tape = dlmread(Tape_global) ; 
     
     
     
     % Use previous Gas inside
     Use_inGasLast=get(handles.radiobutton_inGasLast,'Value'); 
     
     if Use_inGasLast==1
         
     findID_Inner_gas=  fopen('.\cooling\Temp_Inner_gas.txt','r');
     if findID_Inner_gas>0
     
     Temp_gas = textscan(findID_Inner_gas,'%f','delimiter',',','HeaderLines',1);
     fclose(findID_Inner_gas);
     
    
    
     T_in=mean(Temp_gas{:});
      if isnan(T_in)
          
          T_in=20;  
     end
    
     else
         % if the file not existed!
       T_in=20;   
    
     end
     
     else
     % if no use of previous gas was chosen!
      T_in=20; 
     
     end
     
      T_out=20;
     T_out_ingas=[T_out T_in];
     
     
    
     fileID_90=fopen(Analysislog_global,'r'); 
     
Time_last_layer = textscan(fileID_90,'  %f  ','Delimiter',',','HeaderLines',2) ;
Time_last_layer=cell2mat(Time_last_layer);
 fclose(fileID_90);
 
 
 fid13 = fopen('.\Supp_files\Transient_ID.txt','r');
    Tr_ID = textscan(fid13,'%f','delimiter',',');
 fclose(fid13);
 
 
 fid10 = fopen('.\Geo_Parameter.txt');
  out = textscan(fid10,'%s','delimiter',',');
    defaultans_geo = {out{1}{2},out{1}{4},out{1}{6},out{1}{8},out{1}{10}...
        out{1}{12},out{1}{14},out{1}{16},out{1}{18},...
        out{1}{20},out{1}{22},out{1}{24},out{1}{26},...
        out{1}{28},out{1}{30},out{1}{32},out{1}{34}};
    
    fclose(fid10);
    
    
    
fid11 = fopen('.\Process_Parameter.txt');

    out = textscan(fid11,'%s','delimiter',',');
    defaultans_process = {out{1}{2},out{1}{4},out{1}{6},out{1}{8},out{1}{10}...
        out{1}{12},out{1}{14},out{1}{16},out{1}{18},...
        out{1}{20},out{1}{22}};
    
    fclose(fid11);
    
    %% 22 Aug 2019
    V_ave=str2num(defaultans_process{3});
    
       H=str2num(defaultans_process{9});
 
 Layer_number=1;
R_cyl= str2num(defaultans_geo{11});


%% Needs to be modified !!
materials_tape= str2num(defaultans_process{1}) ;  % material of the Tape
materials= str2num(defaultans_process{2}) ;  % material of the Sub
Lx=str2num( defaultans_geo {2}); % is the width of Tape
Ly= str2num(defaultans_geo{10}); % is the thickness   >> Tape defaultans(3)+
z_cyl_end= str2double(defaultans_geo{12});

th_T=str2num(defaultans_geo{3});

Wind_dir=get(handles.popupmenu_Wdir,'Value'); 

delay_layers=get(handles.edit_delayTime,'Value'); 
Time_transient=Tr_ID{1} (5);



    
    Temp_Cooling_pipe(Layer_number,V_ave,T_nip_Tape,T_nip_Sub,R_cyl(1),materials, Lx,Ly,H,Wind_dir,T_out_ingas,Time_last_layer,delay_layers,Time_transient,z_cyl_end,th_T,Inner_gas_rho_cp,...
                      materials_tape)






% --- Executes on button press in pushbutton_Log.
function pushbutton_Log_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile('.\Analysis_log.txt','Analysis Log of previous layer');
% fileID3 = fopen('Comp_Parmeter.txt','w');
% fileID1 = fopen(strcat(path,file),'w');
set(handles.edit_analysis_log,'String',strcat(path,file));



function edit_analysis_log_Callback(hObject, eventdata, handles)
% hObject    handle to edit_analysis_log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_analysis_log as text
%        str2double(get(hObject,'String')) returns contents of edit_analysis_log as a double



% --- Executes during object creation, after setting all properties.
function edit_analysis_log_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_analysis_log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_Wdir.
function popupmenu_Wdir_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Wdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Wdir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Wdir


% --- Executes during object creation, after setting all properties.
function popupmenu_Wdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_Wdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_inGasLast.
function radiobutton_inGasLast_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_inGasLast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_inGasLast



function edit_air_in_pr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_air_in_pr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_air_in_pr as text
%        str2double(get(hObject,'String')) returns contents of edit_air_in_pr as a double


% --- Executes during object creation, after setting all properties.
function edit_air_in_pr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_air_in_pr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_roller_properties_Callback(hObject, eventdata, handles)
% hObject    handle to edit_roller_properties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_roller_properties as text
%        str2double(get(hObject,'String')) returns contents of edit_roller_properties as a double


% --- Executes during object creation, after setting all properties.
function edit_roller_properties_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_roller_properties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_layerNo_Callback(hObject, eventdata, handles)
% hObject    handle to edit_layerNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_layerNo as text
%        str2double(get(hObject,'String')) returns contents of edit_layerNo as a double


% --- Executes during object creation, after setting all properties.
function edit_layerNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_layerNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
