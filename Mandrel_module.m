function varargout = Mandrel_module(varargin)
% Mandrel_module MATLAB code for Mandrel_module.fig
%      Mandrel_module, by itself, creates a new Mandrel_module or raises the existing
%      singleton*.
%
%      H = Mandrel_module returns the handle to a new Mandrel_module or the handle to
%      the existing singleton*.
%
%      Mandrel_module('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Mandrel_module.M with the given input arguments.
%
%      Mandrel_module('Property','Value',...) creates a new Mandrel_module or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Mandrel_module_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Mandrel_module_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Mandrel_module

% Last Modified by GUIDE v2.5 29-Jan-2018 11:35:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Mandrel_module_OpeningFcn, ...
                   'gui_OutputFcn',  @Mandrel_module_OutputFcn, ...
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

% --- Executes just before Mandrel_module is made visible.
function Mandrel_module_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Mandrel_module (see VARARGIN)

% Choose default command line output for Mandrel_module
handles.output = hObject;


javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

% Rx,Ry,Rz,Ax,Ay,nx,ny,L_xyz0,ID,Gauss_Par_X,Gauss_Par_Y

handles.R_cyl=varargin{1};
handles.z_cyl_end=varargin{2};

set(handles.edit1,'String',num2str(handles.R_cyl));

set(handles.edit2,'String',num2str(handles.z_cyl_end));


% handles.Rz=varargin{3};
% 
% handles.Ax=varargin{4};
% set(handles.edit4,'String',num2str(handles.Ax));
% 
% handles.Ay=varargin{5};
% set(handles.edit3,'String',num2str(handles.Ay));
% 
% handles.nx=varargin{6};
% set(handles.edit6,'String',num2str(handles.nx));
% 
% handles.ny=varargin{7};
% set(handles.edit5,'String',num2str(handles.ny));
% 
% handles.L_xyz0=varargin{8};
% handles.ID=varargin{9};





% set(handles.edit1,'Visible','off');
% set(handles.edit2,'Visible','off');
% set(handles.edit7,'Visible','off');
% set(handles.edit8,'Visible','off');
% set(handles.edit9,'Visible','off');
% set(handles.edit10,'Visible','off');
% 
% set(handles.text1,'Visible','off');
% set(handles.text2,'Visible','off');
% set(handles.text7,'Visible','off');
% set(handles.text8,'Visible','off');
% set(handles.text9,'Visible','off');
% set(handles.text10,'Visible','off');

set(gcf,'color','w');
% Update handles structure
guidata(hObject, handles);

% set(gcf,'color','w');

% This sets up the initial plot - only do when we are invisible
% so window can get raised using Mandrel_module.
% if strcmp(get(hObject,'Visible'),'off')
%     plot(rand(5));
% end

% UIWAIT makes Mandrel_module wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Mandrel_module_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% varargout{1}=handles.Ax;
% varargout{2}=handles.Ay;
% varargout{3}=handles.Gauss_Par_X;
% varargout{4}=handles.Gauss_Par_X;

guidata(hObject, handles);

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;

ID=zeros(1,3);

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
       Type=1;
    case 2
   Type=2;
    case 3
          Type=3;
    case 4
         Type=4;
    case 5
          Type=5;
          
              case 6
          Type=6;
end

R_cyl=str2double(get(handles.edit1,'String'));
z_cyl_end=str2double(get(handles.edit2,'String'));



c1=str2double(get(handles.edit3,'String'));
c2=str2double(get(handles.edit4,'String'));

a1=str2double(get(handles.edit5,'String'));
b1=str2double(get(handles.edit6,'String'));


finer=str2double(get(handles.edit7,'String'));
% ID(5)=str2double(get(handles.edit8,'String'));
% 
% ID(6)=str2double(get(handles.edit9,'String'));
% ID(7)=str2double(get(handles.edit10,'String'));

transparency=get(handles.slider3,'Value');
opaqueness=1-transparency;


       RGB=str2num(get(handles.edit11,'String'));         
                handles.Mandrel_plot=Plot_Mandrels(R_cyl,z_cyl_end,Type,finer,c1,c2,a1,b1,RGB,opaqueness);
                
%   handles.Ax=Ax;
%   handles.Ay=Ay;
%    handles.nx=nx;
%   handles.ny=ny;
  
%   handles.Gauss_Par_X=ID(2);
%    handles.Gauss_Par_Y=ID(3);
   
%      handles.Gauss_Par_X2=ID(4);
%    handles.Gauss_Par_Y2=ID(5);
   %      handles.Gauss_Par_X2=ID(6);
%    handles.Gauss_Par_Y2=ID(7);
   
   guidata(hObject, handles);
  

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
           resp=get(handles.edit1,'Visible');
       
           if resp(2)=='n'
           
% set(handles.edit1,'Visible','off');
% set(handles.edit2,'Visible','off');
% set(handles.edit7,'Visible','off');
% set(handles.edit8,'Visible','off');
% set(handles.edit9,'Visible','off');
% set(handles.edit10,'Visible','off');
% 
% set(handles.text1,'Visible','off');
% set(handles.text2,'Visible','off');
% set(handles.text7,'Visible','off');
% set(handles.text8,'Visible','off');
% set(handles.text9,'Visible','off');
% set(handles.text10,'Visible','off');

               
           end
           
    case 2
 resp=get(handles.edit1,'Visible');
       
           if resp(2)=='n'
           
% set(handles.edit1,'Visible','off');
% set(handles.edit2,'Visible','off');
% set(handles.edit7,'Visible','off');
% set(handles.edit8,'Visible','off');
% set(handles.edit9,'Visible','off');
% set(handles.edit10,'Visible','off');
% 
% set(handles.text1,'Visible','off');
% set(handles.text2,'Visible','off');
% set(handles.text7,'Visible','off');
% set(handles.text8,'Visible','off');
% set(handles.text9,'Visible','off');
% set(handles.text10,'Visible','off');

               
           end
    case 3
        
%     set(handles.edit1,'Visible','on');
% set(handles.edit2,'Visible','on');
% set(handles.edit7,'Visible','on');
% set(handles.edit8,'Visible','on');
% set(handles.edit9,'Visible','on');
% set(handles.edit10,'Visible','on');
% 
% set(handles.text1,'Visible','on');
% set(handles.text2,'Visible','on');
% set(handles.text7,'Visible','on');
% set(handles.text8,'Visible','on');
% set(handles.text9,'Visible','on');
% set(handles.text10,'Visible','on');


end


% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
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

set(hObject, 'String', {'Cylinder', 'Cylinder+Spherical Domes','Cylinder+heighted Domes','Ellipsoidical Cylinder','Ellipsoidical Cylinder+Domes','Domes'});



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% function edit8_Callback(hObject, eventdata, handles)
% % hObject    handle to edit8 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edit8 as text
% %        str2double(get(hObject,'String')) returns contents of edit8 as a double
% 
% 
% % --- Executes during object creation, after setting all properties.
% function edit8_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit8 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% get(gca,'Xdata')
% 
% ch=get(gca,'ch');

f=handles.Mandrel_plot;
% ch=get(handles.Mandrel_plot,'ch');

% CD=get(ch,'CData');
fig=figure('Name','3D Mandrel','NumberTitle','off');
set(fig,'color','w');
hold on;

for ii=1:length(f)

x=get(f(ii),'XData');
y=get(f(ii),'YData');
z=get(f(ii),'ZData');
alpha=get(f(ii),'facealpha');

if isvector(x)
    C=z*0;
    fill3(x,y,z,C,'FaceAlpha',0.75);
else

surf(x,y,z,'LineStyle','none','FaceAlpha',alpha);
end

end
axis tight;
axis equal;
% colorbar;
view([30 20]);
colormap summer;
% camlight;
h = camlight('left');
h2 = camlight('right');
for i = 1:90;
   camorbit(2,0);
%    camlight(h,'left');
   lightangle(h,-i,-i);
   pause(.02);
end



% function edit9_Callback(hObject, eventdata, handles)
% % hObject    handle to edit9 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
% function edit9_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit9 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% 
% function edit10_Callback(hObject, eventdata, handles)
% % hObject    handle to edit10 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edit10 as text
% %        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
% function edit10_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit10 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
