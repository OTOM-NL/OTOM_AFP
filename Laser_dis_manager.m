function varargout = Laser_dis_manager(varargin)
% LASER_DIS_MANAGER MATLAB code for Laser_dis_manager.fig
%      LASER_DIS_MANAGER, by itself, creates a new LASER_DIS_MANAGER or raises the existing
%      singleton*.
%
%      H = LASER_DIS_MANAGER returns the handle to a new LASER_DIS_MANAGER or the handle to
%      the existing singleton*.
%
%      LASER_DIS_MANAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LASER_DIS_MANAGER.M with the given input arguments.
%
%      LASER_DIS_MANAGER('Property','Value',...) creates a new LASER_DIS_MANAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Laser_dis_manager_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Laser_dis_manager_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Laser_dis_manager

% Last Modified by GUIDE v2.5 13-Apr-2019 18:11:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Laser_dis_manager_OpeningFcn, ...
                   'gui_OutputFcn',  @Laser_dis_manager_OutputFcn, ...
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

% --- Executes just before Laser_dis_manager is made visible.
function Laser_dis_manager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Laser_dis_manager (see VARARGIN)

% Choose default command line output for Laser_dis_manager
handles.output = hObject;


javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

% Rx,Ry,Rz,Ax,Ay,nx,ny,L_xyz0,ID,Gauss_Par_X,Gauss_Par_Y

% handles.Rx=varargin{1};
% handles.Ry=varargin{2};
% handles.Rz=varargin{3};
handles.Divergence_factor=0;



handles.Ax=varargin{1};
set(handles.edit4,'String',num2str(handles.Ax));

handles.Ay=varargin{2};
set(handles.edit3,'String',num2str(handles.Ay));

handles.nx=varargin{3};
set(handles.edit6,'String',num2str(handles.nx));

handles.ny=varargin{4};
set(handles.edit5,'String',num2str(handles.ny));

% handles.L_xyz0=varargin{5};
handles.ID=varargin{5};

set(handles.edit1,'String',num2str(handles.ID(2)));

set(handles.edit2,'String',num2str(handles.ID(3)));

% handles.Gauss_Par_X=varargin{10};
% handles.Gauss_Par_Y=varargin{11};

% set(handles.edit7,'String',num2str(handles.ID(4)));
% 
% set(handles.edit8,'String',num2str(handles.ID(5)));


% set(handles.edit1,'Visible','off');
% set(handles.edit2,'Visible','off');
set(handles.edit7,'Visible','off');
set(handles.edit8,'Visible','off');
set(handles.edit9,'Visible','off');
set(handles.edit10,'Visible','off');
% 
% set(handles.text1,'Visible','off');
% set(handles.text2,'Visible','off');
set(handles.text7,'Visible','off');
set(handles.text8,'Visible','off');
set(handles.text9,'Visible','off');
set(handles.text10,'Visible','off');


% Update handles structure
guidata(hObject, handles);

% set(gcf,'color','w');

% This sets up the initial plot - only do when we are invisible
% so window can get raised using Laser_dis_manager.
% if strcmp(get(hObject,'Visible'),'off')
%     plot(rand(5));
% end

% UIWAIT makes Laser_dis_manager wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Laser_dis_manager_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% varargout{1}=
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

ID=zeros(1,7);

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        ID(1)=0;
    case 2
%         plot(sin(1:0.01:25.99));
          ID(1)=2;
    case 3
        ID(1)=3;
     case 4
          ID(1)=1;
%         plot(membrane);
%     case 5
%         surf(peaks);
end

ID(2)=str2double(get(handles.edit1,'String'));
ID(3)=str2double(get(handles.edit2,'String'));



Ax=str2double(get(handles.edit4,'String'));
Ay=str2double(get(handles.edit3,'String'));

nx=str2double(get(handles.edit6,'String'));
ny=str2double(get(handles.edit5,'String'));

ID(4)=str2double(get(handles.edit7,'String'));
ID(5)=str2double(get(handles.edit8,'String'));

ID(6)=str2double(get(handles.edit9,'String'));
ID(7)=str2double(get(handles.edit10,'String'));


laser_dis_management (Ax,Ay,nx,ny,...
                    ID);
                
  handles.Ax=Ax;
  handles.Ay=Ay;
   handles.nx=nx;
  handles.ny=ny;
  
  handles.Gauss_Par_X=ID(2);
   handles.Gauss_Par_Y=ID(3);
   
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
set(handles.edit7,'Visible','off');
set(handles.edit8,'Visible','off');
set(handles.edit9,'Visible','off');
set(handles.edit10,'Visible','off');
% 
% set(handles.text1,'Visible','off');
% set(handles.text2,'Visible','off');
set(handles.text7,'Visible','off');
set(handles.text8,'Visible','off');
set(handles.text9,'Visible','off');
set(handles.text10,'Visible','off');

               
           end
           
    case 2
 resp=get(handles.edit1,'Visible');
       
           if resp(2)=='n'
           
% set(handles.edit1,'Visible','off');
% set(handles.edit2,'Visible','off');
set(handles.edit7,'Visible','off');
set(handles.edit8,'Visible','off');
set(handles.edit9,'Visible','off');
set(handles.edit10,'Visible','off');
% 
% set(handles.text1,'Visible','off');
% set(handles.text2,'Visible','off');
set(handles.text7,'Visible','off');
set(handles.text8,'Visible','off');
set(handles.text9,'Visible','off');
set(handles.text10,'Visible','off');

               
           end
    case 3
        
%     set(handles.edit1,'Visible','on');
% set(handles.edit2,'Visible','on');
set(handles.edit7,'Visible','on');
set(handles.edit8,'Visible','on');
set(handles.edit9,'Visible','on');
set(handles.edit10,'Visible','on');
% 
% set(handles.text1,'Visible','on');
% set(handles.text2,'Visible','on');
set(handles.text7,'Visible','on');
set(handles.text8,'Visible','on');
set(handles.text9,'Visible','on');
set(handles.text10,'Visible','on');
   case 4
        
%     set(handles.edit1,'Visible','on');
% set(handles.edit2,'Visible','on');
set(handles.edit7,'Visible','on');
set(handles.edit8,'Visible','on');
set(handles.edit9,'Visible','on');
set(handles.edit10,'Visible','on');
% 
% set(handles.text1,'Visible','on');
% set(handles.text2,'Visible','on');
set(handles.text7,'Visible','on');
set(handles.text8,'Visible','on');
set(handles.text9,'Visible','on');
set(handles.text10,'Visible','on');

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

set(hObject, 'String', {'Uniform', 'Gaussian','Gaussian-4D','Linear'});



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



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
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
% get(gca,'Xdata')
% 
% ch=get(gca,'ch');
ch=get(gca,'ch');
fig=figure('Name','3D Laser Source Window','NumberTitle','off');


for kk=1:length(ch)
    Str=get(ch(kk),'Type');
    if Str(1:4)=='line'
        
    else
        CD=get(ch(kk),'CData');
        
    end
    x=get(ch(kk),'XData');
    y=get(ch(kk),'YData');
    z=get(ch(kk),'ZData');
    
    set(fig,'color','w');
    
    
    if Str(1:4)=='line'
        
        plot3(x,y,z,'c--');
        hold on;
        
    else
        
        
        if iscell(x)
            for ii=1:length(x)
                surf(x{ii},y{ii},CD{ii},'LineStyle','None');
                hold on;
            end
        else
            surf(x,y,CD,'LineStyle','None');
        end
        
        
    end
    
end
axis equal;
colorbar;



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_pattern_maker.
function pushbutton_pattern_maker_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pattern_maker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% pattern from  a side







ch=get(gca,'ch');

CD=get(ch,'CData');
x=get(ch,'XData');
y=get(ch,'YData');
z=get(ch,'ZData');





x_max=max(max(x));
y_max=max(max(y));
% [z_max,index]=max(z(:,:));


% which quarter has the maximum value 
Q3=CD(1,1);
Q4=CD(1,end);
Q2=CD(end,1);
Q1=CD(end,end);

[Q_max, Q_max_index]=max([Q1 Q2 Q3 Q4]);


switch Q_max_index
    case 1
        x=x-x_max;
        y=y-y_max;
        
    case 2
        x=x+x_max;
        y=y-y_max;
        
    case 3
        x=x+x_max;
        y=y+y_max;
        
    case 4
        x=x-x_max;
        y=y+y_max;
end


% Total energy should be devided by 4
CD=CD/4;

% making pattern based on the refference

cla;
mesh(x,y,CD);
hold on;
mesh(-x,y,CD);

hold on;
mesh(-x,-y,CD);
hold on;
mesh(x,-y,CD);


% --- Executes on button press in Load_intensity.
function Load_intensity_Callback(hObject, eventdata, handles)
% hObject    handle to Load_intensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fileName,PathName] = uigetfile('*.txt','Select the intensity text file');
fid1 = fopen(strcat(PathName,fileName));

% fid1 = fopen('Geo_Parmeter.txt');
out = textscan(fid1,'%s','delimiter',',');

% [X,Y]=meshgrid(linspace (-handles.Ax,handles.Ax,handles.nx),linspace (-handles.Ax,handles.Ax,handles.nx) );



% for ii=1:length(out{1,1})/2
%     handles.Geometrical_parameters{ii}=out{1,1}{2*ii,1};
% end
% celldisp(handles.Geometrical_parameters)
guidata(hObject, handles);
fclose(fid1);


% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% Data to be written
% Fig_handles.Laser_direction
%  Fig_handles.Laser_Head_L_xyz0
% Fig_handles.nip_point_M=nip_point_M;


fid10 = fopen('.\Process_Parameter.txt','r');
  out = textscan(fid10,'%s','delimiter',',');
  
  fclose (fid10);
  
out_copy=out;
   

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        No=0;
    case 2
%         plot(sin(1:0.01:25.99));
         No=2;
    case 3
        No=3;

          case 4
        No=1;
        
end

% find index to replace
     ID = contains(out{1},'ID');
     
     No=num2str(No); %num2str(get(handles.popupmenu1,'value'));
   

      
      % remove first index to adjust writing new data
      
      out{1}(1)=[];
      
      
%        Fig_handles = guidata(gcbo);
      
       New_laser_ID=[No,' ',get(handles.edit1,'String'),' ',get(handles.edit2,'String'),' ',get(handles.edit7,'String'),' ',...
          get(handles.edit8,'String'),' ',get(handles.edit9,'String'),' ',get(handles.edit10,'String')];
      
       out{1}{ID}= New_laser_ID;
%          out{1}{L_xyz0_ind}= num2str(Fig_handles.Laser_Head_L_xyz0);
%               out{1}{tv_ind}=num2str(Fig_handles.nip_point_M);
              
%                 rotation_deg=get(sld,'value');
%         set(text_Rot,'String',sprintf('Laser Head Rotation=%d°',rotation_deg));
                     
%                out{1}{Laser_head_Rot}=num2str(rotation_deg);
              
              fid11 = fopen('.\Process_Parameter.txt','w');
              % >> regexprep  % to replace in txt file
              
              for ii=1:2:length(out{1})  % number of lines
      
         fprintf(fid11,'%s\n', strcat(out_copy{1}{ii},', ' ,out{1}{ii}) );

    
              end
              
              %%
           
              % Save the Calculated laser dovergence Factor
fid13 = fopen('.\Supp_files\Laser_characteristic.txt','w');


fprintf(fid13,'Laser_Divergence(Degree) , %f \r\n',handles.Divergence_factor);

fclose(fid13);
              
              


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_Div_Angle_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Div_Angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Div_Angle as text
%        str2double(get(hObject,'String')) returns contents of edit_Div_Angle as a double


% --- Executes during object creation, after setting all properties.
function edit_Div_Angle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Div_Angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_Divergence.
function pushbutton_Divergence_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Divergence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



axes(handles.axes1);
cla;

% ID=zeros(1,7);

% popup_sel_index = get(handles.popupmenu1, 'Value');
% switch popup_sel_index
%     case 1
%         ID(1)=0;
%     case 2
% %         plot(sin(1:0.01:25.99));
%           ID(1)=2;
%     case 3
%         ID(1)=3;
% %     case 4
% %         plot(membrane);
% %     case 5
% %         surf(peaks);
% end

% ID(2)=str2double(get(handles.edit1,'String'));
% ID(3)=str2double(get(handles.edit2,'String'));



Max_Div_Angle=str2double(get(handles.edit_Div_Angle,'String'));

Ax=str2double(get(handles.edit4,'String'));
Ay=str2double(get(handles.edit3,'String'));

nx=str2double(get(handles.edit6,'String'));
ny=str2double(get(handles.edit5,'String'));

% ID(4)=str2double(get(handles.edit7,'String'));
% ID(5)=str2double(get(handles.edit8,'String'));
% 
% ID(6)=str2double(get(handles.edit9,'String'));
% ID(7)=str2double(get(handles.edit10,'String'));


handles.Divergence_factor=laser_Divergence_management  (Ax,Ay,nx,ny,...
                    Max_Div_Angle);
                
  handles.Ax=Ax;
  handles.Ay=Ay;
   handles.nx=nx;
  handles.ny=ny;
  
%   handles.Gauss_Par_X=ID(2);
%    handles.Gauss_Par_Y=ID(3);
   
%      handles.Gauss_Par_X2=ID(4);
%    handles.Gauss_Par_Y2=ID(5);
   %      handles.Gauss_Par_X2=ID(6);
%    handles.Gauss_Par_Y2=ID(7);
   
   guidata(hObject, handles);
