function varargout = LATW_amblifibre(varargin)
%LATW_AMBLIFIBRE MATLAB code file for LATW_amblifibre.fig
%      LATW_AMBLIFIBRE, by itself, creates a new LATW_AMBLIFIBRE or raises the existing
%      singleton*.
%
%      H = LATW_AMBLIFIBRE returns the handle to a new LATW_AMBLIFIBRE or the handle to
%      the existing singleton*.
%
%      LATW_AMBLIFIBRE('Property','Value',...) creates a new LATW_AMBLIFIBRE using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to LATW_amblifibre_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      LATW_AMBLIFIBRE('CALLBACK') and LATW_AMBLIFIBRE('CALLBACK',hObject,...) call the
%      local function named CALLBACK in LATW_AMBLIFIBRE.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LATW_amblifibre

% Last Modified by GUIDE v2.5 30-Dec-2019 16:28:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @LATW_amblifibre_OpeningFcn, ...
    'gui_OutputFcn',  @LATW_amblifibre_OutputFcn, ...
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


% --- Executes just before LATW_amblifibre is made visible.
function LATW_amblifibre_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for LATW_amblifibre
% V.1.11 >> Main changes: Fix bugs , In-out-surface detection of points, tranparency

handles.output = hObject;
set(handles.figure1, 'Name', 'OTOM V.1.20');
% get(handles);

% handles.interpreter='tex';
% Update handles structure
% handles.handlesArray = [findobj('Style','edit') ,findobj('Style', 'text')];

% handles.Optimization(1) = uicontrol(gcf,'Style','edit','Horizontalalignment','right','visible','on','BackgroundColor','white','enable','on','Position',[400,400,100,20]);
% handles.Optimization(2) = uicontrol(gcf,'Style','Text','String','What is the desired Temperature of Substrate nip-point?','Horizontalalignment','Center','visible','on','enable','on','Position',[350,420,200,40]);
% handles.Optimization(3) = uicontrol(gcf,'Style','edit','Horizontalalignment','right','visible','on','BackgroundColor','white','enable','on','Position',[400,340,100,20]);
% handles.Optimization(4) = uicontrol(gcf,'Style','Text','String','What is the desired Temperature of Tape nip-point?','Horizontalalignment','Center','visible','on','enable','on','Position',[350,360,200,40]);

% set(handles.Optimization,'Visible','Off');


set(handles.uipanel7,'Visible','Off');

set(handles.OP_edit1,'Visible','Off');
set(handles.OP_text1,'Visible','Off');
set(handles.OP_edit2,'Visible','Off');
set(handles.OP_text2,'Visible','Off');

% gcf.ToolBar = 'none';
% fig.MenuBar = 'none';

handles.fit_Func=[];
handles.checkbox1=[];
handles.checkbox_transient_code=0;
handles.Input_reader_mode=0;
handles.checkbox_Inline_code=0;
handles.BRDF_mode=0;

handles.UOTmode=0;





handles.ADS_Vars=[];
handles.TwinCat_file_name=[];

handles.manufacturing_type=cell(1,7);

handles.step_UOT_steadystate=[];
handles.nonlinearMaterials_Tape_Sub=[];



url = 'https://www.otomcomposite.eu/wp-content/uploads/OTOM_Main_window.png';
rgb = webread(url);
imshow(rgb);
axis image
% C = imread('UT_Logo_2400_White_EN.png');
% image(C);
axis off;


% C = imread('UT_Logo_2400_White_EN.png');
% imagesc(C);
% axis image;
% axis off;

set(gcf,'color','w');


javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

handles.Node_z_3D_thermal=0;

guidata(hObject, handles);




 




% UIWAIT makes LATW_amblifibre wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LATW_amblifibre_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
options.Interpreter = 'tex';
varargout{1} = handles.output;


% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fclose all;


% Turn off warning
[a, MSGID] = lastwarn();
warning('off', MSGID);
%% return struc data into numeric and Run

%% Geometrical Parameters
th_y=str2double(handles.Geometrical_parameters{1});
W_tape = str2double(handles.Geometrical_parameters{2}); % width of the tape
thick_T=str2double(handles.Geometrical_parameters{3}); % thickness of the tape
R_tape=str2double(handles.Geometrical_parameters{4});   % without thickness, it will be added in general 3D tape cylinder
L_flat=str2double(handles.Geometrical_parameters{5});
deg_tape=str2double(handles.Geometrical_parameters{6});
sui=str2double(handles.Geometrical_parameters{7}) ; %pi/8;  % tape_direction in radian, effect of rotation and sui
L_prim=str2double(handles.Geometrical_parameters{8});   % tape long
% th_1=str2double(handles.Geometrical_parameters{9}); % starting angle in degree
w=str2double(handles.Geometrical_parameters{9}); % width of the substrate tape
thick_sub=str2double(handles.Geometrical_parameters{10}); % thiness of substrate
R_cyl= str2num(handles.Geometrical_parameters{11}); %80; %input('The radius of cylinder =');
z_cyl_end= str2double(handles.Geometrical_parameters{12}); %input('Enetr the end point of cylinder =');
Roller_Pos_TV=str2num(handles.Geometrical_parameters{13});  % for the tape
W_R=str2double(handles.Geometrical_parameters{14});  % width of the Roller
Rxyz=str2num(handles.Geometrical_parameters{15});
Laser_head=str2num(handles.Geometrical_parameters{16});
L_xyz0=str2num(handles.Geometrical_parameters{17});
Laser_head_Rot=str2double(handles.Geometrical_parameters{18});


%% Process Parameters
materials_Tape=str2num(handles.Process_parameters{1});
materials_sub=str2num(handles.Process_parameters{2});
Velocity=str2double(handles.Process_parameters{3});
Total_energy=str2double(handles.Process_parameters{4});
ID=str2num(handles.Process_parameters{5});   % laser distribution pattern
absorbtion_Refractive_ind=str2num(handles.Process_parameters{6});
Temp_Right_sub=str2num(handles.Process_parameters{7}); % Includes also mandrel Temperature
Temp_Right_T_Roller=str2num(handles.Process_parameters{8});
h_conv_Sub=str2double(handles.Process_parameters{9});
h_conv_T=str2num(handles.Process_parameters{10});
Roller_Force=str2double(handles.Process_parameters{11});

assignin('base','Consolidation_Force',Roller_Force);

if Roller_Force ~=0
    
    if isempty(handles.fit_Func)
        
        uiwait(warndlg('No Force-displacement data, Please load the data.'));
        
        %          javaFrame    = get(warn_dlg,'JavaFrame');
        % iconFilePath = 'OTOM-icon.png';
        % javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
        
        H_indentation=0;
        %         Roller_def_Callback(hObject, eventdata, handles);
    else
        
        fitresult= handles.fit_Func;
        H_indentation=fitresult(Roller_Force);
        figure(50);
        plot(Roller_Force,H_indentation,'ko');
        text (Roller_Force,H_indentation,' Maximum normal deformation');
        javaFrame    = get(gcf,'JavaFrame');
        iconFilePath = 'OTOM-icon.png';
        javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
        
    end
    
else
    H_indentation=0;
    
end

if H_indentation > R_tape
    h=errordlg('Normal displacement is more than Roller radius!');
    javaFrame    = get(h,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    H_indentation=R_tape;
end



%% Computational parameters
step_size_angle=str2double(handles.Computational_parameters{1});  % in radian
No_dev=str2double(handles.Computational_parameters{2});
node_space=str2double(handles.Computational_parameters{3});  % change to the number of devision
Angular_space=str2double(handles.Computational_parameters{4});
L_flat_space=str2double(handles.Computational_parameters{5});  % change to the number of devision

%%
% 2D/3D thermal model


Node_z_3D_thermal=handles.Node_z_3D_thermal;

%
% Run loading

% dos('animator_loading\animator_19_frame_camtes.exe -i  &')
dos('animator_19_frame_camtes.exe -i  &');

%% run the Main Optical-thermal program



%% Temperature of the mandrel
if length (Temp_Right_sub) ==2
    T_amb_mandrel=Temp_Right_sub(2);
else
    T_amb_mandrel=30;
end
%%  Check graphical outputs
if isempty(handles.checkbox1)
    Graphic_chekbox=ones(1,9);
else
    Graphic_chekbox=[handles.checkbox1,handles.checkbox2,handles.checkbox3,handles.checkbox4,...
        handles.checkbox5,handles.checkbox6,handles.checkbox7,handles.checkbox8,...
        handles.checkbox9];
end

if handles.checkbox_transient_code
    Transient_ID=[handles.checkbox_transient_code, handles.checkbox_Live_output,...
        str2double({handles.edit_init_Temp,...
        handles.edit_inc, handles.edit_Time,handles.checkbox_Video})];
else
      Transient_ID=[0, 0,20,        1, 1,0];
end


% Should be modified later
if handles.checkbox_Inline_code
    Inline_ID=[] ;
end


% use or not use divergence of the laser
Laser_div_OnOff=get(handles.checkbox_laser_divergence,'Value');

if Laser_div_OnOff

 fid23 = fopen('.\Supp_files\Laser_characteristic.txt','r');
 out = textscan(fid23,'%s ','delimiter',',');
 fclose(fid23);
 
Divergence_factor=str2double(out{1}{2});

else
    Divergence_factor=0;
    
end



if handles.UOTmode
    
    switch  handles.UOTmode
        
        case 1
    str= handles.UOT_optical;
    load(str);
    
    
    
    N_tape = 360; % EVEN number of points, should not be changed !!
    Tape_Sp=[N_tape;W_tape;R_tape;L_flat;thick_T;deg_tape];  %Tape_Specification
    
    [counter_ray,Rot_Roller_axis,tv]=General_3D_Optical_UOT(th_y,Tape_Sp,...
        Rxyz,R_cyl,z_cyl_end,Roller_Pos_TV,ID,Laser_head,absorbtion_Refractive_ind,W_R,H_indentation,...
        Graphic_chekbox,...
        handles.BRDF_mode,jobname{:},...
        CV_mesh,LaserH_mesh,delta_new_all,pos_new_all,nip_point_M_all,Rot_Roller_axis_all,...
        handles.UOT_pathfile,handles.text_status,Divergence_factor);
    
        case 2
    
    str= handles.UOT_Thermal;
    str2= handles.UOT_PreThermal;
    load(str);
    load(str2);
    
     
     
     
    % pay aqttention that Geometrical parameters from analysed optical UOT
        % should be considered not New parameters
if isempty (handles.step_UOT_steadystate)
    
    % Transient 
    thermal_domain_generator_UOT(th_y,W_tape,thick_T,R_tape,L_flat,deg_tape,...
    sui,L_prim,w,thick_sub,No_dev,R_cyl,...
    W_R,materials_Tape,materials_sub,Total_energy,Velocity,...
       node_space,Angular_space,L_flat_space,...
    Temp_Right_sub,Temp_Right_T_Roller,h_conv_Sub,h_conv_T,H_indentation,Node_z_3D_thermal,T_amb_mandrel,....
    Graphic_chekbox,Transient_ID,handles.manufacturing_type,...
     nip_point_M_all,Rot_Roller_axis_all,...
        handles.UOT_pathfile,CV_mesh,Laser_head,ID,handles.text_status,...
        Points_in_domain_all_sub,Points_in_domain_all_Tape,Tape_points_Data)
    
else
    % Steadystate
    UOT_step_SS=str2num(handles.step_UOT_steadystate{:});
    
        thermal_domain_generator_UOT_SteadyState(th_y,W_tape,thick_T,R_tape,L_flat,deg_tape,...
    sui,L_prim,w,thick_sub,No_dev,R_cyl,...
    W_R,materials_Tape,materials_sub,Total_energy,Velocity,...
       node_space,Angular_space,L_flat_space,...
    Temp_Right_sub,Temp_Right_T_Roller,h_conv_Sub,h_conv_T,H_indentation,Node_z_3D_thermal,T_amb_mandrel,....
    Graphic_chekbox,Transient_ID,handles.manufacturing_type,...
     nip_point_M_all,Rot_Roller_axis_all,...
        handles.UOT_pathfile,CV_mesh,Laser_head,ID,handles.text_status,...
        Points_in_domain_all_sub,Points_in_domain_all_Tape,Tape_points_Data,...
       UOT_step_SS )
    
end
    
    
         case 5
    
    str= handles.UOT_Thermal;
    load(str);
    
    
    
   
     
     
    
    Pre_Thermal_UOT(W_tape,thick_T,R_tape,L_flat,deg_tape,...
    R_cyl, W_R,node_space,Angular_space,L_flat_space,...
    H_indentation,nip_point_M_all,Rot_Roller_axis_all,...
    handles.UOT_pathfile,CV_mesh,handles.text_status)
    
    
    end
    
    
else
    
    
    if handles.checkbox_transient_code
        %% Transinet Model solver
        
        %          manufacturing_type{1}=handles.edit_globalSub;
        % manufacturing_type{2}=handles.edit_globalTape;
        % manufacturing_type{3}=handles.edit_initial_sub;
        % manufacturing_type{4}=handles.edit_initial_Tape;
        % manufacturing_type{5}=handles.edit_LocationNo;
        %   manufacturing_type{6}=handles.edit_pitchAngle;
        %   manufacturing_type{7}=handles.edit_delayTime;
        %   manufacturing_type{8}=handles.radiobutton_Wdir;
        
        thermal_domain_generator_Tr(th_y,W_tape,thick_T,R_tape,L_flat,deg_tape,...
            sui,L_prim,step_size_angle,w,thick_sub,No_dev,R_cyl,z_cyl_end,...
            Roller_Pos_TV,W_R,materials_Tape,materials_sub,Velocity,Total_energy,...
            Rxyz,ID,Laser_head,L_xyz0,absorbtion_Refractive_ind,...
            node_space,Angular_space,L_flat_space,...
            Temp_Right_sub,Temp_Right_T_Roller,h_conv_Sub,h_conv_T,H_indentation,Node_z_3D_thermal,T_amb_mandrel,...
            Graphic_chekbox,handles.Input_reader_mode,Transient_ID,handles.manufacturing_type,Laser_head_Rot,...
            handles.BRDF_mode,handles.text_status,Divergence_factor,...
            handles.nonlinearMaterials_Tape_Sub);
        
        
    elseif handles.checkbox_Inline_code
        
        if  handles.TwinCat_file_name
            
            Inline_thermal_domain_generator_Tr(th_y,W_tape,thick_T,R_tape,L_flat,deg_tape,...
                sui,L_prim,step_size_angle,w,thick_sub,No_dev,R_cyl,z_cyl_end,...
                Roller_Pos_TV,W_R,materials_Tape,materials_sub,...
                Rxyz,ID,Laser_head,L_xyz0,absorbtion_Refractive_ind,...
                node_space,Angular_space,L_flat_space,...
                Temp_Right_sub,Temp_Right_T_Roller,h_conv_Sub,h_conv_T,H_indentation,Node_z_3D_thermal,T_amb_mandrel,...
                Graphic_chekbox,Inline_ID,handles.manufacturing_type,Laser_head_Rot,...
                handles.ADS_Vars,handles.TwinCat_file_name,handles.BRDF_mode ,Divergence_factor);
        else
            error('Program terminated >> No ADS initialization specific reason');
            %  drawnow();
            
            
            
        end
        
        
    else
        [Nip_point_Temp_T]=thermal_domain_generator_opt(th_y,W_tape,thick_T,R_tape,L_flat,deg_tape,...
            sui,L_prim,step_size_angle,w,thick_sub,No_dev,R_cyl,z_cyl_end,...
            Roller_Pos_TV,W_R,materials_Tape,materials_sub,Velocity,Total_energy,...
            Rxyz,ID,Laser_head,L_xyz0,absorbtion_Refractive_ind,...
            node_space,Angular_space,L_flat_space,...
            Temp_Right_sub,Temp_Right_T_Roller,h_conv_Sub,h_conv_T,H_indentation,Node_z_3D_thermal,T_amb_mandrel,...
            Graphic_chekbox,handles.Input_reader_mode,Laser_head_Rot,handles.BRDF_mode,Divergence_factor);
        
        
        
    end
end

%Kill loading show success!
dos('taskkill /im animator_19_frame_camtes.exe');

[cdata,map] = imread('complete.png');

h2=msgbox('Operation Completed!',...
    'Success','custom',cdata,map);

javaFrame    = get(h2,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));







% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



close(figure(1));
close(figure(2));
close(figure(3));

close(figure(21));
close(figure(100));

% 3D thermal model figures

close(figure(36));
close(figure(46));
close(figure(56));

% Fitting figures
close(figure(61));
close(figure(62));
close(figure(63));
close(figure(64));
fclose all;

dos('taskkill /im animator_19_frame_camtes.exe');





% [cdata,map] = imread('under_construction.png');
%
% h2=msgbox('Under Development!',...
%     'Notice','custom',cdata,map);
% set(h2,'color','w');





% --------------------------------------------------------------------
function PostP_Callback(hObject, eventdata, handles)
% hObject    handle to PostP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function Optimization_Callback(hObject, eventdata, handles)
% hObject    handle to Optimization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% hideButtonGroup(thermalFrame); hideButtonGroup(thermalInput);
% handlesArray = [gca];

% % Set them all disabled.

% set(handles.handlesArray, 'Visible', 'off');
% set(handles.Optimization,'Visible','On');

s=get(handles.uipanel7,'Visible');
if s(2)=='f'
    set(handles.uipanel7,'Visible','On');
    
    set(handles.OP_edit1,'Visible','On');
    set(handles.OP_text1,'Visible','On');
    set(handles.OP_edit2,'Visible','On');
    set(handles.OP_text2,'Visible','On');
    
    
    
s=cell(1,8);
s{1,1}='Winding Angle(\theta_y)';
s{1,2}='Velocity';
s{1,3}=    'Total Laser Power ';
s{1,4}='Laser-Location (x,y,z)';
s{1,5}=    'Laser Direction (Rx,Ry,Rz)';  
s{1,6}='Laser Head Size';
s{1,7}=    'Laser-ID Parameters';  
s{1,8}='Mandrel Radius';

set(handles.OP_popupmenu1,'String',s) ;
    % set(handles.Optimization, 'Visible', 'On');
        handles.UOTmode=0;
    handles.UOT_Thermal=0;
    handles.UOT_PreThermal=0;
 
    
    set(handles.OP_Run,'String','Run');
    
else
    set(handles.uipanel7,'Visible','Off');
    
    set(handles.OP_edit1,'Visible','Off');
    set(handles.OP_text1,'Visible','Off');
    set(handles.OP_edit2,'Visible','Off');
    set(handles.OP_text2,'Visible','Off');
end
guidata(hObject, handles);



% --------------------------------------------------------------------
function init_Callback(hObject, eventdata, handles)
% hObject    handle to init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handlesArray = [findobj('Style','edit') ,findobj('Style', 'text')];

% % Set them all disabled.

% set(handles.handlesArray , 'Visible', 'on');

% set(handles.uipanel7,'Visible','Off');
%
% set(handles.OP_edit1,'Visible','Off');
% set(handles.OP_text1,'Visible','Off');
% set(handles.OP_edit2,'Visible','Off');
% set(handles.OP_text2,'Visible','Off');

% set(handles.Optimization, 'Visible', 'off');



function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



prompt = {'\theta_y:','Width-Tape','Thick-Tape','Radius-Tape-Roller','Length-flat-Tape','deg-tangent-Tape','\Phi-Sub','Length-Sub',...
    'Half width-Sub','Thick-Sub','Radius-Mandrel, c1, c2','Length-Mandrel / Plate ',...
    'Roller-Pos-TV','Width-Roller','Rxyz[Rx Ry Rz] (Rays orientation)','Laser-head [Ax,Ay,nx,ny] (semi-long,semi-width,No.x, No.y)',...
    'L-xyz0 [x0 y0 z0] (Laser position)', 'Laser Head Rotation (Deg)'}; % 18 numbers

dlg_title = 'Geometrical Parameters';
num_lines = 1;

fid10 = fopen('.\Geo_Parameter.txt');
if fid10 ~=-1
    out = textscan(fid10,'%s','delimiter',',');
    defaultans = {out{1}{2},out{1}{4},out{1}{6},out{1}{8},out{1}{10}...
        out{1}{12},out{1}{14},out{1}{16},out{1}{18},...
        out{1}{20},out{1}{22},out{1}{24},out{1}{26},...
        out{1}{28},out{1}{30},out{1}{32},out{1}{34},out{1}{36}};
    
else
    defaultans = {'-182','20e-3','0.15e-3','0.07','0.25','40','0','16e-2',...
        '4e-2','0.15e-3','30e-2','1',...
        '0.136000 0.100000 0.120000   ','5e-2','-0.811050 -0.584280 -0.028569  ',' 0.014000 0.031500 30 40  ','0.081866 0.212740 0.253960  ', '180'};
end

options.Interpreter = 'tex';
handles.Geometrical_parameters = inputdlg(prompt,dlg_title,num_lines,defaultans,options);
% javaFrame    = get(handles.Geometrical_parameters,'JavaFrame');
% iconFilePath = 'OTOM-icon.png';
% javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));



% fclose(fid10);



if ~isempty(handles.Geometrical_parameters)
    
    th_y=str2double(handles.Geometrical_parameters{1});
    W_tape = str2double(handles.Geometrical_parameters{2}); % width of the tape
    thick_T=str2double(handles.Geometrical_parameters{3}); % thickness of the tape
    R_tape=str2double(handles.Geometrical_parameters{4});   % without thickness, it will be added in general 3D tape cylinder
    L_flat=str2double(handles.Geometrical_parameters{5});
    deg_tape=str2double(handles.Geometrical_parameters{6});
    sui=str2double(handles.Geometrical_parameters{7}) ; %pi/8;  % tape_direction in radian, effect of rotation and sui
    L_prim=str2double(handles.Geometrical_parameters{8});   % tape long
    % th_1=str2double(handles.Geometrical_parameters{9}); % starting angle in degree
    % w=str2double(handles.Geometrical_parameters{10}); % width of the substrate tape
    % z1=str2double(handles.Geometrical_parameters{11}); % position of the nip-point z
    % R_cyl= str2num(handles.Geometrical_parameters{12}); %80; %input('The radius of cylinder =');
    % z_cyl_end= str2double(handles.Geometrical_parameters{13}); %input('Enetr the end point of cylinder =');
    % Roller_Pos_TV=str2num(handles.Geometrical_parameters{14});  % for the tape
    % W_R=str2double(handles.Geometrical_parameters{15});  % width of the Roller
    % Rxyz=str2num(handles.Geometrical_parameters{16});
    % Laser_head=str2num(handles.Geometrical_parameters{17});
    % L_xyz0=str2num(handles.Geometrical_parameters{18});
    % th_1=str2double(handles.Geometrical_parameters{9}); % starting angle in degree
    w=str2double(handles.Geometrical_parameters{9}); % width of the substrate tape
    thick_sub=str2double(handles.Geometrical_parameters{10}); % position of the nip-point z
    R_cyl= str2num(handles.Geometrical_parameters{11}); %80; %input('The radius of cylinder =');
    z_cyl_end= str2double(handles.Geometrical_parameters{12}); %input('Enetr the end point of cylinder =');
    Roller_Pos_TV=str2num(handles.Geometrical_parameters{13});  % for the tape
    W_R=str2double(handles.Geometrical_parameters{14});  % width of the Roller
    Rxyz=str2num(handles.Geometrical_parameters{15});
    Laser_head=str2num(handles.Geometrical_parameters{16});
    L_xyz0=str2num(handles.Geometrical_parameters{17});
    Laser_head_Rot=str2double(handles.Geometrical_parameters{18});
    
    fileID1 = fopen('.\Geo_Parameter.txt','w');
    
    
    if length (R_cyl)==3
        
        fprintf(fileID1,' th_y, %f  \r\n W_tape, %f  \r\n thick_T, %f  \r\n R_tape, %f  \r\n L_flat, %f  \r\n deg_tape, %f  \r\n sui, %f  \r\n L_prim, %f  \r\n w, %f  \r\n thick-sub, %f  \r\n R_cyl, %f %f %f  \r\n z_cyl_end, %f  \r\n Roller_Pos_TV, %f %f %f  \r\n W_R, %f  \r\n Rxyz, %f %f %f   \r\n Laser_head, %f %f %d %d   \r\n L_xyz0, %f %f %f  \r\n Laser_head_Rot, %f \r\n',th_y,W_tape,thick_T,R_tape,L_flat,deg_tape...
            ,sui,L_prim,w,thick_sub,R_cyl,z_cyl_end,Roller_Pos_TV,W_R,Rxyz,Laser_head,L_xyz0,Laser_head_Rot);
        
    else
        fprintf(fileID1,' th_y, %f  \r\n W_tape, %f  \r\n thick_T, %f  \r\n R_tape, %f  \r\n L_flat, %f  \r\n deg_tape, %f  \r\n sui, %f  \r\n L_prim, %f  \r\n w, %f  \r\n thick-sub, %f  \r\n R_cyl, %f   \r\n z_cyl_end, %f  \r\n Roller_Pos_TV, %f %f %f  \r\n W_R, %f  \r\n Rxyz, %f %f %f   \r\n Laser_head, %f %f %d %d   \r\n L_xyz0, %f %f %f  \r\n Laser_head_Rot, %f \r\n',th_y,W_tape,thick_T,R_tape,L_flat,deg_tape...
            ,sui,L_prim,w,thick_sub,R_cyl,z_cyl_end,Roller_Pos_TV,W_R,Rxyz,Laser_head,L_xyz0,Laser_head_Rot);
    end
    fclose(fileID1);
    
end
guidata(hObject, handles);

% % % Roller_Pos_TV should be always half of the cyliner radious !!
%%% Width of substarte always should be in the middle of the Roller



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



prompt = {'Materials-Tape [K, \rho, C_p, ky (optional),kz (optional)]','Materials-Sub [K, \rho, C_p, ky (optional),kz (optional)]','Velocity Tape feed (m/s)','Total-energy (W)','Laser-ID (distribution)',...
    'Refrective index(Sub, Tape, Roller)','In-Temp-sub, Temp Mandrel (Optional) \circ C','Incoming Temp-Tape, Temp-Roller \circ C','h-conv-Sub','h-conv-Tape with (Air, Roller)','Roller Force (N)'};
dlg_title = 'Process parameters';
num_lines = 1;


fid11 = fopen('.\Process_Parameter.txt');
if fid11 ~=-1
    out = textscan(fid11,'%s','delimiter',',');
    defaultans = {out{1}{2},out{1}{4},out{1}{6},out{1}{8},out{1}{10}...
        out{1}{12},out{1}{14},out{1}{16},out{1}{18},...
        out{1}{20},out{1}{22}};
    
else
    defaultans = {'0.72;1560;1425; 5 ; 0.7','0.72;1560;1425; 5 ; 0.7','0.1','5e2','0 ;1; 1','1.95 1.95 1.43',...
        '25 , 30','25 32','100','10 ; 100','0.1'};
end



options.Interpreter = 'tex';
handles.Process_parameters = inputdlg(prompt,dlg_title,num_lines,defaultans,options);
% javaFrame    = get(gcf,'JavaFrame');
% iconFilePath = 'OTOM-icon.png';
% javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));


if ~isempty(handles.Process_parameters )
    
    %% Process Parameters
    materials_Tape=str2num(handles.Process_parameters{1});
    materials_sub=str2num(handles.Process_parameters{2});
    Velocity=str2double(handles.Process_parameters{3});
    Total_energy=str2double(handles.Process_parameters{4});
    ID=str2num(handles.Process_parameters{5});   % laser distribution pattern
    absorbtion_Refractive_ind=str2num(handles.Process_parameters{6});
    Temp_Right_sub=str2num(handles.Process_parameters{7});
    Temp_Right_T_Roller=str2num(handles.Process_parameters{8});
    h_conv_Sub=str2double(handles.Process_parameters{9});
    h_conv_T=str2num(handles.Process_parameters{10});
    Roller_Force=str2double(handles.Process_parameters{11});
    
    fileID2 = fopen('.\Process_Parameter.txt','w');
    
    if length(ID)==3
    fprintf(fileID2,' materials_Tape, %f %f %f %f %f  \r\n materials_sub, %f %f %f %f %f \r\n Velocity, %f  \r\n Total_energy, %f  \r\n ID, %f %f %f  \r\n absorbtion_Refractive_index, %f %f %f \r\n Temp_Right_sub_Mandrel, %f %f  \r\n Temp_Right_Tape-Temp-Roller, %f %f  \r\n   h_conv_Sub, %f  \r\n h_conv_T, %f %f  \r\n  Roller Force , %f  \r\n ',materials_Tape,materials_sub,Velocity,Total_energy,...
        ID,absorbtion_Refractive_ind,Temp_Right_sub,Temp_Right_T_Roller,h_conv_Sub,h_conv_T,Roller_Force);
    
    elseif length(ID)==7
    fprintf(fileID2,' materials_Tape, %f %f %f %f %f  \r\n materials_sub, %f %f %f %f %f \r\n Velocity, %f  \r\n Total_energy, %f  \r\n ID, %f %f %f %f %f %f %f \r\n absorbtion_Refractive_index, %f %f %f  \r\n Temp_Right_sub_Mandrel, %f %f  \r\n Temp_Right_Tape-Temp-Roller, %f %f  \r\n   h_conv_Sub, %f  \r\n h_conv_T, %f %f  \r\n  Roller Force , %f  \r\n ',materials_Tape,materials_sub,Velocity,Total_energy,...
        ID,absorbtion_Refractive_ind,Temp_Right_sub,Temp_Right_T_Roller,h_conv_Sub,h_conv_T,Roller_Force);
   
    end
    
    fclose(fileID2);
    end
    
 
guidata(hObject, handles);


%  cell2mat(answer{1})



% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


prompt = {'No-length-Sub','half-side-Width-Node-Sub','Width-Nodes-Tape','Angular-Nodes-Tape','L-flat-Nodes-Tape'};
dlg_title = 'Computational parameters';
num_lines = 1;

fid12 = fopen('.\Comp_Parameter.txt');
if fid12 ~=-1
    out = textscan(fid12,'%s','delimiter',',');
    defaultans = {out{1}{2},out{1}{4},out{1}{6},out{1}{8},out{1}{10}};
    
else
    defaultans = {'40','7','13','28','46'};
end


options.Interpreter = 'tex';
handles.Computational_parameters = inputdlg(prompt,dlg_title,num_lines,defaultans,options);
% javaFrame    = get(gcf,'JavaFrame');
% iconFilePath = 'OTOM-icon.png';
% javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));




if ~isempty(handles.Computational_parameters)
    
    % Computational parameters
    No_length=str2double(handles.Computational_parameters{1});  % in radian
    No_dev=str2double(handles.Computational_parameters{2});
    Width_node=str2double(handles.Computational_parameters{3});  % change to the number of devision
    Angular_Node=str2double(handles.Computational_parameters{4});
    L_flat_Node=str2double(handles.Computational_parameters{5});  % change to the number of devision
    
    
    fileID3 = fopen('.\Comp_Parameter.txt','w');
    
    fprintf(fileID3,' No-length, %d  \r\n node_num_Z, %d  \r\n Width_node, %d  \r\n Angular_Node, %d  \r\n L_flat_space, %d  \r\n',No_length,No_dev,Width_node,Angular_Node,L_flat_Node);
    
    fclose(fileID3);
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
% C = imread('UT_Logo_2400_White_EN.png');
% image(C);
% axis off;


% --- Executes on selection change in OP_popupmenu1.
function OP_popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to OP_popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns OP_popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from OP_popupmenu1


% --- Executes during object creation, after setting all properties.
function OP_popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OP_popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OP_Range.
function OP_Range_Callback(hObject, eventdata, handles)
% hObject    handle to OP_Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





prompt = {'Lower Value','Upper Value'};
dlg_title = 'Range';
num_lines = [1 30];




if handles.UOTmode==10
    
    
     str= handles.UOT_Thermal;
%     str2= handles.UOT_PreThermal;
    load(str);
%     load(str2);
    
    
   Ax=Laser_head(1);
Ay=Laser_head(2);
    Total_energy=str2double(handles.Process_parameters{4});

set(handles.OP_popupmenu1,'String','Laser-ID 3- Power') ;


       
        defaultans = {[' 0 0 0 0 ',  num2str(-Ax*0.9) ,' ' , num2str(-Ay*0.9),' ',num2str(0.4*Total_energy)],[' 50 50 50 50 ', num2str(Ax*0.9) ,' ' , num2str(Ay*0.9),' ',num2str(2*Total_energy)]};

else 
   
    ID_parameter=get(handles.OP_popupmenu1,'Value') ;
        
switch  ID_parameter
    case   1
        %         th_y=xin;
        defaultans = {'-5','50'};
    case 2
        %         Velocity=xin;
        defaultans = {'0.05','1'};
    case 3
        %         Total_energy=xin;
        defaultans = {'1e2','9e2'};
    case 4
        %         L_xyz0=xin;
        defaultans = {'0.100000, 0.150000, 0.250000  ','0.900000, 0.950000, 2.250000 '};
    case 5
        %Laser_direction (Rx,Ry,Rz)
        %       Rxyz  =xin;
        defaultans = {'-1, -1, 1  ','1, 1, 1 '};
    case 6
        % Laser Head Size (Ax,Ay)
        %       Laser_head(1:2)  =xin;
        defaultans = {'0.100000, 0.150000  ','0.100000, 0.20000'};
        
    case 7
        %Laser-ID parameter  (Par_x, Par_Y)
        %       ID (2:3)  =xin;
        defaultans = {'10, 20  ','20, 60'};
    case 8
        % Mandrel radius R_cyl
        %      R_cyl (1)   =xin;
        defaultans = {'0.100000 ','0.600000'};
        
        
        
        
end


end
% defaultans = {'-180','-110'};

options.Interpreter = 'tex';
handles.Range = inputdlg(prompt,dlg_title,num_lines,defaultans,options);
% javaFrame    = get(gcf,'JavaFrame');
% iconFilePath = 'OTOM-icon.png';
% javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

guidata(hObject, handles);


% --- Executes on button press in OP_Run.
function OP_Run_Callback(hObject, eventdata, handles)
% hObject    handle to OP_Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% the optimization subroutine should get all inout variables and
%    setting for the optimization >> Thus we pass the inputs to this
%    function
if handles.UOTmode==10
    
%     if handles.checkbox_transient_code
%     Transient_ID=[handles.checkbox_transient_code, handles.checkbox_Live_output,...
%         str2double({handles.edit_init_Temp,...
%         handles.edit_inc, handles.edit_Time,handles.checkbox_Video})];
% else
%       Transient_ID=[0, 0,20,        1, 1,0];
% end

% if isempty(handles.checkbox1)
%     Graphic_chekbox=ones(1,9);
% else
%     Graphic_chekbox=[handles.checkbox1,handles.checkbox2,handles.checkbox3,handles.checkbox4,...
%         handles.checkbox5,handles.checkbox6,handles.checkbox7,handles.checkbox8,...
%         handles.checkbox9];
% end

%% Temperature of the mandrel
% if length (Temp_Right_sub) ==2
%     T_amb_mandrel=Temp_Right_sub(2);
% else
%     T_amb_mandrel=30;
% end

%     handles.step_UOT_optimization=Answer;
    
%     handles.UOT_pathfile=path;
    
  str= handles.UOT_Thermal;
str2= handles.UOT_PreThermal;
load(str);
load(str2);

 [x, Fval,exitFlag,Output] =    callObjConstr_4Gui_UOT(handles,...
    th_y,W_tape,thick_T,R_tape,L_flat,deg_tape,...
           R_cyl,...
            W_R,...
            H_indentation,...
                nip_point_M_all,Rot_Roller_axis_all,...
            CV_mesh,Laser_head,...
            Points_in_domain_all_sub,Points_in_domain_all_Tape,Tape_points_Data);
else

[x, Fval,exitFlag,Output] = callObjConstr_4Gui(handles);

end

%Kill loading show success!
dos('taskkill /im animator_19_frame_camtes.exe');

[cdata,map] = imread('complete.png');

h2=msgbox('Operation Completed!',...
    'Success','custom',cdata,map);

javaFrame    = get(h2,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));





% --- Executes on button press in OP_stop.
function OP_stop_Callback(hObject, eventdata, handles)
% hObject    handle to OP_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
drawnow;
% quit force;
error('Program terminated for a specific reason');
%  drawnow();
%  for i=1:2
return;

%  end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function About_Callback(hObject, eventdata, handles)
% hObject    handle to About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[cdata,map] = imread('OTOM-QRcode.png','BackgroundColor',[1 1 1]);

fig=figure('Name','otomcomposite.eu','NumberTitle','off');
fig.ToolBar = 'none';
fig.MenuBar = 'none';
title('https://www.otomcomposite.eu','Color','k');
set(fig,'color','w');

javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

imshow(cdata,map);
axis off;
title('https://www.otomcomposite.eu','Color','k','fontsize',20);

% h2=msgbox('Under Development!',...
%          'Notice','custom',cdata,map);
%      set(h2,'color','w');



% --------------------------------------------------------------------
function Documents_Callback(hObject, eventdata, handles)
% hObject    handle to Documents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

open Docs\How-to-Run-OTOM-V109.pdf;


[cdata,map] = imread('under_construction.png');

h2=msgbox('Under Development!',...
    'Notice','custom',cdata,map);
set(h2,'color','w');



% --------------------------------------------------------------------
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function Exit_Callback(hObject, eventdata, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(gcf);



% --- Executes on button press in OP_checkbox1.
function OP_checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to OP_checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OP_checkbox1


% --- Executes on button press in OP_checkbox2.
function OP_checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to OP_checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OP_checkbox2


% --- Executes when entered data in editable cell(s) in OP_uitable1.


% --------------------------------------------------------------------
function Geo_Parmeter_save_Callback(hObject, eventdata, handles)
% hObject    handle to Geo_Parmeter_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Geometrical Parameters
th_y=str2double(handles.Geometrical_parameters{1});
W_tape = str2double(handles.Geometrical_parameters{2}); % width of the tape
thick_T=str2double(handles.Geometrical_parameters{3}); % thickness of the tape
R_tape=str2double(handles.Geometrical_parameters{4});   % without thickness, it will be added in general 3D tape cylinder
L_flat=str2double(handles.Geometrical_parameters{5});
deg_tape=str2double(handles.Geometrical_parameters{6});
sui=str2double(handles.Geometrical_parameters{7}) ; %pi/8;  % tape_direction in radian, effect of rotation and sui
L_prim=str2double(handles.Geometrical_parameters{8});   % tape long
% th_1=str2double(handles.Geometrical_parameters{9}); % starting angle in degree
% w=str2double(handles.Geometrical_parameters{10}); % width of the substrate tape
% z1=str2double(handles.Geometrical_parameters{11}); % position of the nip-point z
% R_cyl= str2num(handles.Geometrical_parameters{12}); %80; %input('The radius of cylinder =');
% z_cyl_end= str2double(handles.Geometrical_parameters{13}); %input('Enetr the end point of cylinder =');
% Roller_Pos_TV=str2num(handles.Geometrical_parameters{14});  % for the tape
% W_R=str2double(handles.Geometrical_parameters{15});  % width of the Roller
% Rxyz=str2num(handles.Geometrical_parameters{16});
% Laser_head=str2num(handles.Geometrical_parameters{17});
% L_xyz0=str2num(handles.Geometrical_parameters{18});
% th_1=str2double(handles.Geometrical_parameters{9}); % starting angle in degree
w=str2double(handles.Geometrical_parameters{9}); % width of the substrate tape
thick_sub=str2double(handles.Geometrical_parameters{10}); % position of the nip-point z
R_cyl= str2num(handles.Geometrical_parameters{11}); %80; %input('The radius of cylinder =');
z_cyl_end= str2double(handles.Geometrical_parameters{12}); %input('Enetr the end point of cylinder =');
Roller_Pos_TV=str2num(handles.Geometrical_parameters{13});  % for the tape and Roller- position of the Nip-point
W_R=str2double(handles.Geometrical_parameters{14});  % width of the Roller
Rxyz=str2num(handles.Geometrical_parameters{15});
Laser_head=str2num(handles.Geometrical_parameters{16});
L_xyz0=str2num(handles.Geometrical_parameters{17});
Laser_head_Rot=str2double(handles.Geometrical_parameters{18});


[file,path] = uiputfile('.\Geo_Parameter.txt','Save file name');
% fileID3 = fopen('Comp_Parmeter.txt','w');
fileID1 = fopen(strcat(path,file),'w');
% fileID1 = fopen('Geo_Parmeter.txt','w');
%
% if length (R_cyl)==3
%
%     fprintf(fileID1,' th_y, %f  \r\n W_tape, %f  \r\n thick_T, %f  \r\n R_tape, %f  \r\n L_flat, %f  \r\n deg_tape, %f  \r\n sui, %f  \r\n L_prim, %f  \r\n w, %f  \r\n thick-sub, %f  \r\n R_cyl, %f %f %f  \r\n z_cyl_end, %f  \r\n Roller_Pos_TV, %f %f %f  \r\n W_R, %f  \r\n Rxyz, %f %f %f   \r\n Laser_head, %f %f %d %d   \r\n L_xyz0, %f %f %f  \r\n',th_y,W_tape,thick_T,R_tape,L_flat,deg_tape...
%         ,sui,L_prim,w,thick_sub,R_cyl,z_cyl_end,Roller_Pos_TV,W_R,Rxyz,Laser_head,L_xyz0);
%
% else
%     fprintf(fileID1,' th_y, %f  \r\n W_tape, %f  \r\n thick_T, %f  \r\n R_tape, %f  \r\n L_flat, %f  \r\n deg_tape, %f  \r\n sui, %f  \r\n L_prim, %f  \r\n w, %f  \r\n thick-sub, %f  \r\n R_cyl, %f  \r\n z_cyl_end, %f  \r\n Roller_Pos_TV, %f %f %f  \r\n W_R, %f  \r\n Rxyz, %f %f %f   \r\n Laser_head, %f %f %d %d   \r\n L_xyz0, %f %f %f  \r\n',th_y,W_tape,thick_T,R_tape,L_flat,deg_tape...
%         ,sui,L_prim,w,thick_sub,R_cyl,z_cyl_end,Roller_Pos_TV,W_R,Rxyz,Laser_head,L_xyz0);
% end

if length (R_cyl)==3
    
    fprintf(fileID1,' th_y, %f  \r\n W_tape, %f  \r\n thick_T, %f  \r\n R_tape, %f  \r\n L_flat, %f  \r\n deg_tape, %f  \r\n sui, %f  \r\n L_prim, %f  \r\n w, %f  \r\n thick-sub, %f  \r\n R_cyl, %f %f %f  \r\n z_cyl_end, %f  \r\n Roller_Pos_TV, %f %f %f  \r\n W_R, %f  \r\n Rxyz, %f %f %f   \r\n Laser_head, %f %f %d %d   \r\n L_xyz0, %f %f %f  \r\n  Laser_head_Rot, %f \r\n',th_y,W_tape,thick_T,R_tape,L_flat,deg_tape...
        ,sui,L_prim,w,thick_sub,R_cyl,z_cyl_end,Roller_Pos_TV,W_R,Rxyz,Laser_head,L_xyz0,Laser_head_Rot);
    
else
    fprintf(fileID1,' th_y, %f  \r\n W_tape, %f  \r\n thick_T, %f  \r\n R_tape, %f  \r\n L_flat, %f  \r\n deg_tape, %f  \r\n sui, %f  \r\n L_prim, %f  \r\n w, %f  \r\n thick-sub, %f  \r\n R_cyl, %f   \r\n z_cyl_end, %f  \r\n Roller_Pos_TV, %f %f %f  \r\n W_R, %f  \r\n Rxyz, %f %f %f   \r\n Laser_head, %f %f %d %d   \r\n L_xyz0, %f %f %f  \r\n Laser_head_Rot, %f \r\n',th_y,W_tape,thick_T,R_tape,L_flat,deg_tape...
        ,sui,L_prim,w,thick_sub,R_cyl,z_cyl_end,Roller_Pos_TV,W_R,Rxyz,Laser_head,L_xyz0,Laser_head_Rot);
end

fclose(fileID1);



% --------------------------------------------------------------------
function Process_Parmeter_save_Callback(hObject, eventdata, handles)
% hObject    handle to Process_Parmeter_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Process Parameters
materials_Tape=str2num(handles.Process_parameters{1});
materials_sub=str2num(handles.Process_parameters{2});
Velocity=str2double(handles.Process_parameters{3});
Total_energy=str2double(handles.Process_parameters{4});
ID=str2num(handles.Process_parameters{5});   % laser distribution pattern
absorbtion_Refractive_ind=str2num(handles.Process_parameters{6});
Temp_Right_sub=str2num(handles.Process_parameters{7});
Temp_Right_T_Roller=str2num(handles.Process_parameters{8});
h_conv_Sub=str2double(handles.Process_parameters{9});
h_conv_T=str2num(handles.Process_parameters{10});
Roller_Force=str2double(handles.Process_parameters{11});

[file,path] = uiputfile('.\Process_Parameter.txt','Save file name');
% fileID3 = fopen('Comp_Parmeter.txt','w');
fileID2 = fopen(strcat(path,file),'w');
% fileID2 = fopen('Process_Parmeter.txt','w');
if length(ID)>3
    fprintf(fileID2,' materials_Tape, %f %f %f  %f %f \r\n materials_sub, %f %f %f %f %f \r\n Velocity, %f  \r\n Total_energy, %f  \r\n ID, %f %f %f %f %f %f %f \r\n absorbtion_Refractive index, %f %f %f  \r\n Temp_Right_sub-Temp Mandrel, %f %f  \r\n Temp_Right_Tape-Temp-Roller, %f %f  \r\n   h_conv_Sub, %f  \r\n h_conv_T, %f %f  \r\n Roller Force , %f  \r\n',materials_Tape,materials_sub,Velocity,Total_energy,...
    ID,absorbtion_Refractive_ind,Temp_Right_sub,Temp_Right_T_Roller,h_conv_Sub,h_conv_T,Roller_Force);
else
fprintf(fileID2,' materials_Tape, %f %f %f  %f %f \r\n materials_sub, %f %f %f %f %f \r\n Velocity, %f  \r\n Total_energy, %f  \r\n ID, %f %f %f  \r\n absorbtion_Refractive index, %f %f %f   \r\n Temp_Right_sub-Temp Mandrel, %f %f  \r\n Temp_Right_Tape-Temp-Roller, %f %f  \r\n   h_conv_Sub, %f  \r\n h_conv_T, %f %f  \r\n Roller Force , %f  \r\n',materials_Tape,materials_sub,Velocity,Total_energy,...
    ID,absorbtion_Refractive_ind,Temp_Right_sub,Temp_Right_T_Roller,h_conv_Sub,h_conv_T,Roller_Force);
end

fclose(fileID2);

% --------------------------------------------------------------------
function Comp_Parmeter_save_Callback(hObject, eventdata, handles)
% hObject    handle to Comp_Parmeter_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% Computational parameters
No_length=str2double(handles.Computational_parameters{1});  % in radian
No_dev=str2double(handles.Computational_parameters{2});
Width_node=str2double(handles.Computational_parameters{3});  % change to the number of devision
Angular_Node=str2double(handles.Computational_parameters{4});
L_flat_Node=str2double(handles.Computational_parameters{5});  % change to the number of devision


[file,path] = uiputfile('.\Comp_Parameter.txt','Save file name');
% fileID3 = fopen('Comp_Parmeter.txt','w');
fileID3 = fopen(strcat(path,file),'w');
fprintf(fileID3,' No-length, %d  \r\n node_num_Z, %d  \r\n Width_node, %d  \r\n Angular_Node, %d  \r\n L_flat_space, %d  \r\n',No_length,No_dev,Width_node,Angular_Node,L_flat_Node);

fclose(fileID3);



% --------------------------------------------------------------------
function Geo_Parmeter_load_Callback(hObject, eventdata, handles)
% hObject    handle to Geo_Parmeter_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




[fileName,PathName] = uigetfile('*.txt','Select the text file');
fid1 = fopen(strcat(PathName,fileName));

% fid1 = fopen('Geo_Parmeter.txt');
out = textscan(fid1,'%s','delimiter',',');

for ii=1:length(out{1,1})/2
    handles.Geometrical_parameters{ii}=out{1,1}{2*ii,1};
end
% celldisp(handles.Geometrical_parameters)
guidata(hObject, handles);
fclose(fid1);

copyfile(strcat(PathName,fileName),'.\Geo_Parameter.txt');




% --------------------------------------------------------------------
function Process_Parmeter_load_Callback(hObject, eventdata, handles)
% hObject    handle to Process_Parmeter_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



[fileName,PathName] = uigetfile('*.txt','Select the text file');
% fid3 = fopen('Comp_Parmeter.txt');
fid2 = fopen(strcat(PathName,fileName));

% fid2 = fopen('Process_Parmeter.txt');
out = textscan(fid2,'%s','delimiter',',');

for ii=1:length(out{1,1})/2
    handles.Process_parameters{ii}=out{1,1}{2*ii,1};
end
guidata(hObject, handles);

fclose(fid2);


copyfile(strcat(PathName,fileName),'.\Process_Parameter.txt');


% --------------------------------------------------------------------
function Comp_Parmeter_load_Callback(hObject, eventdata, handles)
% hObject    handle to Comp_Parmeter_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName,PathName] = uigetfile('*.txt','Select the text file');
% fid3 = fopen('Comp_Parmeter.txt');
fid3 = fopen(strcat(PathName,fileName));
out = textscan(fid3,'%s','delimiter',',');

for ii=1:length(out{1,1})/2
    handles.Computational_parameters{ii}=out{1,1}{2*ii,1};
end
guidata(hObject, handles);

fclose(fid3);
copyfile(strcat(PathName,fileName),'.\Comp_Parameter.txt');



% --------------------------------------------------------------------
function Untitled_16_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_24_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function New_Objects_Callback(hObject, eventdata, handles)
% hObject    handle to New_Objects (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





% --------------------------------------------------------------------
function Tutorial_Callback(hObject, eventdata, handles)
% hObject    handle to Tutorial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% open Docs\User-manual.pdf;
open Docs\OTOM-V1.06-Guide.pdf;

web('https://www.otomcomposite.eu');


% [cdata,map] = imread('under_construction.png');
%
% h2=msgbox('Under Development!',...
%          'Notice','custom',cdata,map);
%      set(h2,'color','w');


% --------------------------------------------------------------------
function Proxy_connection_Callback(hObject, eventdata, handles)
% hObject    handle to Proxy_connection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[cdata,map] = imread('under_construction.png');

h2=msgbox('Under Development!',...
    'Notice','custom',cdata,map);
set(h2,'color','w');
web('https://www.otomcomposite.eu');


% --------------------------------------------------------------------
function licensing_Callback(hObject, eventdata, handles)
% hObject    handle to licensing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[cdata,map] = imread('under_construction.png');

h2=msgbox('Under Development!',...
    'Notice','custom',cdata,map);
set(h2,'color','w');
web('https://www.otomcomposite.eu');


% --------------------------------------------------------------------
function Mandrel_definition_Callback(hObject, eventdata, handles)
% hObject    handle to Mandrel_definition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% fig=figure('Name','Mandrel_3D_Oval_1','NumberTitle','off');
% fig.ToolBar = 'none';
% fig.MenuBar = 'none';
% set(fig,'color','k');
% [cdata,map] = imread('Docs\latw_mandrel_oval_1.png');
% image(cdata);
%
% fig=figure('Name','Mandrel_3D_Oval_2','NumberTitle','off');
% fig.ToolBar = 'none';
% fig.MenuBar = 'none';
% set(fig,'color','k');
% [cdata,map] = imread('Docs\latw_mandrel_oval_2.jpg');
% image(cdata);
%
% fig=figure('Name','Mandrel_3D_Oval_3','NumberTitle','off');
% fig.ToolBar = 'none';
% fig.MenuBar = 'none';
% set(fig,'color','k');
% [cdata,map] = imread('Docs\latw_mandrel_oval_3.jpg');
% image(cdata);



% [cdata,map] = imread('under_construction.png');
% h2=msgbox('Under Development!',...
%     'Notice','custom',cdata,map);
% set(h2,'color','w');

if ~isempty(handles.Geometrical_parameters)
    % function to display graphic
    
    
    R_cyl= str2num(handles.Geometrical_parameters{11}); %80; %input('The radius of cylinder =');
    z_cyl_end= str2double(handles.Geometrical_parameters{12});
    
    
    % ID=str2num(handles.Process_parameters{5});   % laser distribution pattern
    
    
    
    varargout =  Mandrel_module(R_cyl(1),z_cyl_end);
    
    disp(varargout)
    
    
else
    err_fig=errordlg('there is no geometrical parameters data');
    javaFrame    = get(err_fig,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
end





% --------------------------------------------------------------------
function Report_Callback(hObject, eventdata, handles)
% hObject    handle to Report (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[cdata,map] = imread('under_construction.png');

h2=msgbox('Under Development!',...
    'Notice','custom',cdata,map);
javaFrame    = get(h2,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

set(h2,'color','w');


% --------------------------------------------------------------------
function animation_Callback(hObject, eventdata, handles)
% hObject    handle to animation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% [cdata,map] = imread('under_construction.png');
% 
% h2=msgbox('Under Development!',...
%     'Notice','custom',cdata,map);
% javaFrame    = get(h2,'JavaFrame');
% iconFilePath = 'OTOM-icon.png';
% javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
% 
% set(h2,'color','w');



% --------------------------------------------------------------------
function Rendering_Image_Callback(hObject, eventdata, handles)
% hObject    handle to Rendering_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% [cdata,map] = imread('under_construction.png');
%
% h2=msgbox('Under Development!',...
%     'Notice','custom',cdata,map);
% set(h2,'color','w');

uiwait(warndlg('Print all figures and Save as .png files in \Results\Pics folder!'));

mkdir Results\Pics

% The user should change the image quality

% print('Setup','-dsvg')

print(figure(1),'Results\Pics\Setup.png','-dpng','-r300')


print(figure(2),'Results\Pics\Intensity.png','-dpng','-r300');
print(figure(3),'Results\Pics\Temperature on surfaces.png','-dpng','-r300');

print(figure(21),'Results\Pics\Temperature profile.png','-dpng','-r300');
print(figure(100),'Results\Pics\Temperature 2D.png','-dpng','-r300');

% 3D thermal model figures

Node_z_3D_thermal=handles.Node_z_3D_thermal;

if Node_z_3D_thermal > 1
    print(figure(36),'Results\Pics\3D model Temperature Substrate.png','-dpng','-r300');
    print(figure(46),'Results\Pics\3D model Temperature Up-Bott Substrate.png','-dpng','-r300');
    print(figure(56),'Results\Pics\3D model Temperature profile.png','-dpng','-r300');
    
end




% --------------------------------------------------------------------
function Process_history_Callback(hObject, eventdata, handles)
% hObject    handle to Process_history (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[cdata,map] = imread('under_construction.png');

h2=msgbox('Under Development!',...
    'Notice','custom',cdata,map);
javaFrame    = get(h2,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));


set(h2,'color','w');


% --------------------------------------------------------------------
function Untitled_17_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% h2=msgbox('Under Development!',...
%          'Notice','custom',cdata,map);
%      set(h2,'color','w');


% --------------------------------------------------------------------
function OPC_Communication_Callback(hObject, eventdata, handles)
% hObject    handle to OPC_Communication (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[cdata,map] = imread('under_construction.png');

h2=msgbox('Under Development!',...
    'Notice','custom',cdata,map);
javaFrame    = get(h2,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));


set(h2,'color','w');


% --------------------------------------------------------------------
function Other_Communications_Callback(hObject, eventdata, handles)
% hObject    handle to Other_Communications (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[cdata,map] = imread('under_construction.png');

h2=msgbox('Under Development!',...
    'Notice','custom',cdata,map);
javaFrame    = get(h2,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

set(h2,'color','w');


% --------------------------------------------------------------------
function Thermo_mechanical_Callback(hObject, eventdata, handles)
% hObject    handle to Thermo_mechanical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[cdata,map] = imread('under_construction.png');

h2=msgbox('Under Development!',...
    'Notice','custom',cdata,map);
javaFrame    = get(h2,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

set(h2,'color','w');


% --------------------------------------------------------------------
function Transient_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to Transient_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% [cdata,map] = imread('under_construction.png');
%
% h2=msgbox('Under Development!',...
%     'Notice','custom',cdata,map);
% set(h2,'color','w');
openfig ('.\Transient_control_window');

btn1 = uicontrol('Style', 'pushbutton', 'String', 'Accept',...
    'Position', [10 10 50 20],...
    'Callback', @(btn1,event) OKe_Transient_mode (btn1,hObject, handles));


function OKe_Transient_mode( btn1,hObject, handles)

objects=get(gcf,'Children');

eval([  strcat('handles.',objects(2).Tag),'=objects(2).Value;']);
eval([  strcat('handles.',objects(3).Tag),'=objects(3).Value;']);

eval([  strcat('handles.',objects(5).Tag),'=objects(5).String;']);

eval([  strcat('handles.',objects(6).Tag),'=objects(6).Value;']);

eval([  strcat('handles.',objects(9).Tag),'=objects(9).String;']);
eval([  strcat('handles.',objects(10).Tag),'=objects(10).String;']);


fid13 = fopen('.\Supp_files\Transient_ID.txt','w');
Transient_ID=[handles.checkbox_transient_code, handles.checkbox_Live_output,...
    str2double({handles.edit_init_Temp,...
    handles.edit_inc, handles.edit_Time,handles.checkbox_Video})];

fprintf(fid13,'%f \r\n',Transient_ID);

fclose(fid13);

% for ii=2:length(objects)
%
% eval([  strcat('handles.',objects(ii).Tag),'=objects(ii).Value;']);
%
%
% end

delete(btn1)
guidata(hObject, handles);
savefig('.\Transient_control_window.fig');
close(gcf);
if handles.checkbox_transient_code
    h2=warndlg('In-line monitoring mode is disabled. Input file mode can be Switched ON!');
    javaFrame    = get(h2,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
else
    h2=warndlg('In-line monitoring mode can be Switched ON!');
    javaFrame    = get(h2,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
end







% --------------------------------------------------------------------
function Account_Callback(hObject, eventdata, handles)
% hObject    handle to Account (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% [cdata,map] = imread('under_construction.png');
%
% web('https://www.otomcomposite.eu');
% h2=msgbox('Under Development!',...
%     'Notice','custom',cdata,map);
%   javaFrame    = get(h2,'JavaFrame');
% iconFilePath = 'OTOM-icon.png';
% javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
%
%
% set(h2,'color','w');




% --------------------------------------------------------------------
function Preference_Callback(hObject, eventdata, handles)
% hObject    handle to Preference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[cdata,map] = imread('under_construction.png');
web('https://www.otomcomposite.eu');

h2=msgbox('Under Development!',...
    'Notice','custom',cdata,map);
javaFrame    = get(h2,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));


set(h2,'color','w');


% --------------------------------------------------------------------
function Share_Callback(hObject, eventdata, handles)
% hObject    handle to Share (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% should be placed in Share section of OTOM

  
    




% --------------------------------------------------------------------
function Measuring_Box_Definition_Callback(hObject, eventdata, handles)
% hObject    handle to Measuring_Box_Definition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% [cdata,map] = imread('under_construction.png');


openfig ('.\Measuring_Box_Win');
%       javaFrame    = get(gcf,'JavaFrame');
% iconFilePath = 'OTOM-icon.png';
% javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));



btn1 = uicontrol('Style', 'pushbutton', 'String', 'Accept',...
    'Position', [10 10 50 20],...
    'Callback', @(btn1,event) OKe_Measuring_Box_Definition (btn1,hObject, handles));


function OKe_Measuring_Box_Definition( btn1,hObject, handles)

objects=get(gcf,'Children');



Bly_Sub=objects(3).String;    % Use cooling
Blx_Sub=objects(4).String;  % Sub Initial
c_Sub=objects(5).String;  % Tape inital

Bly_Tape=objects(10).String;    % Use cooling
Blx_Tape=objects(11).String;  % Sub Initial
c_Tape=objects(12).String;

fileID_Post_par = fopen('.\Supp_files\Postprocess_Parameter.txt','w');

fprintf(fileID_Post_par,'Measurement_Box_Tape_c_Blx_Bly, %s %s  %s  \r\n',c_Tape,Blx_Tape,Bly_Tape);
fprintf(fileID_Post_par,'Measurement_Box_Sub_c_Blx_Bly, %s %s  %s  \r\n',c_Sub,Blx_Sub,Bly_Sub);

fclose(fileID_Post_par);

delete(btn1);
%     guidata(hObject, handles);
savefig('.\Measuring_Box_Win.fig')
close(gcf);


% --------------------------------------------------------------------
function Free_winding_Callback(hObject, eventdata, handles)
% hObject    handle to Free_winding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



if ~isempty(handles.Geometrical_parameters)
    % function to display graphic
    th_y=str2double(handles.Geometrical_parameters{1});
    W_tape = str2double(handles.Geometrical_parameters{2}); % width of the tape
    thick_T=str2double(handles.Geometrical_parameters{3}); % thickness of the tape
    R_tape=str2double(handles.Geometrical_parameters{4});   % without thickness, it will be added in general 3D tape cylinder
    L_flat=str2double(handles.Geometrical_parameters{5});
    deg_tape=str2double(handles.Geometrical_parameters{6});
    sui=str2double(handles.Geometrical_parameters{7}) ; %pi/8;  % tape_direction in radian, effect of rotation and sui
    L_prim=str2double(handles.Geometrical_parameters{8});   % tape long
    w=str2double(handles.Geometrical_parameters{9}); % width of the substrate tape
    thick_sub=str2double(handles.Geometrical_parameters{10}); % position of the nip-point z
    R_cyl= str2num(handles.Geometrical_parameters{11}); %80; %input('The radius of cylinder =');
    z_cyl_end= str2double(handles.Geometrical_parameters{12}); %input('Enetr the end point of cylinder =');
    Roller_Pos_TV=str2num(handles.Geometrical_parameters{13});  % for the tape
    W_R=str2double(handles.Geometrical_parameters{14});  % width of the Roller
    Rxyz=str2num(handles.Geometrical_parameters{15});
    Laser_head=str2num(handles.Geometrical_parameters{16});
    L_xyz0=str2num(handles.Geometrical_parameters{17});
    Laser_head_Rot=str2double(handles.Geometrical_parameters{18});
    
    
    N_tape = 360;
    
    Tape_Sp=[N_tape;W_tape;R_tape;L_flat;thick_T;deg_tape];  %Tape_Specification
    
    Roller_Force=str2double(handles.Process_parameters{11});
    
    
    if Roller_Force ~=0
        
        if isempty(handles.fit_Func)
            
            uiwait(warndlg('No Force-displacement data, Please load the data.'));
            
            %                javaFrame    = get(war_fig,'JavaFrame');
            % iconFilePath = 'OTOM-icon.png';
            % javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath))
            
            H_indentation=0;
            %         Roller_def_Callback(hObject, eventdata, handles);
        else
            
            fitresult= handles.fit_Func;
            H_indentation=fitresult(Roller_Force);
            %         figure(50);
            %         plot(Roller_Force,H_indentation,'ko');
            %         text (Roller_Force,H_indentation,' Maximum normal deformation');
            
            
        end
        
    else
        H_indentation=0;
        
    end
    
    if H_indentation > R_tape
        h=msgbox('Normal displacement is more than Roller radius!');
        H_indentation=.8* R_tape;  % 80% as maximum deformation is assumed !!!
    end
    
    
    
    
    
    
    counter_ray=optical_3D_objects_Free_on_Tape(th_y,Tape_Sp,...
        R_cyl,z_cyl_end,Roller_Pos_TV,W_R,H_indentation,L_prim,w,sui,Laser_head);
    
else
    error_fig=errordlg('there is no geometrical parameters data');
    
    javaFrame    = get(error_fig,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
end




%%
% [cdata,map] = imread('under_construction.png');
%
% h2=msgbox('Under Development!',...
%     'Notice','custom',cdata,map);
% set(h2,'color','w');


% --------------------------------------------------------------------
function Roller_subs_fix_winding_Callback(hObject, eventdata, handles)
% hObject    handle to Roller_subs_fix_winding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[cdata,map] = imread('under_construction.png');

h2=msgbox('Under Development!',...
    'Notice','custom',cdata,map);
javaFrame    = get(h2,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));


set(h2,'color','w');


% --------------------------------------------------------------------
function Roller_Laser_fix_winding_Callback(hObject, eventdata, handles)
% hObject    handle to Roller_Laser_fix_winding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[cdata,map] = imread('under_construction.png');

h2=msgbox('Under Development!',...
    'Notice','custom',cdata,map);
javaFrame    = get(h2,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));


set(h2,'color','w');

%%
%
% if ~isempty(handles.Geometrical_parameters)
%     % function to display graphic
%     th_y=str2double(handles.Geometrical_parameters{1});
%     W_tape = str2double(handles.Geometrical_parameters{2}); % width of the tape
%     thick_T=str2double(handles.Geometrical_parameters{3}); % thickness of the tape
%     R_tape=str2double(handles.Geometrical_parameters{4});   % without thickness, it will be added in general 3D tape cylinder
%     L_flat=str2double(handles.Geometrical_parameters{5});
%     deg_tape=str2double(handles.Geometrical_parameters{6});
%     sui=str2double(handles.Geometrical_parameters{7}) ; %pi/8;  % tape_direction in radian, effect of rotation and sui
%     L_prim=str2double(handles.Geometrical_parameters{8});   % tape long
%     w=str2double(handles.Geometrical_parameters{9}); % width of the substrate tape
%     thick_sub=str2double(handles.Geometrical_parameters{10}); % position of the nip-point z
%     R_cyl= str2num(handles.Geometrical_parameters{11}); %80; %input('The radius of cylinder =');
%     z_cyl_end= str2double(handles.Geometrical_parameters{12}); %input('Enetr the end point of cylinder =');
%     Roller_Pos_TV=str2num(handles.Geometrical_parameters{13});  % for the tape
%     W_R=str2double(handles.Geometrical_parameters{14});  % width of the Roller
%     Rxyz=str2num(handles.Geometrical_parameters{15});
%     Laser_head=str2num(handles.Geometrical_parameters{16});
%     L_xyz0=str2num(handles.Geometrical_parameters{17});
%
%     N_tape = 360;
%
%     Tape_Sp=[N_tape;W_tape;R_tape;L_flat;thick_T;deg_tape];  %Tape_Specification
%
%     Roller_Force=str2double(handles.Process_parameters{11});
%
%
%     if Roller_Force ~=0
%
%         if isempty(handles.fit_Func)
%
%             uiwait(warndlg('No Force-displacement data, Please load the data.'));
%             H_indentation=0;
%             %         Roller_def_Callback(hObject, eventdata, handles);
%         else
%
%             fitresult= handles.fit_Func;
%             H_indentation=fitresult(Roller_Force);
%             %         figure(50);
%             %         plot(Roller_Force,H_indentation,'ko');
%             %         text (Roller_Force,H_indentation,' Maximum normal deformation');
%
%
%         end
%
%     else
%         H_indentation=0;
%
%     end
%
%     if H_indentation > R_tape
%         h=msgbox('Normal displacement is more than Roller radius!');
%         H_indentation=.8* R_tape;  % 80% as maximum deformation is assumed !!!
%     end
%
%
%
%
%     counter_ray=optical_3D_objects(th_y,Tape_Sp,...
%         R_cyl,z_cyl_end,Roller_Pos_TV,W_R,H_indentation);
%
% else
%     errordlg('there is no geometrical parameters data');
%
% end


%%


% [cdata,map] = imread('under_construction.png');
%
% h2=msgbox('Under Development!',...
%     'Notice','custom',cdata,map);
% set(h2,'color','w');

% --------------------------------------------------------------------
function Roller_def_Callback(hObject, eventdata, handles)
% hObject    handle to Roller_def (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName,PathName] = uigetfile('*.txt','Select the Force-y_dis text file');
% fid3 = fopen('Comp_Parmeter.txt');
if fileName
    
    fid3 = fopen(strcat(PathName,fileName));
    out = textscan(fid3,'%s','delimiter',',');
    
    length_force_dis=length(out{1,1});
    force=zeros(length_force_dis,1);
    normal_dis=zeros(length_force_dis,1);
    
    for ii=1:length_force_dis
        data=str2num(out{1,1}{ii,1});
        force(ii)=data(1);
        normal_dis(ii)=data(2);
    end
    
    figure(50);
    javaFrame    = get(gcf,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    prompt = {strcat('Polynomial degree- Not more than ',num2str(length_force_dis-1))};
    
    if length_force_dis >8
        prompt = {strcat('Polynomial degree- Not more than ',7)};
        
    end
    
    
    dlg_title = 'fitting option';
    num_lines = 1;
    defaultans = {'3'};
    
    fit_degree = cell2mat(inputdlg(prompt,dlg_title,num_lines,defaultans));
    ft = fittype( strcat('poly',fit_degree) );
    javaFrame    = get(gcf,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    
    % Fit model to data.
    % [fitresult, gof] = fit( force, normal_dis, ft )
    [fitresult] = fit( force, normal_dis, ft );
    
    % plot(fitresult,force, normal_dis,'predfunc');
    plot(fitresult,force, normal_dis,'fit','residuals');
    legend Location SouthWest;
    title('residuals');
    subplot(2,1,1);
    legend Location NorthWest;
    title(formula(fitresult));
    hold on;
    
    % for ii=1:double(fit_degree)
    %     p(ii)=strcat('fitresult.p',num2str(ii));
    %
    % end
    
    % p1 = fitresult.p1;
    % p2 = fitresult.p2;
    % p3 = fitresult.p3;
    % p4 = fitresult.p4;
    
    fclose(fid3);
    handles.fit_Func=fitresult;
    
    guidata(hObject, handles);
    
else
    %     	disp('Error: No file selection ');
    uiwait(msgbox('Error: No file selection','Not complete','modal'));
end


% --------------------------------------------------------------------
function ADS_Client_File_Callback(hObject, eventdata, handles)
% hObject    handle to ADS_Client_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% [ADS_fileName,ADS_PathName] = uigetfile('*.dll','Select ADS Client file');
% fid3 = fopen('Comp_Parmeter.txt');


% try
% if ADS_fileName

%     handles.ADS_fileName=ADS_fileName;
%     handles.ADS_PathName=ADS_PathName;

%     guidata(hObject, handles);


% else

%     uiwait(msgbox('Error: No file selection','Not complete','modal'));
% end
openfig ('.\Inlinemonitoring_control_window');

btn1 = uicontrol('Style', 'pushbutton', 'String', 'Accept',...
    'Position', [10 10 50 20],...
    'Callback', @(btn1,event) OKe_Inlinemonitoring_mode (btn1,hObject, handles));


function OKe_Inlinemonitoring_mode( btn1,hObject, handles)

objects=get(gcf,'Children');

eval([  strcat('handles.',objects(6).Tag),'=objects(6).Value;']);
% eval([  strcat('handles.',objects(3).Tag),'=objects(3).Value;']);
%
% eval([  strcat('handles.',objects(5).Tag),'=objects(5).String;']);
%
% eval([  strcat('handles.',objects(6).Tag),'=objects(6).Value;']);
%
% eval([  strcat('handles.',objects(9).Tag),'=objects(9).String;']);
% eval([  strcat('handles.',objects(10).Tag),'=objects(10).String;']);


% fid13 = fopen('.\Supp_files\Transient_ID.txt','w');
% Transient_ID=[handles.checkbox_transient_code, handles.checkbox_Live_output,...
%                str2double({handles.edit_init_Temp,...
%                handles.edit_inc, handles.edit_Time,handles.checkbox_Video})];

% fprintf(fid13,'%f \r\n',Transient_ID);

% fclose(fid11);

% for ii=2:length(objects)
%
% eval([  strcat('handles.',objects(ii).Tag),'=objects(ii).Value;']);
%
%
% end

delete(btn1)
guidata(hObject, handles);
savefig('.\OKe_Inlinemonitoring_mode.fig')
close(gcf);


% --------------------------------------------------------------------
function Live_Monitoring_Callback(hObject, eventdata, handles)
% hObject    handle to Live_Monitoring (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[cdata,map] = imread('under_construction.png');

web('https://www.otomcomposite.eu/inline-monitoring');
h2=msgbox('Under Development!',...
    'Notice','custom',cdata,map);
javaFrame    = get(h2,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));


set(h2,'color','w');


% show varibles and plot them both in TwinCAT and matlab (input/output) and
% update them very simultaneously
% figure(31);
% fig=figure('Name','Live monitoring','NumberTitle','off');
% % fig = gcf; % current figure handle
% fig.ToolBar = 'none';
% fig.MenuBar = 'none';
% set(fig,'color','k');
% % fig.Title='Live monitoring';
%
%
%
% subplot(2,2,1);
% ax=gca;
% set(ax,'color',[0.2 0.2 0.1]);
% ax.GridColor = 'w';
% ax.XColor = 'y';
% ax.YColor = 'y';
% grid(ax,'on');
% title('Input variable 1','Color','w');
% hold(ax);
% % ax.OuterPosition=[0 0 0.1 0.1];
% % ax.ActivePositionProperty ='outerposition';
%
% %
%
% %
%
% subplot(2,2,2);
% ax=gca;
% set(ax,'color',[0.2 0.2 0.1]);
% ax.GridColor = 'w';
% ax.XColor = 'y';
% ax.YColor = 'y';
% grid(ax,'on');
% title('Output variable 1','Color','w');
% hold(ax);
%
% %
%
% %
%
% subplot(2,2,3);
% ax=gca;
% set(ax,'color',[0.2 0.2 0.1]);
% ax.GridColor = 'w';
% ax.XColor = 'y';
% ax.YColor = 'y';
% grid(ax,'on');
% title('Input variable 2','Color','w');
% hold(ax);
% %
%
% %
%
%
% subplot(2,2,4);
% ax=gca;
% set(ax,'color',[0.2 0.2 0.1]);
% ax.GridColor = 'w';
% ax.XColor = 'y';
% ax.YColor = 'y';
% grid(ax,'on');
% title('Output variable 2','Color','w');
% hold(ax);
%
% %
%
% %
%
%
% % Creat a GIF
% axes('Position' , [0.475 0.92 0.05 0.05]);
% gifplayer;


%%

% fileID1 = fopen('.\Tape_Temp_RedBox.txt');
% fileID2 = fopen('.\Sub3D_Temp_RedBox.txt');
% fileID3 = fopen('.\Inputs.txt');
%
%
% % C = textscan(fileID1,'%d %d %d %d','HeaderLines',1,'Delimiter',',') ;
% Temp_Tape = textscan(fileID1,' %f  ','Delimiter',',','HeaderLines',0) ;
% Temp_Tape_txt=cell2mat(Temp_Tape);
%
% Temp_Sub = textscan(fileID2,' %f  ','Delimiter',',','HeaderLines',0) ;
% Temp_Sub_txt=cell2mat(Temp_Sub);
%
% Inputs = textscan(fileID3,' %f %f  ','Delimiter',',','HeaderLines',1) ;
% Inputs_txt=cell2mat(Inputs);
%
%
% fclose(fileID1);
% fclose(fileID2);
% fclose(fileID3);
%
% Div=40;
%
% if fileID1 >0 && fileID3 >0
%
% y=ones(1,Div)*20;
% y2=y;
%
% t=1:Div;
% % figure;
%  subplot(2,2,1);
% xlabel('Time');
% ylabel('Temperature');
% title('Demo-Tape');
%  PL_Tape=plot(t,y);
%
%  subplot(2,2,2);
%  xlabel('Time');
% ylabel('Temperature');
% title('Demo-Substrate');
%  PL_Sub=plot(t,y2);
%
%
%
%
%  Y_in1=zeros(1,Div);
%   Y_in2=zeros(1,Div);
%
%   subplot(2,2,3);
% xlabel('Time');
% ylabel('Velocity');
% title('Inputs');
%  PL_Vel=plot(t,Y_in1);
%
%  subplot(2,2,4);
%  xlabel('Time');
% ylabel('Power');
% title('Inputs');
%  Pl_Pow=plot(t,Y_in2);
%
%
%
%
%
%
%
%
%
% for ii=1:length(Temp_Tape_txt)
%      subplot(2,2,1);
%      delete (PL_Tape);
%     PL_Tape=plot(t,y,'w');
% %          drawnow;
%     y=[y(2:end) Temp_Tape_txt(ii) ];
%
%      subplot(2,2,2);
%     y2=[y2(2:end) Temp_Sub_txt(ii) ];
%     delete (PL_Sub);
%       PL_Sub=plot(t,y2,'w');
% %       drawnow;
%
%
%        subplot(2,2,3);
%     Y_in1=[Y_in1(2:end) Inputs_txt(ii,1) ];
%     delete (PL_Vel);
%       PL_Vel=plot(t,Y_in1,'w');
% %       drawnow;
%
%
%           subplot(2,2,4);
%     Y_in2=[Y_in2(2:end) Inputs_txt(ii,2) ];
%     delete (Pl_Pow);
%       Pl_Pow=plot(t,Y_in2,'w');
%       drawnow;
%
%
%
%
%
%
%     pause(0.01);
%
% end
% end




% --------------------------------------------------------------------
function ADS_Variables_Callback(hObject, eventdata, handles)
% hObject    handle to ADS_Variables (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



fig=figure('Name','Varibles Connection window','NumberTitle','off');
javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

fig.ToolBar = 'none';
fig.MenuBar = 'none';
set(fig,'color','w');

set(gcf, 'Resize', 'off');



s=cell(10,2);

s{1,1}= 'Process Run ID [INT]'; % indicate whether the process running or not
s{2,1}= 'Velocity [LREAL]';
s{3,1}='Power [LREAL]';
s{4,1}='Tape Temperature [LREAL]';
s{5,1}='Sub Temperature [LREAL]';

s{1,2}='g_o_Process_Run';
s{2,2}='G_XO_PortalSpeed_C1';
s{3,2}='G_XO_LaserPowerActualValue';

s{4,2}='g_o_Tape_Temperature';
s{5,2}='g_o_Sub_Temperature';
s{6,2}='g_o_consolidation_Force';
s{7,2}='g_o_tapeTension';
s{8,2}='g_o_tape_Temp';
s{9,2}= 'g_o_sub_Temp';
s{10,2}='g_o_nip_Temp';




t = uitable(fig,'Data',s,'Position',[10 10 180 220]);
t.ColumnName = {'   OTOM Variables (Fixed)   ','   TwinCat Variables name   '};
t.ColumnEditable(:,2) = true;
t.BackgroundColor(:,2)='g';

% 'Units','normalized'

% fig = gcf;
hTable = get(fig,'children');

hTableExtent = get(hTable,'Extent');
hTablePosition = get(hTable,'Position');
set(hTable,'Position',[50 50 round(1.0*hTableExtent(3)) round(1.5*hTableExtent(4))]);
set(fig,'position',[100 100 round(1.4*hTableExtent(3)) round(2*hTableExtent(4))]);




handles.ADS_table=t;
guidata(hObject, handles);


Edit_text = uicontrol('Style','edit',...
    'String','Enter TwinCat file name',...
    'Position',[200 370 220 30],'Visible','on','Fontsize',12,'BackgroundColor',[0.8 0.8 0.95]);

btn = uicontrol('Style', 'pushbutton', 'String', 'Accept',...
    'Position', [10 10 50 20],...
    'Callback', @(btn,event) Link_variables(btn,t,Edit_text,hObject, handles));



function Link_variables(btn,tbl,Edit_text,hObject, handles)

%     disp(tbl)
% open('ADS_com_Sensors.m');
% choose the path

% [fileName,PathName] = uigetfile('ADS_com_Sensors.m','Select ADS Cconnector file');
%
% fileName='ADS_Com_InlineMonitor.m';
% if fileName
%
%
% addpath(PathName);


%     tmp = PathName;
%     mkdir(tmp);


% runtmp = fullfile(tmp,fileName);
% demodir = fullfile(matlabroot,'toolbox','matlab',...
%     'demos',fileName);
% copyfile(demodir,runtmp);
%     open(strcat(fileName,'(tbl.Data,handles.ADS_fileName,handles.ADS_PathName,Edit_text.String)'));
% else
%
%     uiwait(msgbox('Error: No file selection','Not complete','modal'));
% end

% 3 Tasks should be done in this function :
% 1- recieve variable from TwinCat  >> input collection
% 2-Compute the temperature
% 3- Send back variables and new Temperature to TwinCat >> feedback section



handles.ADS_Vars=tbl.Data;
handles.TwinCat_file_name=Edit_text.String;


guidata(hObject, handles);

%  err=ADS_Com_InlineMonitor(tbl.Data,handles.ADS_fileName,handles.ADS_PathName,Edit_text.String);

close(gcf);
% t=msgbox('Linked!');
% if isempty(err)
uiwait(msgbox('Varaibles Linked for In-line Monitoring!','Success','modal'));
% else
%     msgbox (err.message,'Error','error');
% end












% --------------------------------------------------------------------
function Sub_2D_3D_Callback(hObject, eventdata, handles)
% hObject    handle to Sub_2D_3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%   btn = uicontrol('Style', 'pushbutton', 'String', 'Accept',...
%         'Position', [10 10 50 20],...
%         'Callback', @(btn,event) Link_variables(btn,tbl.Data));
%       btn = uicontrol('Style', 'pushbutton', 'String', 'Accept',...
%         'Position', [10 10 50 20],...
%         'Callback', @(btn,event) Link_variables(btn,tbl.Data));

fig=figure('Name','2D/3D Substrate','NumberTitle','off');

javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

% fig = gcf; % current figure handle
fig.ToolBar = 'none';
fig.MenuBar = 'none';


set(gcf, 'Resize', 'off');


% set(fig,'color','w');
% fig.Title='Live monitoring';


C = imread('3D_slice_4Top.png');
image(C);
axis off;



% subplot(2,2,2);
ax=gca;

pos = get(ax, 'Position');
pos(1) = 0.0;
pos(3) = 1;
set(gca, 'Position', pos)

% set(ax,'color',[0.2 0.2 0.1]);
%  ax.GridColor = 'w';
%  ax.XColor = 'y';
%   ax.YColor = 'y';
%  grid(ax,'on');
title('3D slice and 2D upper surface','Color','b');






text1 = uicontrol('Style','text',...
    'String','Nodes in z-direction?',...
    'Position',[80 80 120 20],'Visible','off');

rQ = uicontrol('Style','edit',...
    'String','5',...
    'Position',[80 55 120 20],'Visible','off');



popup = uicontrol('Style', 'popup',...
    'String', {'2D','3D','2D and 3D'},...
    'Position', [80 280 100 100],...
    'Callback', @(popup, event) popup_func (popup, rQ, text1) );


btn1 = uicontrol('Style', 'pushbutton', 'String', 'Accept',...
    'Position', [250 10 50 20],...
    'Callback', @(btn1,event) Accept (btn1,popup,rQ,hObject, handles));





function popup_func(popup,rQ, text1)
%         x = linspace(0,2*pi,100);
%         y = sin(x);
%         plot(ax,x,y)
if popup.Value > 1
    
    set(text1,'Visible','on');
    set(rQ,'Visible','on');
else
    set(text1,'Visible','off');
    set(rQ,'Visible','off');
    
end


function Accept( btn1,popup,rQ,hObject, handles)

if popup.Value >1
    
    handles.Node_z_3D_thermal=str2double(get(rQ,'String'));
    guidata(hObject, handles);
    
    close(gcf);
else
    handles.Node_z_3D_thermal=0;
        guidata(hObject, handles);
    
    close(gcf);
end


% --------------------------------------------------------------------
function Free_winding_Sub_Callback(hObject, eventdata, handles)
% hObject    handle to Free_winding_Sub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if ~isempty(handles.Geometrical_parameters)
    % function to display graphic
    th_y=str2double(handles.Geometrical_parameters{1});
    W_tape = str2double(handles.Geometrical_parameters{2}); % width of the tape
    thick_T=str2double(handles.Geometrical_parameters{3}); % thickness of the tape
    R_tape=str2double(handles.Geometrical_parameters{4});   % without thickness, it will be added in general 3D tape cylinder
    L_flat=str2double(handles.Geometrical_parameters{5});
    deg_tape=str2double(handles.Geometrical_parameters{6});
    sui=str2double(handles.Geometrical_parameters{7}) ; %pi/8;  % tape_direction in radian, effect of rotation and sui
    L_prim=str2double(handles.Geometrical_parameters{8});   % tape long
    w=str2double(handles.Geometrical_parameters{9}); % width of the substrate tape
    thick_sub=str2double(handles.Geometrical_parameters{10}); % position of the nip-point z
    R_cyl= str2num(handles.Geometrical_parameters{11}); %80; %input('The radius of cylinder =');
    z_cyl_end= str2double(handles.Geometrical_parameters{12}); %input('Enetr the end point of cylinder =');
    Roller_Pos_TV=str2num(handles.Geometrical_parameters{13});  % for the tape
    W_R=str2double(handles.Geometrical_parameters{14});  % width of the Roller
    Rxyz=str2num(handles.Geometrical_parameters{15});
    Laser_head=str2num(handles.Geometrical_parameters{16});
    L_xyz0=str2num(handles.Geometrical_parameters{17});
    Laser_head_Rot=str2double(handles.Geometrical_parameters{18});
    
    N_tape = 360;
    
    Tape_Sp=[N_tape;W_tape;R_tape;L_flat;thick_T;deg_tape];  %Tape_Specification
    
    Roller_Force=str2double(handles.Process_parameters{11});
    
    
    if Roller_Force ~=0
        
        if isempty(handles.fit_Func)
            
            uiwait(warndlg('No Force-displacement data, Please load the data.'));
            %                            javaFrame    = get(wait_fig,'JavaFrame');
            % iconFilePath = 'OTOM-icon.png';
            % javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
            
            H_indentation=0;
            %         Roller_def_Callback(hObject, eventdata, handles);
        else
            
            fitresult= handles.fit_Func;
            H_indentation=fitresult(Roller_Force);
            %         figure(50);
            %         plot(Roller_Force,H_indentation,'ko');
            %         text (Roller_Force,H_indentation,' Maximum normal deformation');
            
            
        end
        
    else
        H_indentation=0;
        
    end
    
    if H_indentation > R_tape
        h=msgbox('Normal displacement is more than Roller radius!');
        H_indentation=.8* R_tape;  % 80% as maximum deformation is assumed !!!
    end
    
    
    
    
    counter_ray=optical_3D_objects_Free_on_Sub(th_y,Tape_Sp,...
        R_cyl,z_cyl_end,Roller_Pos_TV,W_R,H_indentation,L_prim,w,sui,Laser_head);
    
else
    error_fig=errordlg('there is no geometrical parameters data');
    javaFrame    = get(error_fig,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    
end




%%
% [cdata,map] = imread('under_construction.png');
%
% h2=msgbox('Under Development!',...
%     'Notice','custom',cdata,map);
% set(h2,'color','w');


% --------------------------------------------------------------------
function Laser_distribution_control_Callback(hObject, eventdata, handles)
% hObject    handle to Laser_distribution_control (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.Geometrical_parameters)
    % function to display graphic
    
    Rxyz=str2num(handles.Geometrical_parameters{15});
    Laser_head=str2num(handles.Geometrical_parameters{16});
    L_xyz0=str2num(handles.Geometrical_parameters{17});
    
    Rx= Rxyz(1); %-0.8239; %input('direction of ray in x-direction =');
    Ry= Rxyz(2); %-0.1972;% input('direction of ray in y-direction =');
    Rz=Rxyz(3); %-0.3455; %input('direction of ray in z-direction =');
    normR=norm([Rx Ry Rz]); % normalize the Rx,Ry, Rz
    Rx=Rx/normR;
    Ry=Ry/normR;
    Rz=Rz/normR;
    
    
    Ax=Laser_head(1);
    Ay=Laser_head(2);
    nx=Laser_head(3);
    ny=Laser_head(4);
    
    ID=str2num(handles.Process_parameters{5});   % laser distribution pattern
    
    % Gauss_Par_X=10;
    % Gauss_Par_Y=10;
    
    
    % laser_dis_management (Rx,Ry,Rz,Ax,Ay,nx,ny,L_xyz0,ID,Gauss_Par_X,Gauss_Par_Y)
    
    varargout = Laser_dis_manager(Ax,Ay,nx,ny,ID);
    
%     disp(varargout);
    
    % Laser_head(1)=varargout{1};
    % Laser_head(2)=varargout{2};
    % Laser_head(3)=varargout{3};
    % Laser_head(4)=varargout{4};
    
    
    % str2num(handles.Geometrical_parameters{16})=Laser_head
    
    % Gauss_Par_XY=laser_dis_management (Rx,Ry,Rz,Ax,Ay,nx,ny,L_xyz0,ID,Gauss_Par_X,Gauss_Par_Y);
    
else
    error_fig=errordlg('there is no geometrical parameters data');
    javaFrame    = get(error_fig,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
end


% --------------------------------------------------------------------
function Fig_results_Callback(hObject, eventdata, handles)
% hObject    handle to Fig_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% saveas(fig,'MySimulinkDiagram.bmp');

% define the output directory

mkdir Results

saveas(figure(1),'Results\Setup.fig');
saveas(figure(2),'Results\Intensity.fig');
saveas(figure(3),'Results\Temperature on surfaces.fig');

saveas(figure(21),'Results\Temperature profile.fig');
saveas(figure(100),'Results\Temperature 2D.fig');

% 3D thermal model figures

Node_z_3D_thermal=handles.Node_z_3D_thermal;

if Node_z_3D_thermal > 1
    saveas(figure(36),'Results\3D model Temperature Substrate.fig');
    saveas(figure(46),'Results\3D model Temperature Up-Bott Substrate.fig');
    saveas(figure(56),'Results\3D model Temperature profile.fig');
    
end


% --------------------------------------------------------------------
function Figure_open_Callback(hObject, eventdata, handles)
% hObject    handle to Figure_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName,PathName] = uigetfile('*.fig','Select the figure file');
javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

fid1 = openfig(strcat(PathName,fileName));

javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
cameratoolbar('Show');


% --------------------------------------------------------------------
function Input_file_mode_Callback(hObject, eventdata, handles)
% hObject    handle to Input_file_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in OutPut_Control_Butt.
function OutPut_Control_Butt_Callback(hObject, eventdata, handles)
% hObject    handle to OutPut_Control_Butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h_outC=openfig ('.\Output_control');

javaFrame    = get(h_outC,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));


% Check_MB=get(Output_checks.checkbox1,'Value');
% Check_Tecplot=get(Output_checks.checkbox2,'Value');
% Check_Tape_profile=get(Output_checks.checkbox3,'Value');
% Check_Tape_contour=get(Output_checks.checkbox4,'Value');
% % Check_Sub_Profile=get(Output_checks.checkbox5,'Value');
% Check_Sub_contour=get(Output_checks.checkbox6,'Value');
% Check_Config=get(Output_checks.checkbox7,'Value');
% Check_combined_intensity=get(Output_checks.checkbox8,'Value');
%
% Check_Combined_Temp=get(Output_checks.checkbox9,'Value');



btn1 = uicontrol('Style', 'pushbutton', 'String', 'Accept',...
    'Position', [10 10 50 20],...
    'Callback', @(btn1,event) OKe (btn1,hObject, handles));






function OKe( btn1,hObject, handles)

objects=get(gcf,'Children');


for ii=2:length(objects)
    
    eval([  strcat('handles.',objects(ii).Tag),'=objects(ii).Value;']);
    
    
end

delete(btn1)
guidata(hObject, handles);
savefig('.\Output_control.fig')
close(gcf);



% --------------------------------------------------------------------
function Manufacturing_type_Callback(hObject, eventdata, handles)
% hObject    handle to Manufacturing_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% openfig ('.\Manufacturing_Type');
Manufacturing_Type;

javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));



btn1 = uicontrol('Style', 'pushbutton', 'String', 'Accept',...
    'Position', [10 10 50 20],...
    'Callback', @(btn1,event) OKe_Manufacturing_type (btn1,hObject, handles));


function OKe_Manufacturing_type( btn1,hObject, handles)

objects=get(gcf,'Children');



% eval([  strcat('handles.',objects(2).Tag),'=objects(2).String;']); %8 global temp sub
% eval([  strcat('handles.',objects(3).Tag),'=objects(3).String;']); %7 global Temp Tape
% eval([  strcat('handles.',objects(4).Tag),'=objects(4).String;']); %6  initial sub
% eval([  strcat('handles.',objects(5).Tag),'=objects(5).String;']); %5 initial Tape
% eval([  strcat('handles.',objects(10).Tag),'=objects(10).String;']); % Number of measurement
%
% eval([  strcat('handles.',objects(14).Tag),'=objects(14).String;']); %pitch_angle
% eval([  strcat('handles.',objects(17).Tag),'=objects(17).String;']); % delay_time
%
%
% eval([  strcat('handles.',objects(8).Tag),'=objects(8).Value;']); %Radiobottom winding direction, + for 1, 0 for -

% for ii=2:length(objects)
%
% eval([  strcat('handles.',objects(ii).Tag),'=objects(ii).Value;']);
%
%
% end
% manufacturing_type{1}=objects(5).Value;    % Use cooling
% manufacturing_type{2}=objects(11).String;  % Sub Initial
% manufacturing_type{3}=objects(12).String;  % Tape inital
% manufacturing_type{4}=objects(16).String;  % No locations

tags=get(objects,'Tag');
index_radiobutton_use_cooling=find(strcmp(tags, 'radiobutton_use_cooling'));
index_edit_initial_sub=find(strcmp(tags, 'edit_initial_sub'));
index_edit_initial_Tape=find(strcmp(tags, 'edit_initial_Tape'));
index_edit_LocationNo=find(strcmp(tags, 'edit_LocationNo'));

index_edit_roller_properties=find(strcmp(tags, 'edit_roller_properties'));

index_radiobutton_inGasLast=find(strcmp(tags, 'radiobutton_inGasLast'));



index_edit_layerNo=find(strcmp(tags, 'edit_layerNo'));




manufacturing_type{1}=objects(index_radiobutton_use_cooling).Value;    % Use cooling
manufacturing_type{2}=objects(index_edit_initial_sub).String;  % Sub Initial
manufacturing_type{3}=objects(index_edit_initial_Tape).String;  % Tape inital
manufacturing_type{4}=objects(index_edit_LocationNo).String;  % No locations

% Roller properties
manufacturing_type{5}=objects(index_edit_roller_properties).String;  % No locations

manufacturing_type{6}=objects(index_radiobutton_inGasLast).Value;  % No locations


%% 22Aug 2019
manufacturing_type{7}=objects(index_edit_layerNo).String; 




% manufacturing_type{5}=objects(10).String;
%   manufacturing_type{6}=objects(14).String;
%   manufacturing_type{7}=objects(17).String;
%   manufacturing_type{8}=objects(8).Value;

handles.manufacturing_type=manufacturing_type;

delete(btn1)
guidata(hObject, handles);
savefig('.\Manufacturing_Type.fig');
close(gcf);


% --------------------------------------------------------------------
function BRDF_Callback(hObject, eventdata, handles)
% hObject    handle to BRDF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% [cdata,map] = imread('under_construction.png');
%
% h2=msgbox('Under Development!',...
%     'Notice','custom',cdata,map);
%   javaFrame    = get(h2,'JavaFrame');
% iconFilePath = 'OTOM-icon.png';
% javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
%
% set(h2,'color','w');



selection = questdlg('Optical model BRDF mode?',...
    'Warning', ...
    'ON','OFF','ON');

%             javaFrame    = get(selection,'JavaFrame');
% iconFilePath = 'OTOM-icon.png';
% javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));



if selection(1:2)=='ON'
    
    handles.BRDF_mode=1;
    % [fileName,PathName] = uigetfile({'*.xlsx';'*.txt';'*.xslx';'*.csv';'*.xls'},'Select the Data file');
    
    prompt = {'Fiber Rotation (deg)> Substrate,Tape','Division> Substrate,Tape'...
        ,'sig_T> Substrate,Tape','sig_F> Substrate,Tape','Amplitude> Substrate,Tape','Threshold> Substrate,Tape'};
    title = 'BRDF parameters';
    dims = [1 40];
    definput = {'0,90','10 10','0.1 0.1','0.5 0.5','1 1','0.2 0.2' };
    answer = inputdlg(prompt,title,dims,definput);
    
    if ~isempty(answer)
        
        fileID_BRDF_par = fopen('.\Supp_files\BRDF_Par_sub_tape.txt','w');
        
        fprintf(fileID_BRDF_par,'Fiber_Rot_(deg), %s   \r\n',answer{1});
        fprintf(fileID_BRDF_par,'Div, %s   \r\n',answer{2});
        fprintf(fileID_BRDF_par,'sig_T, %s    \r\n',answer{3});
        fprintf(fileID_BRDF_par,'sig_F, %s   \r\n',answer{4});
        fprintf(fileID_BRDF_par,'Amp, %s     \r\n',answer{5});
        fprintf(fileID_BRDF_par,'threshold, %s    \r\n',answer{6});
        
        fclose(fileID_BRDF_par);
        
        %% for more general BRDF
%          prompt = {'Fiber Rotation (deg)> Substrate,Tape','Division_sig_T> Substrate,Tape'...
%         ,'Division_sig_F> Substrate,Tape'...
%         ,'sig_T> Substrate,Tape','sig_F> Substrate,Tape','Amplitude> Substrate,Tape','Threshold> Substrate,Tape'};
%     title = 'BRDF parameters';
%     dims = [1 40];
%     definput = {'0,90','3 3','10 10','0.1 0.1','0.5 0.5','1 1','0.5 0.5' };
%     answer = inputdlg(prompt,title,dims,definput);
%     
%     if ~isempty(answer)
%         
%         fileID_BRDF_par = fopen('.\Supp_files\BRDF_Par_sub_tape.txt','w');
%         
%         fprintf(fileID_BRDF_par,'Fiber_Rot_(deg), %s   \r\n',answer{1});
%         fprintf(fileID_BRDF_par,'Div_sig_T, %s   \r\n',answer{2});
%         fprintf(fileID_BRDF_par,'Div_sig_F, %s   \r\n',answer{3});
%         
%         
%         fprintf(fileID_BRDF_par,'sig_T, %s    \r\n',answer{4});
%         fprintf(fileID_BRDF_par,'sig_F, %s   \r\n',answer{5});
%         fprintf(fileID_BRDF_par,'Amp, %s     \r\n',answer{6});
%         fprintf(fileID_BRDF_par,'threshold, %s    \r\n',answer{7});
%         
        
        %%
        
        
    end
    
    
    
    warn_dlg=warndlg('BRDF mode is turned ON.');
    
    javaFrame    = get(warn_dlg,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    set(handles.uipanel5,'Title','Optical (BRDF) - Thermal')
    
else
    handles.BRDF_mode=0;
    warn_dlg=warndlg('BRDF mode is turned OFF.');
    javaFrame    = get(warn_dlg,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    set(handles.uipanel5,'Title','Optical (Specular) - Thermal')
end

guidata(hObject, handles);


% --------------------------------------------------------------------
function UOT_Callback(hObject, eventdata, handles)
% hObject    handle to UOT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Kinematic_Callback(hObject, eventdata, handles)
% hObject    handle to Kinematic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.Geometrical_parameters)
    % function to display graphic
    th_y=str2double(handles.Geometrical_parameters{1});
    W_tape = str2double(handles.Geometrical_parameters{2}); % width of the tape
    thick_T=str2double(handles.Geometrical_parameters{3}); % thickness of the tape
    R_tape=str2double(handles.Geometrical_parameters{4});   % without thickness, it will be added in general 3D tape cylinder
    L_flat=str2double(handles.Geometrical_parameters{5});
    deg_tape=str2double(handles.Geometrical_parameters{6});
    sui=str2double(handles.Geometrical_parameters{7}) ; %pi/8;  % tape_direction in radian, effect of rotation and sui
    L_prim=str2double(handles.Geometrical_parameters{8});   % tape long
    w=str2double(handles.Geometrical_parameters{9}); % width of the substrate tape
    thick_sub=str2double(handles.Geometrical_parameters{10}); % position of the nip-point z
    R_cyl= str2num(handles.Geometrical_parameters{11}); %80; %input('The radius of cylinder =');
    z_cyl_end= str2double(handles.Geometrical_parameters{12}); %input('Enetr the end point of cylinder =');
    Roller_Pos_TV=str2num(handles.Geometrical_parameters{13});  % for the tape
    W_R=str2double(handles.Geometrical_parameters{14});  % width of the Roller
    Rxyz=str2num(handles.Geometrical_parameters{15});
    Laser_head=str2num(handles.Geometrical_parameters{16});
    L_xyz0=str2num(handles.Geometrical_parameters{17});
    Laser_head_Rot=str2double(handles.Geometrical_parameters{18});
    
    N_tape = 360;
    
    Tape_Sp=[N_tape;W_tape;R_tape;L_flat;thick_T;deg_tape];  %Tape_Specification
    
    Roller_Force=str2double(handles.Process_parameters{11});
    
    
    if Roller_Force ~=0
        
        if isempty(handles.fit_Func)
            
            uiwait(warndlg('No Force-displacement data, Please load the data.'));
            %                            javaFrame    = get(wait_fig,'JavaFrame');
            % iconFilePath = 'OTOM-icon.png';
            % javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
            
            H_indentation=0;
            %         Roller_def_Callback(hObject, eventdata, handles);
        else
            
            fitresult= handles.fit_Func;
            H_indentation=fitresult(Roller_Force);
            %         figure(50);
            %         plot(Roller_Force,H_indentation,'ko');
            %         text (Roller_Force,H_indentation,' Maximum normal deformation');
            
            
        end
        
    else
        H_indentation=0;
        
    end
    
    if H_indentation > R_tape
        h=msgbox('Normal displacement is more than Roller radius!');
        H_indentation=.8* R_tape;  % 80% as maximum deformation is assumed !!!
    end
    
    
    if ~isempty(handles.Computational_parameters)
        No_dev=str2double(handles.Computational_parameters{2});
    else
        No_dev=10;
    end
    
    
    counter_ray=Kinematic_3D_optical(th_y,Tape_Sp,...
        R_cyl,z_cyl_end,Roller_Pos_TV,W_R,H_indentation,L_prim,w,sui,Laser_head,No_dev);
    
else
    error_fig=errordlg('there is no geometrical parameters data');
    javaFrame    = get(error_fig,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    
end


% --------------------------------------------------------------------
function Optical_UOT_Callback(hObject, eventdata, handles)
% hObject    handle to Optical_UOT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selection = questdlg('Unsteady Optical  mode?',...
    'Warning', ...
    'ON','OFF','ON');

%             javaFrame    = get(selection,'JavaFrame');
% iconFilePath = 'OTOM-icon.png';
% javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));



if selection(1:2)=='ON'
    
    [file,path] = uigetfile('.\Kinematic_UOT.mat','Open Kinematics data');
    % fileID3 = fopen('Comp_Parmeter.txt','w');
    str=strcat(path,file);
    % fprintf(fileID3,' No-length, %d  \r\n node_num_Z, %d  \r\n Width_node, %d  \r\n Angular_Node, %d  \r\n L_flat_space, %d  \r\n',No_length,No_dev,Width_node,Angular_Node,L_flat_Node);
    warn_dlg=warndlg('Unsteady Optical  mode is ready to use!, press Calculate');
    javaFrame    = get(warn_dlg,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    
    
    
%         prompt={'Enter a job analysis name'};
% name = 'Analysis Name';
% defaultans = {'Example0'};
% options.Interpreter = 'tex';
% jobname = inputdlg(prompt,name,[1 40],defaultans,options);
    
    
%  jobname='Example0';
% dir=strcat('.\Analysis_UOT\', jobname{:});
% mkdir(dir);   
    
    
    
    handles.UOTmode=1;
%     handles.UOTjobname=jobname{:};
    handles.UOT_optical=str;
    handles.UOT_pathfile=path;
    
    
    set(handles.calculate,'String','Calculate UOT-Optical');
    
    
else
    warn_dlg=warndlg('Unsteady Optical  mode is off.');
    javaFrame    = get(warn_dlg,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    
    handles.UOTmode=0;
    handles.UOT_optical=0;
    set(handles.calculate,'String','Calculate');
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function Thermal_UOT_Callback(hObject, eventdata, handles)
% hObject    handle to Thermal_UOT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function Post_Processing_Callback(hObject, eventdata, handles)
% hObject    handle to Post_Processing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Optical_post_process_Callback(hObject, eventdata, handles)
% hObject    handle to Optical_post_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    
    [file,path_UOT] = uigetfile('.\Optical_UOT.mat','Open Optical Analysis data');
    % fileID3 = fopen('Comp_Parmeter.txt','w');
    str=strcat(path_UOT,file);
  
    load(str);
    Total_energy=str2double(handles.Process_parameters{4});
    
    Optical_post_processing( jobname,R_cyl,z_cyl_end,nip_point_M_all,Rot_Roller_axis_all,...
N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,H_indentation,tv_def_all,...
CV_mesh,path_UOT,Laser_head,Total_energy)

    
    

% --------------------------------------------------------------------
function Thermal_post_process_Callback(hObject, eventdata, handles)
% hObject    handle to Thermal_post_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%   [file,path_UOT] = uigetfile('.\Optical_UOT.mat','Open Optical Analysis data');
    % fileID3 = fopen('Comp_Parmeter.txt','w');
%     str=strcat(path_UOT,file);
%   
%     load(str);
   [file,path_UOT] = uigetfile('.\Optical_UOT.mat','Open Optical Analysis data');
%     [file2,path_UOT2] = uigetfile('.\Pre_Thermal_UOT.mat','Open Pre-thermal Analysis data');
   file2='Pre_Thermal_UOT.mat';
   
    % fileID3 = fopen('Comp_Parmeter.txt','w');
    str=strcat(path_UOT,file);
      str2=strcat(path_UOT,file2);
      
    load(str);
           load(str2);
               Total_energy=str2double(handles.Process_parameters{4});
               
               
            
               
    
    Thermal_post_processing( jobname,R_cyl,z_cyl_end,nip_point_M_all,Rot_Roller_axis_all,...
N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,H_indentation,tv_def_all,...
CV_mesh,path_UOT,Laser_head,...
Points_in_domain_all_sub,Points_in_domain_all_Tape,Total_energy)

    
   


% --------------------------------------------------------------------
function Pre_thermal_UOT_Callback(hObject, eventdata, handles)
% hObject    handle to Pre_thermal_UOT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 selection = questdlg('Pre Thermal Analysis of UOT model?',...
    'Warning', ...
    'ON','OFF','ON');




if selection(1:2)=='ON'
    
    [file,path] = uigetfile('.\Optical_UOT.mat','Open Optical Analysis  data');
 
    str=strcat(path,file);
    
 warn_dlg=warndlg('UOT Pre Thermal  mode is ready to calculate!, Press Calculate');
    javaFrame    = get(warn_dlg,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    
    handles.UOTmode=5;
    handles.UOT_Thermal=str;
    handles.UOT_pathfile=path;
    set(handles.calculate,'String','Calculate UOT-PreThermal');
    
    
else
    warn_dlg=warndlg('Unsteady Pre-Thermal  mode is off.');
    javaFrame    = get(warn_dlg,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    
    handles.UOTmode=0;
    handles.UOT_Thermal=0;
    set(handles.calculate,'String','Calculate');
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function Inputfile_Read_data_Callback(hObject, eventdata, handles)
% hObject    handle to Inputfile_Read_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% select .xls Excel file to read data


% fig=figure('Name','On/off','NumberTitle','off','Position',[10 20 10 200]);
%
%    btn1 = uicontrol('Style', 'radiobutton', 'String', 'Accept',...
%     'Position', [10 10 50 20],...
%     'Callback', @(btn1,event) OKe (btn1,hObject, handles));



selection = questdlg('Input file mode?',...
    'Warning', ...
    'ON','OFF','ON');

%             javaFrame    = get(selection,'JavaFrame');
% iconFilePath = 'OTOM-icon.png';
% javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));



if selection(1:2)=='ON'
    handles.Input_reader_mode=1;
    
    
    selection_read = questdlg('Read new (Velocity, Power) from File?',...
        'Warning', ...
        'Yes','No','Yes');
    
    if selection_read(1:2)=='Ye'
        
        [fileName,PathName] = uigetfile({'*.xlsx';'*.txt';'*.xslx';'*.csv';'*.xls'},'Select the Data file');
        
        prompt = {'Enter sheet number','Start Row','End Row','Start Letter','End Letter','Skip step data'};
        title = 'Input data reader';
        dims = [1 35];
        definput = {'1','1','100','A','B','0'};
        answer = inputdlg(prompt,title,dims,definput);
        
        %         javaFrame    = get(answer,'JavaFrame');
        % iconFilePath = 'OTOM-icon.png';
        % javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
        
        % javaFrame    = get(gcf,'JavaFrame');
        % iconFilePath = 'OTOM-icon.png';
        % javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
        
        sheetNo=str2double(answer{1});
        Start_row=str2double(answer{2});
        end_row=str2double(answer{3});
        Start_Col=answer{4};
        end_Col=answer{5};
        
        workbookFile=strcat(PathName,fileName);
        
        % [VarName1,VarName2] = importfile_Inputs(workbookFile,sheetNo,Start_row,end_row);
        
        [VarName1,VarName2] = importfile_Inputs_Excel(workbookFile,sheetNo,Start_row,end_row,Start_Col,end_Col);
        
        % for large number of data
        skip_step=1+str2double(answer{6});
        
        VarName1=cell2mat(VarName1(1:skip_step:end));
        VarName2=cell2mat(VarName2(1:skip_step:end));
        
        fileID_input_reader = fopen('.\Inputs_mode.txt','w');
        
        fprintf(fileID_input_reader,'%f %f \r\n',[VarName1';VarName2']);
        
        fclose(fileID_input_reader);
        
    end
    
    guidata(hObject, handles);
    warn_dlg=warndlg('Data Input reader mode is turned ON.');
    
    javaFrame    = get(warn_dlg,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
else
    handles.Input_reader_mode=0;
    warn_dlg=warndlg('Data Input reader mode is turned OFF.');
    javaFrame    = get(warn_dlg,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    guidata(hObject, handles);
end



% --------------------------------------------------------------------
function Inputfile_FFT_Callback(hObject, eventdata, handles)
% hObject    handle to Inputfile_FFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_laser_divergence.
function checkbox_laser_divergence_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_laser_divergence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_laser_divergence


% --------------------------------------------------------------------
function UOT_Thermal_Transient_Callback(hObject, eventdata, handles)
% hObject    handle to UOT_Thermal_Transient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 selection = questdlg('Thermal modeling of UOT model?',...
    'Warning', ...
    'ON','OFF','ON');




if selection(1:2)=='ON'
    
    
        handles.step_UOT_steadystate=[];
        
        
        
    [file,path] = uigetfile('.\Optical_UOT.mat','Open Optical Analysis data, Folder should contain Pre-Thermal analysis');
 
    str=strcat(path,file);
    str2=strcat(path,'Pre_Thermal_UOT.mat');
    
    
 warn_dlg=warndlg('UOT Thermal  mode is ready for Transient analysis!, Press Calculate');
    javaFrame    = get(warn_dlg,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    
    handles.UOTmode=2;
    handles.UOT_Thermal=str;
    handles.UOT_PreThermal=str2;
    
    
    handles.UOT_pathfile=path;
    set(handles.calculate,'String','Calculate UOT-Thermal(Tr)');
    
    
else
    warn_dlg=warndlg('Unsteady Thermal  mode is off.');
    javaFrame    = get(warn_dlg,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    
    handles.UOTmode=0;
    handles.UOT_Thermal=0;
    handles.UOT_PreThermal=0;
    
    set(handles.calculate,'String','Calculate');
end
guidata(hObject, handles);



% --------------------------------------------------------------------
function UOT_Thermal_steadystate_Callback(hObject, eventdata, handles)
% hObject    handle to UOT_Thermal_steadystate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 selection = questdlg('Thermal modeling of UOT model?',...
    'Warning', ...
    'ON','OFF','ON');




if selection(1:2)=='ON'
    
    [file,path] = uigetfile('.\Optical_UOT.mat','Open Optical Analysis data, Folder should contain Pre-Thermal analysis');
 
    str=strcat(path,file);
    str2=strcat(path,'Pre_Thermal_UOT.mat');
    
    
 warn_dlg=warndlg('UOT Thermal  mode is ready for Steady-state analysis!, Press Calculate');
    javaFrame    = get(warn_dlg,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    
    handles.UOTmode=2;
    handles.UOT_Thermal=str;
    handles.UOT_PreThermal=str2;
    
    
    prompt = {'Enter your UOT step for steadystate analysis:'};
dlg_title = 'Steps of UOT path';
 defaultans = {'1'};
 Answer= inputdlg(prompt,dlg_title,[1 50],defaultans);
    
    handles.step_UOT_steadystate=Answer;
    
    handles.UOT_pathfile=path;
    set(handles.calculate,'String','Calculate UOT-Thermal(SS)');
    
    
else
    warn_dlg=warndlg('Unsteady Thermal  mode is off.');
    javaFrame    = get(warn_dlg,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
        handles.step_UOT_steadystate=[];
    handles.UOTmode=0;
    handles.UOT_Thermal=0;
    handles.UOT_PreThermal=0;
    
    set(handles.calculate,'String','Calculate');
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function smoothing_Callback(hObject, eventdata, handles)
% hObject    handle to smoothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function smoothing_line_Callback(hObject, eventdata, handles)
% hObject    handle to smoothing_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 warn_dlg=warndlg('Smoothing lines based on Local regression using weighted linear least squares and a 2nd degree polynomial model.');
    javaFrame    = get(warn_dlg,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    

    prompt = {'Enter smoothing degree (more than 0).','Figure number'};
dlg_title = 'Smoothing degree';
 defaultans = {'0.5','1'};
 Answer= inputdlg(prompt,dlg_title,[1 50],defaultans);
 
 smooth_deg=str2double(Answer{1});
  Fig_num=str2double(Answer{2});
%   axes_num=str2double(Answer{3});
  
 h_num=figure(Fig_num);
 
%  axxx=0;
 
 for jj=1:length ((h_num.Children))
     
     s=h_num.Children(jj).Type;
     
     if s(1:2)=='ax'
%          axxx=axxx+1;
%           if axxx==axes_num
         ch=get(h_num.Children(jj),'ch');
         
         for ii=1:length(ch)
             if ch(ii).Type(1:4)=='line'
                 data=ch(ii);
                 
                 % figure;
                 % hold on;
                 % plot(data.XData,data.YData)
                 % [m,n]=size(data);
                 
                 yy1 = smooth(data.XData,data.YData,smooth_deg,'loess');
                 % yy2 = smooth(data.XData,data.YData,10,'rloess');
                 
                 
                 set(ch(ii),'YData',yy1)
                 
             end
         end
%      end
     end
 end
% --------------------------------------------------------------------
function smoothing_contours_Callback(hObject, eventdata, handles)
% hObject    handle to smoothing_contours (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


 warn_dlg=warndlg('Smoothing contours based on Interpolation for 2-D gridded data.');
    javaFrame    = get(warn_dlg,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    

    prompt = {'Figure number','Enter discritization number','contour levels'};
dlg_title = 'Smoothing numbers';
 defaultans = {'1','100','100'};
 Answer= inputdlg(prompt,dlg_title,[1 50],defaultans);
 
  Fig_num=str2double(Answer{1});
 newpoints=str2double(Answer{2});
  smoothing=str2double(Answer{3});

  
 h_num=figure(Fig_num);

 for jj=1:length ((h_num.Children))
     
     s=h_num.Children(jj).Type;
     
     if s(1:2)=='ax'
         
         
         

ch=get(h_num.Children(jj),'ch');
%  CD=get(ch,'CData');
for ii=1:length(ch)
    ss=ch(ii).Type;
    if ss(1:3)=='con' | ss(1:3)=='sur'
x=get(ch(ii),'XData');
y=get(ch(ii),'YData');
z=get(ch(ii),'ZData');

% [m,n]=size(x)
xrow=x;
ycol=y;
BDmatrix=z;

% newpoints = 100;
[xq,yq] = meshgrid(...
            linspace(min(min(xrow,[],2)),max(max(xrow,[],2)),newpoints ),...
            linspace(min(min(ycol,[],1)),max(max(ycol,[],1)),newpoints )...
          );
BDmatrixq = interp2(xrow,ycol,BDmatrix,xq,yq,'spline');
figure

[c,h]=contourf(xq,yq,BDmatrixq,smoothing,'linestyle','none');
colormap(jet(smoothing));
axis equal
    end
end
    javaFrame    = get(gcf,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    
     end
 end


% --------------------------------------------------------------------
function UOT_Optimization_Callback(hObject, eventdata, handles)
% hObject    handle to UOT_Optimization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if (handles.UOTmode)==0

%  selection = questdlg('Thermal modeling of UOT model?',...
%     'Warning', ...
%     'ON','OFF','ON');




% if selection(1:2)=='ON'
    
    [file,path] = uigetfile('.\Optical_UOT.mat','Open Optical Analysis data, Folder should contain Pre-Thermal analysis');
 
    str=strcat(path,file);
    str2=strcat(path,'Pre_Thermal_UOT.mat');
    
    
 warn_dlg=warndlg('UOT Thermal  mode is ready for optimization, Press Run UOT! ');
    javaFrame    = get(warn_dlg,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    
     handles.UOT_Thermal=str;
    handles.UOT_PreThermal=str2;
        handles.UOT_pathfile=path;
end
    % for optimization 
    handles.UOTmode=10;
   
    
    
    prompt = {'Starting step:','Finishing step:'};
dlg_title = 'Steps of UOT optimization';
 defaultans = {'6','10'};
 Answer= inputdlg(prompt,dlg_title,[1 50],defaultans);
    
    handles.step_UOT_optimization=Answer;
    

%     set(handles.calculate,'String','Calculate UOT-Thermal(SS)');
    
        set(handles.OP_Run,'String','Run UOT');
% else
%     warn_dlg=warndlg('UOT optimization  mode is off.');
%     javaFrame    = get(warn_dlg,'JavaFrame');
%     iconFilePath = 'OTOM-icon.png';
%     javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
%     
%          handles.step_UOT_optimization=[];
%     handles.UOTmode=0;
%     handles.UOT_Thermal=0;
%     handles.UOT_PreThermal=0;
%     
%     set(handles.OP_Run,'String','Run');


%      str= handles.UOT_Thermal;
%     str2= handles.UOT_PreThermal;
%     load(str);
%     load(str2);    
%    Ax=Laser_head(1);
% Ay=Laser_head(2);
    


% s=cell(1,8);
% s{1,1}='Laser-ID 3';
set(handles.OP_popupmenu1, 'value', 1);
% set(handles.OP_popupmenu1,'String',s) ;
% 
% s(1,2:8)=[];
% 
% 
% set(handles.OP_popupmenu1,'String',s) ;
set(handles.OP_popupmenu1,'String','Laser-ID 3- Power') ;



guidata(hObject, handles);































% prompt = {'Lower Value','Upper Value'};
% dlg_title = 'Range';
% num_lines = 1;

% Laser_head=str2num(handles.Geometrical_parameters{16});



%  str= handles.UOT_Thermal;
%     str2= handles.UOT_PreThermal;
%     load(str);
%     load(str2);
%     
%     
%    Ax=Laser_head(1);
% Ay=Laser_head(2);
%     
% Ax=Laser_head(1);
% Ay=Laser_head(2);
% 
% 
% ID_parameter=
% set(handles.OP_popupmenu1,'String','Laser-ID') ;
% 
% switch  ID_parameter
%     case   1
%                 th_y=xin;
%        
%         defaultans = {['3 0 0 0 0 ', num2str(-Ax) ,' ' , num2str(-Ay)  ],['3 0 0 0 0 ', num2str(-Ax) ,' ' , num2str(-Ay)  ]};
%     case 2
%                 Velocity=xin;
%         defaultans = {'0.05','1'};
%     case 3
%                 Total_energy=xin;
%         defaultans = {'1e2','9e2'};
%     case 4
%                 L_xyz0=xin;
%         defaultans = {'0.100000, 0.150000, 0.250000  ','0.900000, 0.950000, 2.250000 '};
%     case 5
%         Laser_direction (Rx,Ry,Rz)
%               Rxyz  =xin;
%         defaultans = {'-1, -1, 1  ','1, 1, 1 '};
%     case 6
%         Laser Head Size (Ax,Ay)
%               Laser_head(1:2)  =xin;
%         defaultans = {'0.100000, 0.150000  ','0.100000, 0.20000'};
%         
%     case 7
%         Laser-ID parameter  (Par_x, Par_Y)
%               ID (2:3)  =xin;
%         defaultans = {'10, 20  ','20, 60'};
%     case 8
%         Mandrel radius R_cyl
%              R_cyl (1)   =xin;
%         defaultans = {'0.100000 ','0.600000'};
%         
%         
%         
%         
% end

% defaultans = {'-180','-110'};

% options.Interpreter = 'tex';
% handles.Range = inputdlg(prompt,dlg_title,num_lines,defaultans,options);
% javaFrame    = get(gcf,'JavaFrame');
% iconFilePath = 'OTOM-icon.png';
% javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

% guidata(hObject, handles);



function edit_Absorbed_power_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Absorbed_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Absorbed_power as text
%        str2double(get(hObject,'String')) returns contents of edit_Absorbed_power as a double


% --- Executes during object creation, after setting all properties.
function edit_Absorbed_power_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Absorbed_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Absorbed_power_T_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Absorbed_power_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Absorbed_power_T as text
%        str2double(get(hObject,'String')) returns contents of edit_Absorbed_power_T as a double


% --- Executes during object creation, after setting all properties.
function edit_Absorbed_power_T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Absorbed_power_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Gen_Pareto_Result_Callback(hObject, eventdata, handles)
% hObject    handle to Gen_Pareto_Result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.UOTmode==10
  str= handles.UOT_Thermal;
str2= handles.UOT_PreThermal;

else
   [file,path] = uigetfile('.\Optical_UOT.mat','Open Optical Analysis data, Folder should contain Pareto optimization results');
  str=strcat(path,file);
    str2=strcat(path,'Pre_Thermal_UOT.mat');
    
    
 warn_dlg=warndlg('UOT Thermal  mode is ready for Pareto generation, Press Run UOT! ');
    javaFrame    = get(warn_dlg,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    
     handles.UOT_Thermal=str;
    handles.UOT_PreThermal=str2;
        handles.UOT_pathfile=path;
        
        
          handles.UOTmode=10;
    prompt = {'Starting step:','Finishing step:'};
dlg_title = 'Steps of UOT Pareto generation from  optimization ';
 defaultans = {'6','10'};
 Answer= inputdlg(prompt,dlg_title,[1 50],defaultans);
    
    handles.step_UOT_optimization=Answer;
    

%     set(handles.calculate,'String','Calculate UOT-Thermal(SS)');
    
        set(handles.OP_Run,'String','Run UOT');
end

    prompt = {'Objective number'};
dlg_title = 'Select  Objective number, use 2 Objectives for their average';
 defaultans = {'1'};
 Answer= inputdlg(prompt,dlg_title,[1 20],defaultans);
    
    handles.Objective_num=Answer;

    
guidata(hObject, handles);


load(str);
load(str2);

   callObjConstr_4Gui_UOT_Pareto_Gen(handles,...
    th_y,W_tape,thick_T,R_tape,L_flat,deg_tape,...
           R_cyl,...
            W_R,...
            H_indentation,...
                nip_point_M_all,Rot_Roller_axis_all,...
            CV_mesh,Laser_head,...
            Points_in_domain_all_sub,Points_in_domain_all_Tape,Tape_points_Data);


% --- Executes on button press in OP_exact_power.
function OP_exact_power_Callback(hObject, eventdata, handles)
% hObject    handle to OP_exact_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OP_exact_power


% --------------------------------------------------------------------
function Selected_UOT_Laser_Callback(hObject, eventdata, handles)
% hObject    handle to Selected_UOT_Laser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.UOTmode==10
  str= handles.UOT_Thermal;
str2= handles.UOT_PreThermal;

else
   [file,path] = uigetfile('.\Optical_UOT.mat','Open Optical Analysis data, Folder should contain Pareto optimization results');
  str=strcat(path,file);
    str2=strcat(path,'Pre_Thermal_UOT.mat');
    
    
%  warn_dlg=warndlg('UOT Thermal  mode is ready for Pareto generation, Press Run UOT! ');
%     javaFrame    = get(warn_dlg,'JavaFrame');
%     iconFilePath = 'OTOM-icon.png';
%     javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    
    handles.UOT_Thermal=str;
    handles.UOT_PreThermal=str2;
    handles.UOT_pathfile=path;
        
        
          handles.UOTmode=10;
    prompt = {'Starting step:','Finishing step:'};
dlg_title = 'Steps of UOT Pareto generation from  optimization ';
 defaultans = {'6','10'};
 Answer= inputdlg(prompt,dlg_title,[1 50],defaultans);
    
    handles.step_UOT_optimization=Answer;
    

%     set(handles.calculate,'String','Calculate UOT-Thermal(SS)');
    
        set(handles.OP_Run,'String','Run UOT');
end

%     prompt = {'Objective number'};
% dlg_title = 'Select  Objective number, use 2 Objectives for their average';
%  defaultans = {'1'};
%  Answer= inputdlg(prompt,dlg_title,[1 20],defaultans);
%     
%     handles.Objective_num=Answer;

guidata(hObject, handles);


load(str);
load(str2);

   Rendering_optimization_laser_intensity(handles,...
              Laser_head);


% --- Executes on button press in checkbox_nonlin_thermal.
function checkbox_nonlin_thermal_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_nonlin_thermal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_nonlin_thermal
% [cdata,map] = imread('complete.png');
val=get(handles.checkbox_nonlin_thermal,'Value');

if val
h2=msgbox('Only Active for Transient thermal model!',...
    'Active');

javaFrame    = get(h2,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));


 prompt = {'Tape Material name','Substrate Material name'};
    title = 'Thermal properties (should exist)';
    dims = [1 60];
    definput = {'Mat1','Mat1' };
    answer = inputdlg(prompt,title,dims,definput);

handles.nonlinearMaterials_Tape_Sub=[answer];

guidata(hObject, handles);

else
   h2=msgbox('Deactive for Transient thermal model!',...
    'Deactive');

javaFrame    = get(h2,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath)); 

    handles.nonlinearMaterials_Tape_Sub=[];

guidata(hObject, handles);

end

% --------------------------------------------------------------------
function Non_lin_thermal_Callback(hObject, eventdata, handles)
% hObject    handle to Non_lin_thermal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get the data and do the fitting  and give a flag for non-linear modeling
% !!
   [file,path] = uigetfile('.xlsx','Open Excel file for temperature-dependent thermal properties');
  str=strcat(path,file);
%     str2=strcat(path,'Pre_Thermal_UOT.mat');

% Read the data folder  and show them 

t = readtable(str);
vars = {'Temperature','Cp','Density','Kx','Ky','Kz'};
% t = t(1:9,vars);
% fig = uifigure;

fig=figure('Name','Thermal property Data','NumberTitle','off');
javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

% fig = gcf; % current figure handle
fig.ToolBar = 'none';
fig.MenuBar = 'none';
set(gcf, 'Resize', 'off');
uit = uitable(fig,'Data',t{:,:},'ColumnName',vars,'Units', 'normalized', 'Position',[0, 0, 1, 1]);
% , 'Position',[0, 0, 1, 1]);
uit.ColumnEditable = true;
% fig.Units='normalized';
% , 'Position',[0, 0, 1, 1]


 Poly_fit_degree = uicontrol(fig, 'Style', 'edit', 'String', '3',...
        'Position', [150 20 80 20]);
     Text_poly = uicontrol(fig, 'Style', 'text', 'String', 'Polynomial degree',...
        'Position', [150 40 100 20],'Backgroundcolor','w');
    
    
     Material_name = uicontrol(fig, 'Style', 'edit', 'String', 'Mat1',...
        'Position', [250 20 80 20]);
     Text_Material_name = uicontrol(fig, 'Style', 'text', 'String', 'Material name',...
        'Position', [250 40 100 20],'Backgroundcolor','w');
     
    btn_save_fitting = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Save Fitting',...
        'Position', [10 20 90 20]     ,...
     'Callback', @(btn_save_fitting,event) save_fitting(Poly_fit_degree.String,uit.Data,fig,Material_name.String));
    
% s = uistyle('BackgroundColor','green');
%  addStyle(uit,s,'cell',[1,1]);
