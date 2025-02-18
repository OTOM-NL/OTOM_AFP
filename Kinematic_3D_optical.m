
% This program is able show 3D objects

function counter_ray=Kinematic_3D_optical(th_y,Tape_Sp,...
    R_cyl,z_cyl_end,Roller_Pos_TV,W_R,H_indentation,...
    L_prim,w,sui,Laser_head,...
    No_dev)




counter_ray=1;
% Laser_head,L_xyz0

% Laser_Target=[0,0,0];


N_tape=Tape_Sp(1); % EVEN number of points, should not be changed !!
W_tape=Tape_Sp(2) ; % width of the tape
R_tape=Tape_Sp(3);
L_flat=Tape_Sp(4);
thick_T=Tape_Sp(5); % thickness of the tape
deg_tape=Tape_Sp(6);

Ax=Laser_head(1);
Ay=Laser_head(2);
nx=Laser_head(3);
ny=Laser_head(4);

% Create figure
fig_HO=figure('Name','Kinematic on Dome-Cylinder','NumberTitle','off');
cameratoolbar(fig_HO, 'Show');
javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

% create structure of handles
% Fig_handles = guihandles(fig_HO);
% % Add some additional data as a new field called numberOfErrors
% % Fig_handles.numberOfErrors = 0;
% % Save the structure
% guidata(fig_HO,Fig_handles)



Fig_handles = guidata(gcbo);

% Fig_handles.Temp=3;





% fig = gcf; % current figure handle
% fig.ToolBar = 'none';
% fig.MenuBar = 'none';
hold on;


% Output is a handle to be manupulated
Mandrel_plot=plot3D_cylinder4optical_3D_objects(R_cyl,z_cyl_end);


%%
dcm_obj = datacursormode(fig_HO);
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on');


set(dcm_obj,'UpdateFcn',@NewCallback_datacursor);
hTarget = handle(Mandrel_plot(3));
hDatatip = dcm_obj.createDatatip(hTarget);
set(hDatatip,'Host',hTarget);
datacursormode toggle;







%%




theta_ind=acosd((R_tape-H_indentation)/R_tape);  %theta_ind in degree

Nip_Mov =R_tape*sind(theta_ind);

%%
%  No_discrit=50; % change the step size
%  No_dev=10;
%  No_dev_L=ceil(abs(L_prim*No_discrit)); % for representation
No_dev_L=50; %

Fig_handles.No_dev_L=No_dev_L;



nip_point_M=[0,R_cyl(1),z_cyl_end];
%   if  R_cyl(1) automatically goes to the center of placement
% upper dome >> R_cyl(3)
% Bottom dome >> R_cyl(2)

% Winding direction;
% m=0.05;

[Fig_handles.h_path,Fig_handles.Points_center_path,Fig_handles.Lag_points, Boarders,starting_index]=Winding_path_cartesian_adaptive (R_cyl,z_cyl_end,L_prim,w,-abs(mod(th_y,90))-90,No_dev_L,No_dev,Nip_Mov,nip_point_M);

%  h_path;
Points_center_path=Fig_handles.Points_center_path;
Lag_points=Fig_handles.Lag_points;


% laser direction line for deleting
Fig_handles.h_laser_dir=[];


axis equal;
axis off;
%   nip_point_M=Roller_Pos_TV;

% I should make Eulerian mesh but Lagrangian nodes
% starting_index=floor(No_dev_L*(Nip_Mov/abs(L)));

% Number of node in Eulerian description
%% index of CV to be considered in optical model >> depends on length of
% laser size
Long_CV=20;

xnode=2*No_dev -1;
ynode=Long_CV ; %No_dev_L-starting_index;

X=zeros(ynode,xnode);
Y=zeros(ynode,xnode);
Z=zeros(ynode,xnode);
Intensity=Z;

for kk=1:xnode
    %     index=(((1:ynode)*xnode)-kk+1);
    X(1:end,kk)=Lag_points(1,1:ynode,kk);
    Y(1:end,kk)=Lag_points(2,1:ynode,kk);
    Z(1:end,kk)=Lag_points(3,1:ynode,kk);
    Intensity(1:end,kk)=0; % RHS in the thermal model
end
% pause(0.5);

h_Eulerian=mesh(X,Y,Z,Intensity,'FaceColor','None','EdgeColor','y');
%%


Fig_handles.h_Eulerian=h_Eulerian;
% Fig_handles.numberOfErrors = 0;
% Save the structure
guidata(fig_HO,Fig_handles) ;

nip_point_M=[Points_center_path(1,1),Points_center_path(2,1),Points_center_path(3,1)];
tv=nip_point_M;  % for the tape
% tv_R=[0;R_cyl(1)+R_tape+thick_T-H_indentation;tv3];  % for center of Roller

%    if index < length(Points_center_path)
%      Vec_tan=Points_center_path(1:3,index)-Points_center_path(1:3,index+1);
%
%      else
Vec_tan= Points_center_path(1:3,1)-Points_center_path(1:3,2);
%      end

% indicate 0 or 180 CW or CCW
Dir=0;

[Tape_Roller_Handle,Rot_Roller_axis]=Tape_plot_3D_ALL_objects_Dome_Car(N_tape,W_tape,R_tape,L_flat,tv,-abs(mod(th_y,90))+90,thick_T,deg_tape,W_R,theta_ind,R_cyl,z_cyl_end,H_indentation,Vec_tan,Dir);
Fig_handles.Rot_Roller_axis=Rot_Roller_axis;

Fig_handles.Tape_Roller_Handle=Tape_Roller_Handle;

Fig_handles.Dir=0;


% Fig_handles.numberOfErrors = 0;
% Save the structure
guidata(fig_HO,Fig_handles) ;



camlight(-90,90);

dcm_obj = datacursormode(fig_HO);
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on');


set(dcm_obj,'UpdateFcn',@NewCallback_datacursor);
hTarget = handle(h_Eulerian);
hDatatip = dcm_obj.createDatatip(hTarget);
set(hDatatip,'Host',hTarget);
datacursormode toggle;



c_info = getCursorInfo(dcm_obj);


Laser_Target=c_info.Position;



Fig_handles.hplot_start=[];
Fig_handles.hplot_finish=[];
Fig_handles.finish_index=[];
Fig_handles.start_index=[];

Fig_handles.Laser_direction=[];
Fig_handles.Laser_Head_L_xyz0=[0 0 0];
Fig_handles.Laser_Target=Laser_Target;
Fig_handles.nip_point_M=nip_point_M;

guidata(gcbo,Fig_handles) ;


% plot3(Laser_Target(1),Laser_Target(2),Laser_Target(3),'r*');

% [points_T_G,z,dl,th_v]=Tape_points (deg_tape,R_tape,W_tape,L_flat,tv,th_y,W_R,...
%     2,5,3,theta_ind)
%
%
% node_num_Z_T=length(z);
%
%
% xnode_T=node_num_Z_T;
% ynode_T=length  (th_v)+length(dl);
%
% X_T=zeros(xnode_T,ynode_T);
% Y_T=zeros(xnode_T,ynode_T);
% Z_T=zeros(xnode_T,ynode_T);
% Intensity_T=zeros(xnode_T,ynode_T);
%
% for kk=1:xnode_T
% %     index=((1:ynode_T)*xnode_T)-kk+1;
%     X_T(kk,1:end)=points_T_G(kk,1)';
%     Y_T(kk,1:end)=points_T_G(kk,2)';
%     Z_T(kk,1:end)=points_T_G(kk,3)';
%      Intensity_T(kk,1:end)=1;
% end
%
% surf(X_T,Y_T,Z_T,Intensity_T);




axis equal;

% grid  on;
% axis([-2*R_cyl 2*R_cyl -2*R_cyl 2*R_cyl -z_cyl_end/2 1.5*z_cyl_end]);
view([225 300]);


distance=R_cyl(1);
% r=distance;

% tt=linspace(0,2*pi,50);
% xx=r*sin(tt);
% yy=r*cos(tt);
% zz=yy*0;

% Circle_Rotation=-(pi/180)*th_y ;


% new_xyz=Transformation_Rot_y(Laser_Target,Circle_Rotation,[xx;yy;zz]);


[x,y,z]=sphere(50);
x=x*distance + (Laser_Target(1));
y=y*distance + (Laser_Target(2));
z=z*distance + (Laser_Target(3));

circle_Laser=mesh(x,y,z,'FaceAlpha',0.1,'EdgeAlpha',0.1,'FaceColor','g');

% circle_Laser=plot3(new_xyz(1,:),new_xyz(2,:),new_xyz(3,:),'g--');


dcm_obj = datacursormode(fig_HO);
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on');


set(dcm_obj,'UpdateFcn',@NewCallback_datacursor);

%% make target

hTarget = handle(circle_Laser);
hDatatip = dcm_obj.createDatatip(hTarget);
set(hDatatip,'Host',hTarget);
datacursormode toggle;


%%

a = findall(gcf);
b = findall(a,'ToolTipString','Data Cursor');
set(b,'Visible','Off');



text_orientation = uicontrol('Style','text',...
    'String',sprintf(' nx=%f \n ny=%f \n nz=%f',0 ,0,0),'Fontsize',12,...
    'Position',[10 300 200 100],'Visible','on');



text1 = uicontrol('Style','text',...
    'String','Laser Distance',...
    'Position',[80 80 120 20],'Visible','on');

rQ = uicontrol('Style','edit',...
    'String',num2str(R_cyl(1)),...
    'Position',[80 55 120 20],'Visible','on');



text_Lasertarget = uicontrol('Style','text',...
    'String','Target',...
        'Position', [200 75 120 20],'Visible','on');

rQ_Lasertarget = uicontrol('Style','edit',...
    'String','0 0 0',...
    'Position',[200 55 120 20],'Visible','on');

text_Lasersource = uicontrol('Style','text',...
    'String','Source',...
        'Position', [320 75 120 20],'Visible','on');
    
rQ_Source= uicontrol('Style','edit',...
    'String','0 0 0',...
    'Position',[320 55 120 20],'Visible','on');


text_nip_point = uicontrol('Style','text',...
    'String','Nip point',...
        'Position', [450 75 120 20],'Visible','on');
    
rQ_nip_point= uicontrol('Style','edit',...
    'String','Inf',...
    'Position',[450 55 120 20],'Visible','on');



btn_save_editbox = uicontrol('Style', 'pushbutton', 'String', 'Save editBox',...
    'Position', [580 75 120 20],...
    'Callback', @(btn_save_editbox,event) Save_editbox (btn_save_editbox,rQ_Lasertarget,rQ_Source,rQ_nip_point,rQ));

btn_load_editbox = uicontrol('Style', 'pushbutton', 'String', 'Load editBox',...
    'Position', [580 55 120 20],...
    'Callback', @(btn_load_editbox,event) Load_editbox (btn_load_editbox,rQ_Lasertarget,rQ_Source,rQ_nip_point,rQ));


%% No.discritization along the length
text_dis_path = uicontrol('Style','text',...
    'String','No.On path',...
    'Units','normalized'...
    ,'Position', [0 0.94 0.08 0.02],'Visible','on');

rQ_dis_path = uicontrol('Style','edit',...
    'String',50,...
    'Units','normalized'...
    ,'Position', [0 0.90 0.08 0.02],'Visible','on');


text_Eulerian = uicontrol('Style','text',...
    'String','No.Eulerian',...
    'Units','normalized'...
    ,'Position', [0 0.85 0.08 0.02],'Visible','on');

rQ_Eulerian = uicontrol('Style','edit',...
    'String',20,...
    'Units','normalized'...
    ,'Position', [0 0.81 0.08 0.02],'Visible','on');







btn1 = uicontrol('Style', 'pushbutton', 'String', 'OKe',...
    'Position', [80 10 100 20],...
    'Callback', @(btn1,event) Accept (btn1,rQ,th_y,fig_HO,Laser_Target));

btn2 = uicontrol('Style', 'pushbutton', 'String', 'Put Laser',...
    'Position', [300 10 100 20],...
    'Callback', @(btn2,event) Put (btn2,fig_HO,text_orientation,rQ_Lasertarget,rQ_Source));

btn3 = uicontrol('Style', 'pushbutton', 'String', 'Head Align Target',...
    'Position', [500 10 150 20],...
    'Callback', @(btn3,event) Tape_target (btn3,rQ,th_y,fig_HO,circle_Laser,rQ_Lasertarget,rQ_Source));


btn4 = uicontrol('Style', 'pushbutton', 'String', 'Change point of path','Units','normalized'...
    ,'Position', [0.8 0.8 0.2 0.05],...
    'Callback', @(btn4,event) Select_Nip_point (btn4,fig_HO) );



% for adaptive mode
Laser_Head_Checkbox = uicontrol('Style','checkbox',...
    'String','Adaptive Laser Pos\Dir',...
    'Position',[40 180 180 20],'Visible','on');



text_Rot = uicontrol('Style','text',...
    'String','Laser Head Rotation',...
    'Position',[40 220 180 20],'Visible','on');


sld = uicontrol('Style', 'slider',...
    'Min',0,'Max',180,'Value',90,...
    'Position', [40 200 180 20],...
    'Callback', @(sld,event) Laser_surf (text_Rot,sld,fig_HO,Ax,Ay,rQ_Lasertarget,rQ_Source) );


% save path info between a start point and end point


btn_save_path_info = uicontrol('Style', 'pushbutton', 'String', 'Save path info','Units','normalized'...
    ,'Position', [0.8 0.1 0.1 0.05],...
    'Callback', @(btn_save_path_info,event) Save_intxt (btn_save_path_info, th_y,R_cyl,...
    z_cyl_end,No_dev,Ax,Ay,nx,ny));




max_point=length(Points_center_path);

sld_nipOnPath = uicontrol('Style', 'slider','Units','normalized'...
    ,'Position', [0 0.97 1 0.03],'SliderStep', [1/max_point , 10/max_point ],...
    'Min',1,'Max',max_point,'Value',1,...
    'Callback', @(sld_nipOnPath,event) nipOnPath (sld_nipOnPath,fig_HO,Fig_handles.h_Eulerian,...
    N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,R_cyl,...
    z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,...
    rQ,H_indentation,...
    rQ_Eulerian,Laser_Head_Checkbox) );

btn5 = uicontrol('Style', 'pushbutton', 'String', 'Accept New path','Units','normalized'...
    ,'Position', [0.8 0.7 0.1 0.05],...
    'Callback', @(btn5,event) Accept_Nip_point (btn5,fig_HO,...
    Fig_handles.Tape_Roller_Handle,...
    N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,R_cyl,...
    z_cyl_end,L_prim,w,sui,No_dev,Nip_Mov,...
    rQ,H_indentation,...
    rQ_dis_path,rQ_Eulerian,sld_nipOnPath,rQ_nip_point) );




popup_dir = uicontrol('Style','popup',...
    'String',{'Forward','Backward'},...
    'Units','normalized'...
    ,'Position', [0 0.75 0.08 0.02],'Callback', @(popup_dir,event) Direction(popup_dir));


%Save start and end index based on nip-point path

btn_start = uicontrol('Style', 'pushbutton', 'String', 'Set start','Units','normalized'...
    ,'Position', [0.8 0.15 0.1 0.02],...
    'Callback', @(btn_start,event) set_start_point (btn_start,sld_nipOnPath));

btn_end = uicontrol('Style', 'pushbutton', 'String', 'Set finish','Units','normalized'...
    ,'Position', [0.8 0.17 0.1 0.02],...
    'Callback', @(btn_end,event) set_end_point (btn_end,sld_nipOnPath));





% popup = uicontrol('Style', 'popup',...
%            'String', {'parula','jet','hsv','hot','cool','gray'},...
%            'Position', [20 340 100 50],...
%            'Callback', @setmap);

%    Fig_handles.sld_nipOnPath=sld_nipOnPath;
%    guidata(gcbo,Fig_handles) ;

function Direction(popup_dir)

Fig_handles = guidata(gcbo);
Fig_handles.Dir=180*(get(popup_dir,'Value')-1);
guidata(gcbo,Fig_handles) ;
%          popup_dir=Fig_handles.popup_dir





function nipOnPath (sld_nipOnPath,fig_HO,h_Eulerian,...
    N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,R_cyl,...
    z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,...
    rQ,H_indentation,...
    rQ_Eulerian,Laser_Head_Checkbox)


Fig_handles = guidata(gcbo);
%delete previous laser head
fucked=findobj(fig_HO,'FaceColor','y');
delete(fucked);


%%
% Direction of old axis system, nx,ny,nz
% old_axis=Fig_handles.Tape_Roller_Handle{2};
% X_axis_old=[old_axis.XData(2)-old_axis.XData(1)...
%     old_axis.YData(2)-old_axis.YData(1)...
%     old_axis.ZData(2)-old_axis.ZData(1)];
% Y_axis_old=[old_axis.XData(5)-old_axis.XData(4)...
%     old_axis.YData(5)-old_axis.YData(4)...
%     old_axis.ZData(5)-old_axis.ZData(4)];
% Z_axis_old=[old_axis.XData(8)-old_axis.XData(7)...
%     old_axis.YData(8)-old_axis.YData(7)...
%     old_axis.ZData(8)-old_axis.ZData(7)];

% old_axisXData=old_axis.XData;
% old_axisYData=old_axis.YData;
% old_axisZData=old_axis.ZData;



%%
%center
% Center_Roller=[old_axis.XData(1)...
%     old_axis.YData(1)...
%    old_axis.ZData(1)];



% XinParent=X_axis_old*Fig_handles.Rot_Roller_axis;
% YinParent=Y_axis_old*Fig_handles.Rot_Roller_axis;
% ZinParent=Z_axis_old*Fig_handles.Rot_Roller_axis;

% XinParent=XinParent/norm(XinParent);
% YinParent=YinParent/norm(YinParent);
% ZinParent=ZinParent/norm(ZinParent);
%%

for ii=1:length(Fig_handles.Tape_Roller_Handle)
    delete(Fig_handles.Tape_Roller_Handle{ii});
end
delete(Fig_handles.h_Eulerian);

% indicate 0 or 180 CW or CCW
Dir=Fig_handles.Dir;

Points_center_path=Fig_handles.Points_center_path;
Lag_points=Fig_handles.Lag_points;



index=floor(get(sld_nipOnPath,'Value'));
nip_point_M=Points_center_path(1:3,index);
Fig_handles.sld_index=index;

if index < length(Points_center_path)
    Vec_tan=Points_center_path(1:3,index)-Points_center_path(1:3,index+1);
    
else
    Vec_tan= Points_center_path(1:3,index)-Points_center_path(1:3,index-1);
end

tv=nip_point_M;

[Tape_Roller_Handle,Rot_Roller_axis]=Tape_plot_3D_ALL_objects_Dome_Car(N_tape,W_tape,R_tape,L_flat,tv,-abs(mod(th_y,90))+90,thick_T,deg_tape,W_R,theta_ind,R_cyl,z_cyl_end,H_indentation,Vec_tan,Dir);


% to know current handle
Fig_handles.Tape_Roller_Handle=Tape_Roller_Handle;

% to know current Rot roller axis
Fig_handles.Rot_Roller_axis=Rot_Roller_axis;



New_axis=Tape_Roller_Handle{2};
AxisC_New=[New_axis.XData(1)...
    New_axis.YData(1)...
    New_axis.ZData(1)];

if get(Laser_Head_Checkbox,'Value') ==1
    
    % check whether laser direction is exist as refference
    if ~isempty(Fig_handles.Laser_direction)
        
        delete(Fig_handles.h_laser_dir);
        Center_Roller=Fig_handles.Center_Roller_Ref;
        
        %Direction
        % delta=Fig_handles.Laser_direction;
        % laser Center
        pos=Fig_handles.Laser_Head_L_xyz0;
        %Laser Target
        Laser_Target=Fig_handles.Laser_Target;
        
        % Transformation
        
        %center
        Temp_pos=pos-Center_Roller;
        % back to Refference
        pos_in_parent=Temp_pos*Fig_handles.Rot_Roller_axis_Ref;
        % go to new system
        pos_new=(pos_in_parent/Rot_Roller_axis)+AxisC_New;
        
        
        Temp_Laser_Target=Laser_Target-Center_Roller;
        % back to Refference
        Laser_Target_in_parent=Temp_Laser_Target*Fig_handles.Rot_Roller_axis_Ref;
        % go to new system
        Laser_Target_new=(Laser_Target_in_parent/Rot_Roller_axis)+AxisC_New;
        
        % current laser direction
        % delta_in_parent=delta*Fig_handles.Rot_Roller_axis;
        delta_new=Laser_Target_new-pos_new;
        
        
        h_laser_dir=quiver3(pos_new(1),pos_new(2),pos_new(3),delta_new(1),delta_new(2),delta_new(3),0,'k','MaxHeadSize',0.5,'color','r','LineWidth',2);
        
        
        % set (text_orientation ,'String',sprintf(' nx=%f \n ny=%f \n nz=%f',delta_new(1),delta_new(2),delta_new(3)),'Fontsize',16);
        Fig_handles.h_laser_dir=h_laser_dir;
        
        
        
        %%
        % for rectangular laser head
        X_Actual=Fig_handles.X_Head;
        Y_Actual=Fig_handles.Y_Head;
        Z_Actual=Fig_handles.Z_Head;
        Power_Actual=0*Z_Actual;
        
        % Trnasform the laser head points to the parent axis and then to the
        % current axis
        
        P_head=zeros(3,4);
        % for ii=1:4
        %   P_head(1:3,ii)=[X_Actual(ii) Y_Actual(ii) Z_Actual(ii)] ;
        % end
        P_head(1:3,:)=[X_Actual(1:end); Y_Actual(1:end) ;Z_Actual(1:end)] ;
        
        
        %center
        Temp_P_head=P_head-Center_Roller';
        % back to Refference
        Temp_P_head_in_parent=Temp_P_head'*Fig_handles.Rot_Roller_axis_Ref;
        % go to new system
        P_head_new=((Temp_P_head_in_parent/Rot_Roller_axis)+AxisC_New)';
        
        X_Actual(1:4)=P_head_new(1,1:4);
        
        Y_Actual(1:4)=P_head_new(2,1:4);
        
        Z_Actual(1:4)=P_head_new(3,1:4);
        
        surf(X_Actual,Y_Actual,Z_Actual,Power_Actual,'LineStyle','-', 'FaceColor','y');
        
        
        
        
    end
    
    
    
    %  plot3(P1(1),P1(2),P1(3),'r*');
    %  plot3(P2(1),P2(2),P2(3),'r*');
    %  plot3(P3(1),P3(2),P3(3),'r*');
    
    
end





%%
% New_axis=Tape_Roller_Handle{2};
%
% X_axis_New=[New_axis.XData(2)-New_axis.XData(1)...
%     New_axis.YData(2)-New_axis.YData(1)...
%     New_axis.ZData(2)-New_axis.ZData(1)];
% Y_axis_New=[New_axis.XData(5)-New_axis.XData(4)...
%     New_axis.YData(5)-New_axis.YData(4)...
%     New_axis.ZData(5)-New_axis.ZData(4)];
% Z_axis_New=[New_axis.XData(8)-New_axis.XData(7)...
%     New_axis.YData(8)-New_axis.YData(7)...
%     New_axis.ZData(8)-New_axis.ZData(7)];


% X_axis_New=[New_axis.XData(2)-New_axis.XData(1);
% Y_axis_New=New_axis.YData(5)-New_axis.YData(4);
% Z_axis_New=New_axis.ZData(8)-New_axis.ZData(7);


% v1=X_axis_old'/norm(X_axis_old);   % the initial laser head on xy plane
% v2=X_axis_New'/norm(X_axis_New);
% Rot_plane_X=transpose(fcn_RotationFromTwoVectors(v1, v2)); % find the rotation matrix


% v1=[0; 1 ;0];   % the initial laser head on xy plane
% v2=Y_axis_old'/norm(Y_axis_old);
% Rot_plane_GL2L=transpose(fcn_RotationFromTwoVectors(v1, v2));
%
%
% v1=Y_axis_old'/norm(Y_axis_old);   % the initial laser head on xy plane
% v2=Y_axis_New'/norm(Y_axis_New);
% Rot_plane_Y=transpose(fcn_RotationFromTwoVectors(v1, v2));

% v1=Z_axis_old'/norm(Z_axis_old);   % the initial laser head on xy plane
% v2=Z_axis_New'/norm(Z_axis_New);
% Rot_plane_Z=transpose(fcn_RotationFromTwoVectors(v1, v2));


% Rot_Movement=Rot_plane_X*Rot_plane_Y*Rot_plane_Z;

% X_axis_old=[1 0 0]*1e-2;
% Y_axis_old=[0 1 0]*1e-2;
% Z_axis_old=[0 0 1]*1e-2;
%
%
% X_generated=(Rot_plane_Y*Rot_plane_GL2L\X_axis_old');
%
%  theta_2Refconfig=acosd((X_generated'*X_axis_New')/(norm(X_generated)*norm(X_axis_New)));
%
%  vector_cross=cross(X_generated,X_axis_New)/(norm(X_generated)*norm(X_axis_New));
%  theta_normal_cross_product=acosd((vector_cross*Y_axis_New')/(norm(vector_cross)*norm(Y_axis_New)));
%    th_y_corr=-sign(cosd(theta_normal_cross_product))* theta_2Refconfig;


% Make a Rotation that transform axis to final axis
% Rot_y=[cosd(th_y_corr) 0 sind(th_y_corr) ; 0 1 0; -sind(th_y_corr) 0 cosd(th_y_corr)];
%
% Rot_z=([cosd(-th_y_corr)  -sind(-th_y_corr) 0 ;  sind(-th_y_corr)  cosd(-th_y_corr) 0 ; 0 0 1]);
%
% Rot_x=([1 0 0 ;0 cosd(th_y_corr)  -sind(th_y_corr)  ; 0  sind(th_y_corr)  cosd(th_y_corr) ]);
%
%
% Rot_Movement=1*Rot_plane_Y*Rot_plane_GL2L ;




% Temp=X_axis_old'

% center of previous roller location
% P_0=[old_axisXData(1), old_axisYData(1),old_axisZData(1)]';
% % P_Y=[old_axisXData(5), old_axisYData(5),old_axisZData(5)]';
% % P_Z=[old_axisXData(8), old_axisYData(8),old_axisZData(8)]';
%
% % from last roller to the new roller
% Tr=[New_axis.XData(1)-old_axisXData(1),...
%     New_axis.YData(1)-old_axisYData(1),...
%     New_axis.ZData(1)-old_axisZData(1)]';
%
% P1=(Rot_Movement\X_axis_old')+Tr+P_0;
% P2=(Rot_Movement\Y_axis_old')+Tr+P_0;
% P3=(Rot_Movement\Z_axis_old')+Tr+P_0;

% plot3(P1(1),P1(2),P1(3),'r*');
% plot3(P2(1),P2(2),P2(3),'r*');
% plot3(P3(1),P3(2),P3(3),'r*');

% Y_axis_old=old_axis.YData(5)-old_axis.YData(4);
% Z_axis_old=old_axis.ZData(8)-old_axis.ZData(7);

%%




% Fig_handles.nip_point_M=nip_point_M;


Long_CV=str2double(get(rQ_Eulerian,'string'));


Fig_handles.Long_CV=Long_CV;

% Long_CV=20;

xnode=2*No_dev -1;
ynode=Long_CV ; %No_dev_L-starting_index;

X=zeros(ynode,xnode);
Y=zeros(ynode,xnode);
Z=zeros(ynode,xnode);
Intensity=Z;

% index+(ynode *cosd(Dir)) >1

if cosd(Dir) >0
    if index+ynode > length(Points_center_path)
        index=length(Points_center_path)-ynode;
    end
else
    if index-ynode < 1
        index=ynode+1;
    end
end


for kk=1:xnode
    %     index=(((1:ynode)*xnode)-kk+1);
    X(1:end,kk)=Lag_points(1,index+(1:ynode)*cosd(Dir),kk);
    Y(1:end,kk)=Lag_points(2,index+(1:ynode)*cosd(Dir),kk);
    Z(1:end,kk)=Lag_points(3,index+(1:ynode)*cosd(Dir),kk);
    Intensity(1:end,kk)=0; % RHS in the thermal model
end
% pause(0.5);

h_Eulerian=mesh(X,Y,Z,Intensity,'FaceColor','None','EdgeColor','y');
%%


Fig_handles.h_Eulerian=h_Eulerian;
guidata(gcbo,Fig_handles) ;





% make laser surface for rotation
function  Laser_surf(text_Rot,sld,fig_HO,Ax,Ay,rQ_Lasertarget,rQ_Source)

%delete previous laser head
fucked=findobj(fig_HO,'FaceColor','y');
delete(fucked);


if str2num(rQ_Lasertarget.String) == str2num(rQ_Source.String)

    dcm_obj = datacursormode(fig_HO);
c_info = getCursorInfo(dcm_obj);
    Laser_Target=c_info(2).Position;
pos=c_info(1).Position;

else
    Laser_Target=str2num(rQ_Lasertarget.String);
    pos=str2num(rQ_Source.String);

end



L_xyz0=pos;

delta=Laser_Target-pos;

v1=[0 0 1]';   % the initial laser head on xy plane
v2=(delta/norm(delta))';
Rot_plane=transpose(fcn_RotationFromTwoVectors(v1, v2)); % find the rotation matrix





rotation_deg=get(sld,'value');
set(text_Rot,'String',sprintf('Laser Head Rotation=%d?',rotation_deg));

rotation=rotation_deg*(pi/180);


xp_laser=linspace(-Ax,Ax,2);
yp_laser=linspace(-Ay,Ay,2);
[X,Y]=meshgrid(xp_laser,yp_laser);  % ny * nx = size
Z=zeros(size(X));
Power_Actual=Z;

X_Actual=zeros(size(X));
Y_Actual=zeros(size(X));
Z_Actual=zeros(size(X));

counter=0;
for ii=1:length(xp_laser)
    for jj=1:length(yp_laser)
        counter=counter+1;
        
        new_points=Transformation_Rot_z([0 0 0],-rotation,[xp_laser(ii);yp_laser(jj); 0]);
        laser_point_Actual= [new_points(1),new_points(2), new_points(3)]* (Rot_plane);
        X_Actual(jj,ii) =laser_point_Actual(1);
        Y_Actual(jj,ii) =laser_point_Actual(2);
        Z_Actual(jj,ii) =laser_point_Actual(3);
        %         laser_point_Actual_all(counter,1:3)=laser_point_Actual;
        
    end
end

X_Actual=X_Actual+L_xyz0(1);
Y_Actual=Y_Actual+L_xyz0(2);
Z_Actual=Z_Actual+L_xyz0(3);

Fig_handles = guidata(gcbo);

Fig_handles.X_Head=X_Actual;
Fig_handles.Y_Head=Y_Actual;
Fig_handles.Z_Head=Z_Actual;


Fig_handles.Head_rotation=rotation;


guidata(gcbo,Fig_handles) ;

surf(X_Actual,Y_Actual,Z_Actual,Power_Actual,'LineStyle','-', 'FaceColor','y');

% surf(X,Y,Z)



%% make target

function  circle_Laser=Tape_target(btn3,rQ,th_y,fig_HO,circle_Laser,rQ_Lasertarget,rQ_Source)





if str2num(rQ_Lasertarget.String) == str2num(rQ_Source.String)

    dcm_obj = datacursormode(fig_HO);
c_info = getCursorInfo(dcm_obj);
    Laser_Target=c_info(2).Position;
pos=c_info(1).Position;

else
    Laser_Target=str2num(rQ_Lasertarget.String);
    pos=str2num(rQ_Source.String);

end






dcm_obj = datacursormode(fig_HO);
c_info = getCursorInfo(dcm_obj);
% Laser_Target=c_info(2).Position;


% fucked=findobj(fig_HO,'color','g','LineStyle','--');
fucked=findobj(fig_HO,'FaceColor','g');
delete(fucked);

distance=str2double(get(rQ,'String'));




[x,y,z]=sphere(50);
x=x*distance + (Laser_Target(1));
y=y*distance + (Laser_Target(2));
z=z*distance + (Laser_Target(3));

circle_Laser=mesh(x,y,z,'FaceAlpha',0.1,'EdgeAlpha',0.1,'FaceColor','g');

%%


dcm_obj = datacursormode(fig_HO);
% set(dcm_obj,'DisplayStyle','datatip',...
%     'SnapToDataVertex','off','Enable','on');


% set(dcm_obj,'UpdateFcn',@NewCallback_datacursor);

hTarget = handle(circle_Laser);
hDatatip = dcm_obj.createDatatip(hTarget);
set(hDatatip,'Host',hTarget);
datacursormode toggle;







function  Accept( btn1,rQ,th_y,fig_HO,Laser_Target)



% fucked=findobj(fig_HO,'color','g','LineStyle','--');
fucked=findobj(fig_HO,'FaceColor','g');

delete(fucked);


distance=str2double(get(rQ,'String'));



[x,y,z]=sphere(50);
x=x*distance + (Laser_Target(1));
y=y*distance + (Laser_Target(2));
z=z*distance + (Laser_Target(3));

circle_Laser=mesh(x,y,z,'FaceAlpha',0.1,'EdgeAlpha',0.1,'FaceColor','g');
%%


dcm_obj = datacursormode(fig_HO);
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on');


set(dcm_obj,'UpdateFcn',@NewCallback_datacursor);

hTarget = handle(circle_Laser);
hDatatip = dcm_obj.createDatatip(hTarget);
set(hDatatip,'Host',hTarget);
datacursormode toggle;



function Put ( btn2,fig_HO,text_orientation,rQ_Lasertarget,rQ_Source)


Fig_handles = guidata(gcbo);

old_axis=Fig_handles.Tape_Roller_Handle{2};

%center of Ref
Center_Roller=[old_axis.XData(1)...
    old_axis.YData(1)...
    old_axis.ZData(1)];


% save center_Ref
Fig_handles.Center_Roller_Ref=Center_Roller;

% Rot Roller axis
Fig_handles.Rot_Roller_axis_Ref=Fig_handles.Rot_Roller_axis;








% fucked=findobj(fig_HO,'color','r','LineWidth',2);
% delete(fucked);
delete(Fig_handles.h_laser_dir);



dcm_obj = datacursormode(fig_HO);
c_info = getCursorInfo(dcm_obj);





if str2num(rQ_Lasertarget.String) == str2num(rQ_Source.String)

    % dcm_obj = datacursormode(fig_HO);
% c_info = getCursorInfo(dcm_obj);

    Laser_Target=c_info(2).Position;

pos=c_info(1).Position;


else


    Laser_Target=str2num(rQ_Lasertarget.String);
    pos=str2num(rQ_Source.String);

end




% Laser_Target=c_info(2).Position;

% pos=c_info(1).Position;

delta=Laser_Target-pos;

h_laser_dir=quiver3(pos(1),pos(2),pos(3),delta(1),delta(2),delta(3),0,'k','MaxHeadSize',0.5,'color','r','LineWidth',2);

norm_delta=norm(delta);
delta=delta/norm_delta;
set (text_orientation ,'String',sprintf(' nx=%f \n ny=%f \n nz=%f',delta(1),delta(2),delta(3)),'Fontsize',12);

%Direction
Fig_handles.Laser_direction=delta;
% laser Center
Fig_handles.Laser_Head_L_xyz0=pos;
%Laser Target
Fig_handles.Laser_Target=Laser_Target;

%   Fig_handles.nip_point_M
Fig_handles.h_laser_dir=h_laser_dir;

guidata(gcbo,Fig_handles) ;





%%



function  Select_Nip_point( btn4,fig_HO)


%  Fig_handles.h_path=h_path;

% fucked=findobj(fig_HO,'color','g','LineStyle','--');
fucked=findobj(fig_HO,'FaceColor','g');


delete(fucked);

fucked2=findobj(fig_HO,'color','r','LineWidth',2);
delete(fucked2);



Fig_handles = guidata(gcbo);
% Modify the value of your counter
% myhandles.numberOfErrors = myhandles.numberOfErrors + 1;
% Save the change you made to the structure


Tape_Roller_Handle=Fig_handles.Tape_Roller_Handle;

h_Eulerian=Fig_handles.h_Eulerian;

% Save the structure
% guidata(fig_HO,Fig_handles) ;

delete(h_Eulerian);
delete(Fig_handles.h_path);

for ii=1:length(Tape_Roller_Handle)
    
    delete(Tape_Roller_Handle{ii});
    
end



dcm_obj = datacursormode(fig_HO);
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on');
datacursormode on;








function  Accept_Nip_point( btn5,fig_HO,...
    Tape_Roller_Handle,...
    N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,R_cyl,...
    z_cyl_end,L_prim,w,sui,No_dev,Nip_Mov,...
    rQ,H_indentation,...
    rQ_dis_path,rQ_Eulerian,sld_nipOnPath,rQ_nip_point)




Fig_handles = guidata(gcbo);




No_dev_L=str2double(get(rQ_dis_path,'string'));

Long_CV=str2double(get(rQ_Eulerian,'string'));




dcm_obj = datacursormode(fig_HO);
c_info = getCursorInfo(dcm_obj);


if rQ_nip_point.String(1)=='I'
% center of the nip-point
nip_point_M=c_info.Position;

else
   nip_point_M= str2num(rQ_nip_point.String);
    
end


[h_path,Points_center_path,Lag_points, Boarders,starting_index]=Winding_path_cartesian_adaptive (R_cyl,z_cyl_end,L_prim,w,-abs(mod(th_y,90))-90,No_dev_L,No_dev,Nip_Mov,nip_point_M);



%       sld_nipOnPath=Fig_handles.sld_nipOnPath;
%    guidata(gcbo,Fig_handles) ;

max_point=length(Points_center_path);

set(sld_nipOnPath,'SliderStep', [1/max_point , 10/max_point ],...
    'Min',1,'Max',max_point,'Value',1);





% save in structure to be accessible from another function
Fig_handles.h_path=h_path;
Fig_handles.Points_center_path=Points_center_path;
Fig_handles.Lag_points=Lag_points;



index=1;

% tangent along the path
if index < length(Points_center_path)
    Vec_tan=Points_center_path(1:3,index)-Points_center_path(1:3,index+1);
    
else
    Vec_tan= Points_center_path(1:3,index)-Points_center_path(1:3,index-1);
end

Dir=Fig_handles.Dir;
% start from begining of the path
tv=Points_center_path(:,1);
[Tape_Roller_Handle,Rot_Roller_axis]=Tape_plot_3D_ALL_objects_Dome_Car(N_tape,W_tape,R_tape,L_flat,tv,-abs(mod(th_y,90))+90,thick_T,deg_tape,W_R,theta_ind,R_cyl,z_cyl_end,H_indentation,Vec_tan,Dir);

Fig_handles.Tape_Roller_Handle=Tape_Roller_Handle;
Fig_handles.Rot_Roller_axis=Rot_Roller_axis;


% plot Substrate






% Long_CV=20;

xnode=2*No_dev -1;
ynode=Long_CV ; %No_dev_L-starting_index;

X=zeros(ynode,xnode);
Y=zeros(ynode,xnode);
Z=zeros(ynode,xnode);
Intensity=Z;


% to know the direction
if cosd(Dir) >0
    if index+ynode > length(Points_center_path)
        index=length(Points_center_path)-ynode;
    end
else
    if index-ynode < 1
        index=ynode+1;
    end
end


% filling XYZ matrix with Eulerian points
for kk=1:xnode
    %     index=(((1:ynode)*xnode)-kk+1);
    X(1:end,kk)=Lag_points(1,index+(1:ynode)*cosd(Dir),kk);
    Y(1:end,kk)=Lag_points(2,index+(1:ynode)*cosd(Dir),kk);
    Z(1:end,kk)=Lag_points(3,index+(1:ynode)*cosd(Dir),kk);
    %     Intensity(1:end,kk)=0; % RHS in the thermal model
end
% pause(0.5);

h_Eulerian=mesh(X,Y,Z,Intensity,'FaceColor','None','EdgeColor','y');



Fig_handles.h_Eulerian=h_Eulerian;

% Save the structure


Fig_handles.nip_point_M=nip_point_M;

guidata(gcbo,Fig_handles) ;




%%

distance=str2double(get(rQ,'String'));

% r=distance;

% tt=linspace(0,2*pi,50);
% xx=r*sin(tt);
% yy=r*cos(tt);
% zz=yy*0;
%
% Circle_Rotation=-(pi/180)*th_y ;
%
%
% new_xyz=Transformation_Rot_y(nip_point_M,Circle_Rotation,[xx;yy;zz]);
%
%
% circle_Laser=plot3(new_xyz(1,:),new_xyz(2,:),new_xyz(3,:),'g--');


[x,y,z]=sphere(50);
x=x*distance + (nip_point_M(1));
y=y*distance + (nip_point_M(2));
z=z*distance + (nip_point_M(3));

circle_Laser=mesh(x,y,z,'FaceAlpha',0.1,'EdgeAlpha',0.1,'FaceColor','g');
%%


dcm_obj = datacursormode(fig_HO);
% set(dcm_obj,'DisplayStyle','datatip',...
%     'SnapToDataVertex','off','Enable','on');


% set(dcm_obj,'UpdateFcn',@NewCallback_datacursor);

hTarget = handle(circle_Laser);
hDatatip = dcm_obj.createDatatip(hTarget);
set(hDatatip,'Host',hTarget);
datacursormode toggle;




function set_end_point(btn_end,sld_nipOnPath)

Fig_handles = guidata(gcbo);

if exist('Fig_handles.hplot_finish')
delete(Fig_handles.hplot_finish);
end



index=floor(get(sld_nipOnPath,'Value'));
Fig_handles.finish_index=index;
Points_center_path=Fig_handles.Points_center_path;

nip_point_M=Points_center_path(1:3,index);

% Fig_handles.hplot_finish=plot3(nip_point_M(1),nip_point_M(2),nip_point_M(3),'bo','MarkerSize',10);

Fig_handles.hplot_finish=text(nip_point_M(1),nip_point_M(2),nip_point_M(3),'\leftarrow \fontsize{15} Finish location');




set(btn_end,'String',sprintf('Set finish_i=%d',index));

% Fig_handles.Tape_Roller_Handle{ii}
guidata(gcbo,Fig_handles);

function set_start_point(btn_start,sld_nipOnPath)

Fig_handles = guidata(gcbo);

if exist('Fig_handles.hplot_start')
delete(Fig_handles.hplot_start);
end

index=floor(get(sld_nipOnPath,'Value'));
Fig_handles.start_index=index;
Points_center_path=Fig_handles.Points_center_path;

nip_point_M=Points_center_path(1:3,index);

% Fig_handles.hplot_start=plot3(nip_point_M(1),nip_point_M(2),nip_point_M(3),'bo','MarkerSize',10);

Fig_handles.hplot_start=text(nip_point_M(1),nip_point_M(2),nip_point_M(3),'\leftarrow \fontsize{15} Start location');



set(btn_start,'String',sprintf('Set start_i=%d',index));

% Fig_handles.Tape_Roller_Handle{ii}
guidata(gcbo,Fig_handles);





function  Save_intxt( btn_save_path_info, th_y,R_cyl,...
    z_cyl_end,No_dev,Ax,Ay,nx,ny)


%Data to be written
%     Fig_handles.finish_index==Fig_handles.start_index  >> Steady-state
%     condition
% else is the unsteady optical model  >> UOT
% Assign start point and end point
Fig_handles = guidata(gcbo);





if ~isempty(Fig_handles.finish_index) & ~isempty(Fig_handles.start_index)
    
    % if the range is determined, then try to save data of rays, locations
    % during the path
    
    
    prompt={'Enter a job analysis name'};
name = 'Analysis Name';
defaultans = {'Example0'};
options.Interpreter = 'tex';
jobname = inputdlg(prompt,name,[1 40],defaultans,options);
    
    
%  jobname='Example0';
dir=strcat('.\Analysis_UOT\', jobname{:});
mkdir(dir);   
    
    
    
    
    
    
    
    
    
    %%
    % laser head in refference Config
            Center_Roller=Fig_handles.Center_Roller_Ref;
        
        %Direction
        % delta=Fig_handles.Laser_direction;
        % laser Center
        pos=Fig_handles.Laser_Head_L_xyz0;
        %Laser Target
        Laser_Target=Fig_handles.Laser_Target;
        
        
    %delete previous laser head



L_xyz0=Fig_handles.Laser_Head_L_xyz0;

delta=Fig_handles.Laser_direction;

v1=[0 0 1]';   % the initial laser head on xy plane
v2=(delta/norm(delta))';
Rot_plane=transpose(fcn_RotationFromTwoVectors(v1, v2)); % find the rotation matrix


rotation=Fig_handles.Head_rotation;



xp_laser=linspace(-Ax,Ax,nx);
yp_laser=linspace(-Ay,Ay,ny);
[X,Y]=meshgrid(xp_laser,yp_laser);  % ny * nx = size
% Z=zeros(size(X));
 Power_Actual=X;

X_Actual=zeros(size(X));
Y_Actual=zeros(size(X));
Z_Actual=zeros(size(X));



counter=0;
for ii=1:length(xp_laser)
    for jj=1:length(yp_laser)
        counter=counter+1;
        
        new_points=Transformation_Rot_z([0 0 0],-rotation,[xp_laser(ii);yp_laser(jj); 0]);
        laser_point_Actual= [new_points(1),new_points(2), new_points(3)]* (Rot_plane);
        X_Actual(jj,ii) =laser_point_Actual(1);
        Y_Actual(jj,ii) =laser_point_Actual(2);
        Z_Actual(jj,ii) =laser_point_Actual(3);
  
    end
end


% laser head surface in a refference config!
X_Actual=X_Actual+L_xyz0(1);
Y_Actual=Y_Actual+L_xyz0(2);
Z_Actual=Z_Actual+L_xyz0(3);


% change to parent config

        P_head=zeros(3,nx*ny);       
        P_head(1:3,:)=[X_Actual(1:end); Y_Actual(1:end) ;Z_Actual(1:end)] ;
                %center
        Temp_P_head=P_head-Center_Roller';
        % back to Refference
        Temp_P_head_in_parent=Temp_P_head'*Fig_handles.Rot_Roller_axis_Ref;

    
    
    % indicate 0 or 180 CW or CCW
    Dir=Fig_handles.Dir;
    
    Points_center_path=Fig_handles.Points_center_path;
    Lag_points=Fig_handles.Lag_points;
    
            Long_CV=Fig_handles.Long_CV;
        
        % Long_CV=20;
        
        xnode=2*No_dev -1;
        ynode=Long_CV ; %No_dev_L-starting_index;
        
        X=zeros(ynode,xnode);
        Y=zeros(ynode,xnode);
        Z=zeros(ynode,xnode);
        Intensity=Z;
    
    
    path_No=abs(Fig_handles.start_index-Fig_handles.finish_index);
    
    Rot_Roller_axis_all=zeros(path_No,3,3);
    
      pos_new_all=zeros(path_No,1,3);
        delta_new_all=zeros(path_No,1,3);
        nip_point_M_all=zeros(path_No,1,3);
        
        
      CV_mesh=zeros(path_No,3,ynode,xnode);
    LaserH_mesh=zeros(path_No,3,ny,nx);
    
    kickoff_dir=sign(Fig_handles.finish_index-Fig_handles.start_index);
    
    counter=0;
    for ii=Fig_handles.start_index:kickoff_dir:Fig_handles.finish_index
        counter=counter+1;
        index=ii ;
        nip_point_M=Points_center_path(1:3,index);
        % Fig_handles.sld_index=index;
        
        if index < length(Points_center_path)
            Vec_tan=Points_center_path(1:3,index)-Points_center_path(1:3,index+1);
            
        else
            Vec_tan= Points_center_path(1:3,index)-Points_center_path(1:3,index-1);
        end
        
        tv=nip_point_M;
        
        %to save all nip-pointto put roller-tape system
        nip_point_M_all(counter,:,:)=nip_point_M';
        % [Tape_Roller_Handle,Rot_Roller_axis]=Tape_plot_3D_ALL_objects_Dome_Car(N_tape,W_tape,R_tape,L_flat,tv,-abs(mod(th_y,90))+90,thick_T,deg_tape,W_R,theta_ind,R_cyl,z_cyl_end,H_indentation,Vec_tan,Dir);
        
        [Rot_Roller_axis,~]=Rot_Matrix_Finder_RollerTape_dome_Car(tv,-abs(mod(th_y,90))+90,R_cyl,z_cyl_end,Vec_tan',Dir);
        
  Rot_Roller_axis_all(counter,:,:)=Rot_Roller_axis;
        
       
        AxisC_New=tv';
        
        % it should be adaptive and know
        % if get(Laser_Head_Checkbox,'Value') ==1
        
        % check whether laser direction is exist as refference
        %     if ~isempty(Fig_handles.Laser_direction)
        
%         delete(Fig_handles.h_laser_dir);

        
        % Transformation
        
        %center
        Temp_pos=pos-Center_Roller;
        % back to Refference
        pos_in_parent=Temp_pos*Fig_handles.Rot_Roller_axis_Ref;
        % go to new system
        pos_new=(pos_in_parent/Rot_Roller_axis)+AxisC_New;
        
        %save laser head center
        pos_new_all(counter,:,:)=pos_new;
      
        
        Temp_Laser_Target=Laser_Target-Center_Roller;
        % back to Refference
        Laser_Target_in_parent=Temp_Laser_Target*Fig_handles.Rot_Roller_axis_Ref;
        % go to new system
        Laser_Target_new=(Laser_Target_in_parent/Rot_Roller_axis)+AxisC_New;
        
        % current laser direction
        % delta_in_parent=delta*Fig_handles.Rot_Roller_axis;
        delta_new=Laser_Target_new-pos_new;
        
          %save laser direction
         delta_new_all(counter,:,:)=delta_new;
        
        h_laser_dir=quiver3(pos_new(1),pos_new(2),pos_new(3),delta_new(1),delta_new(2),delta_new(3),0,'k','MaxHeadSize',0.5,'color','r','LineWidth',1.5);
        
        
        % set (text_orientation ,'String',sprintf(' nx=%f \n ny=%f \n nz=%f',delta_new(1),delta_new(2),delta_new(3)),'Fontsize',16);
%         Fig_handles.h_laser_dir=h_laser_dir;
        
        
        
        %%
        % for rectangular laser head
%         X_Actual=Fig_handles.X_Head;
%         Y_Actual=Fig_handles.Y_Head;
%         Z_Actual=Fig_handles.Z_Head;
%         
      
        
%         Power_Actual=0*Z_Actual;
        

        % go to new system
        P_head_new=((Temp_P_head_in_parent/Rot_Roller_axis)+AxisC_New)';
        
        
        % save info of laser head locations
        X_Actual(:)=P_head_new(1,:);
        
        Y_Actual(:)=P_head_new(2,:);
        
        Z_Actual(:)=P_head_new(3,:);
        
      laser_head=surf(X_Actual,Y_Actual,Z_Actual,Power_Actual,'LineStyle','-', 'FaceColor','y');
        
        
        % 4D matrix to save all data
         LaserH_mesh(counter,1,:,:)=X_Actual;
         LaserH_mesh(counter,2,:,:)=Y_Actual;
          LaserH_mesh(counter,3,:,:)=Z_Actual;

        
        % index+(ynode *cosd(Dir)) >1
        
        if cosd(Dir) >0
            if index+ynode > length(Points_center_path)
                index=length(Points_center_path)-ynode;
            end
        else
            if index-ynode < 1
                index=ynode+1;
            end
        end
        
        
        for kk=1:xnode
            %     index=(((1:ynode)*xnode)-kk+1);
            X(1:end,kk)=Lag_points(1,index+(1:ynode)*cosd(Dir),kk);
            Y(1:end,kk)=Lag_points(2,index+(1:ynode)*cosd(Dir),kk);
            Z(1:end,kk)=Lag_points(3,index+(1:ynode)*cosd(Dir),kk);
            Intensity(1:end,kk)=0; % RHS in the thermal model
        end
        
      
        
        
        CV_mesh(counter,1,:,:)=X;
         CV_mesh(counter,2,:,:)=Y;
          CV_mesh(counter,3,:,:)=Z;
          
        % pause(0.5);
        
        % h_Eulerian=mesh(X,Y,Z,Intensity,'FaceColor','None','EdgeColor','y');
        %%
        
        
        % Fig_handles.h_Eulerian=h_Eulerian;
       
        pause(0.01);
          delete(h_laser_dir);
          delete(laser_head);
    end
    
    
    % save data as text file or matlab file
    
     guidata(gcbo,Fig_handles) ;
     
    Laser_head=zeros(1,4);
    
Laser_head(1)=Ax;
Laser_head(2)=Ay;
Laser_head(3)=nx;
Laser_head(4)=ny;
     
     save(strcat(dir,'\Kinematic_UOT.mat'),'Rot_Roller_axis_all','nip_point_M_all',...
         'CV_mesh','LaserH_mesh','delta_new_all','pos_new_all',...
         'jobname','Laser_head','R_cyl','z_cyl_end');
%   
% jobname='Example0';
% 
% mkdir(strcat('Analysis_UOT\', jobname));



     
     
else
    
    errordlg('Not enough inputs!');
    
    
end



function  Save_editbox (btn_save_editbox,rQ_Lasertarget,rQ_Source,rQ_nip_point,rQ)

    prompt={'Enter a  name'};
name = 'Analysis Name';
defaultans = {'EditBox_1'};
options.Interpreter = 'tex';
jobname = inputdlg(prompt,name,[1 40],defaultans,options);
    
    jobname=cell2mat(jobname);
%  jobname='Example0';
dir=strcat('.\Analysis_UOT\');
 mkdir(dir);  

target= str2num(rQ_Lasertarget.String);
source=str2num(rQ_Source.String);
nip= str2num(rQ_nip_point.String);
distance= str2num(rQ.String);
     
     save(strcat(dir,jobname,'.mat'),'target','source','nip','distance');
%   


function  Load_editbox (btn_save_editbox,rQ_Lasertarget,rQ_Source,rQ_nip_point,rQ)


   [file,path] = uigetfile('.\EditBox_1.mat','Open Editbox data');
  str=strcat(path,file);

load(str);

 set(rQ_Lasertarget,'String',num2str(target));
set(rQ_Source,'String',num2str(source));
 set(rQ_nip_point,'String',num2str(nip));
 set(rQ,'String',num2str(distance));



% Fig_handles.Laser_direction
%  Fig_handles.Laser_Head_L_xyz0
% Fig_handles.nip_point_M=nip_point_M;

%
% fid10 = fopen('.\Geo_Parameter.txt','r');
% out = textscan(fid10,'%s','delimiter',',');
% fclose (fid10);
% out_copy=out;
% % find index to replace
% Rxyz_ind = contains(out{1},'Rxyz');
% %   out{1}{Rxyz_ind};
% L_xyz0_ind = contains(out{1},'L_xyz0');
% %    out{1}{L_xyz0_ind};
% tv_ind = contains(out{1},'Roller_Pos_TV');
% Laser_head_Rot = contains(out{1},'Laser_head_Rot');
% %       out{1}{tv_ind};
% % remove first index to adjust writing new data
% out{1}(1)=[];
% Fig_handles = guidata(gcbo);
%
% out{1}{Rxyz_ind}= num2str(Fig_handles.Laser_direction);
% out{1}{L_xyz0_ind}= num2str(Fig_handles.Laser_Head_L_xyz0);
% out{1}{tv_ind}=num2str(Fig_handles.nip_point_M);
% rotation_deg=get(sld,'value');
% %         set(text_Rot,'String',sprintf('Laser Head Rotation=%d?',rotation_deg));
% out{1}{Laser_head_Rot}=num2str(rotation_deg);
%
% fid11 = fopen('.\Geo_Parameter.txt','w');
% % >> regexprep  % to replace in txt file
%
% for ii=1:2:length(out{1})  % number of lines
%
%     fprintf(fid11,'%s\n', strcat(out_copy{1}{ii},', ' ,out{1}{ii}) );
%
%
% end


