
% This program is able show 3D objects

function counter_ray=optical_3D_objects_Free_on_Sub(th_y,Tape_Sp,...
    R_cyl,z_cyl_end,Roller_Pos_TV,W_R,H_indentation,...
    L_prim,w,sui,Laser_head)

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


% Create figure
fig_HO=figure('Name','Laser HO finder','NumberTitle','off');
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


Mandrel_plot=plot3D_cylinder4optical_3D_objects(R_cyl,z_cyl_end);


%%
dcm_obj = datacursormode(fig_HO);
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on');


set(dcm_obj,'UpdateFcn',@NewCallback_datacursor);
hTarget = handle(Mandrel_plot(1));
hDatatip = dcm_obj.createDatatip(hTarget);
set(hDatatip,'Host',hTarget);
datacursormode toggle;


% 
% hCursorbar = graphics.cursorbar(fig_HO); 
% % drawnow
% hCursorbar.TargetMarkerStyle = 'o';


% c_info = getCursorInfo(dcm_obj);

% 
% nip_point_M=c_info.Position;

% 
% btn1 = uicontrol('Style', 'pushbutton', 'String', 'Nip_point Location','Units','normalized'...
%     ,'Position', [0.1 0.9 0.2 0.8],...
%     'Callback', @(btn1,event) Accept_Nip_point (btn1,h_T,Tape_Roller_Handle,fig_HO) );



%%




 theta_ind=acosd((R_tape-H_indentation)/R_tape);  %theta_ind in degree
 
 Nip_Mov =R_tape*sind(theta_ind); 
 
 %%
 No_dev=20;
 No_dev_L=ceil(abs(L_prim*100)); % for representation
 
%    dcm_obj = datacursormode(fig_HO);
% c_info = getCursorInfo(dcm_obj);
% % center of the nip-point
% nip_point_M=c_info.Position;



 nip_point_M=[0,R_cyl(1),0];
%   if  R_cyl(1) automatically goes to the center of placement
 

 
%   nip_point_M=Roller_Pos_TV;
 

if R_cyl(1) ~= 0
 
 if nip_point_M (3) >=0 && nip_point_M (3) <=z_cyl_end
 
[points, Boarders,starting_index]=helical_3D_points_4free_onSub(R_cyl(1),z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,nip_point_M);


 elseif nip_point_M (3) < 0  
     % bottom dome part

[points, Boarders,starting_index]=sub_on_ellips(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,nip_point_M);

 elseif nip_point_M (3) > z_cyl_end
     
     % Upper Dome part
[points, Boarders,starting_index]=sub_on_ellips(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,nip_point_M);


 end
 
else
  
    [points, Boarders,starting_index]=sub_on_flat(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,nip_point_M);

    
end



xnode=2*No_dev -1;
ynode=No_dev_L-starting_index;

X=zeros(xnode,ynode);
Y=zeros(xnode,ynode);
Z=zeros(xnode,ynode);
Intensity=Z;

for kk=1:xnode
    index=(((1:ynode)*xnode)-kk+1);
    X(kk,1:end)=points(1,index);
    Y(kk,1:end)=points(2,index);
    Z(kk,1:end)=points(3,index);
     Intensity(kk,1:end)=0; % RHS in the thermal model
end
% pause(0.5);

h_T=surf(X,Y,Z,Intensity,'LineStyle','None');
%%


Fig_handles.Substrate_Handle=h_T;
% Fig_handles.numberOfErrors = 0; 
% Save the structure
guidata(fig_HO,Fig_handles) ;


tv=nip_point_M;  % for the tape
% tv_R=[0;R_cyl(1)+R_tape+thick_T-H_indentation;tv3];  % for center of Roller



[Tape_Roller_Handle]=Tape_plot_3D_ALL_objects_Free_on_Sub(N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,R_cyl,z_cyl_end,H_indentation);


Fig_handles.Tape_Roller_Handle=Tape_Roller_Handle;
% Fig_handles.numberOfErrors = 0; 
% Save the structure
guidata(fig_HO,Fig_handles) ;



camlight(-90,90);

dcm_obj = datacursormode(fig_HO);
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on');


set(dcm_obj,'UpdateFcn',@NewCallback_datacursor);
hTarget = handle(h_T);
hDatatip = dcm_obj.createDatatip(hTarget);
set(hDatatip,'Host',hTarget);
datacursormode toggle;



c_info = getCursorInfo(dcm_obj);


Laser_Target=c_info.Position;





 Fig_handles.Laser_direction=[0 0 0];
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
    'String',sprintf(' nx=%f \n ny=%f \n nz=%f',0 ,0,0),'Fontsize',20,...
    'Position',[10 300 200 100],'Visible','on','Fontsize',20);



text1 = uicontrol('Style','text',...
    'String','Laser Distance',...
    'Position',[80 80 120 40],'Visible','on');

rQ = uicontrol('Style','edit',...
    'String',num2str(R_cyl(1)),...
    'Position',[80 55 120 40],'Visible','on');



btn1 = uicontrol('Style', 'pushbutton', 'String', 'OKe',...
    'Position', [80 10 100 40],...
    'Callback', @(btn1,event) Accept (btn1,rQ,th_y,fig_HO,Laser_Target));

btn2 = uicontrol('Style', 'pushbutton', 'String', 'Put Laser',...
    'Position', [300 10 100 40],...
    'Callback', @(btn2,event) Put (btn2,fig_HO,text_orientation));

btn3 = uicontrol('Style', 'pushbutton', 'String', 'Head Align Target',...
    'Position', [500 10 150 40],...
    'Callback', @(btn3,event) Tape_target (btn3,rQ,th_y,fig_HO,circle_Laser));


btn4 = uicontrol('Style', 'pushbutton', 'String', 'Nip_point Location','Units','normalized'...
    ,'Position', [0.8 0.8 0.2 0.05],...
    'Callback', @(btn4,event) Select_Nip_point (btn4,fig_HO) );


btn5 = uicontrol('Style', 'pushbutton', 'String', 'Accept Nip_point','Units','normalized'...
    ,'Position', [0.8 0.7 0.1 0.05],...
    'Callback', @(btn5,event) Accept_Nip_point (btn5,fig_HO,...
         Fig_handles.Tape_Roller_Handle,...
       N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,R_cyl,...
       z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,...
        rQ,H_indentation) );



    text_Rot = uicontrol('Style','text',...
    'String','Laser Head Rotation',...
    'Position',[40 240 180 30],'Visible','on');


 sld = uicontrol('Style', 'slider',...
        'Min',0,'Max',180,'Value',90,...
        'Position', [40 200 180 20],...
        'Callback', @(sld,event) Laser_surf (text_Rot,sld,fig_HO,Ax,Ay) ); 
    
btn6 = uicontrol('Style', 'pushbutton', 'String', 'Save Data','Units','normalized'...
    ,'Position', [0.8 0.1 0.1 0.05],...
    'Callback', @(btn6,event) Save_intxt (btn6,sld));
   
   
   
% make laser surface for rotation 
    function  Laser_surf(text_Rot,sld,fig_HO,Ax,Ay)
        
        %delete previous laser head
        fucked=findobj(fig_HO,'FaceColor','y');
    delete(fucked);
        
        
            dcm_obj = datacursormode(fig_HO);
c_info = getCursorInfo(dcm_obj);

 Laser_Target=c_info(2).Position;

 pos=c_info(1).Position;
 L_xyz0=pos;
        
  delta=Laser_Target-pos;
 
 v1=[0 0 1]';   % the initial laser head on xy plane
 v2=(delta/norm(delta))';
 Rot_plane=transpose(fcn_RotationFromTwoVectors(v1, v2)); % find the rotation matrix
 
        
        
        
        
        rotation_deg=get(sld,'value');
        set(text_Rot,'String',sprintf('Laser Head Rotation=%d°',rotation_deg));
        
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
 
  surf(X_Actual,Y_Actual,Z_Actual,Power_Actual,'LineStyle','-', 'FaceColor','y'); 
 
% surf(X,Y,Z)
    
    
    
%% make target

function  circle_Laser=Tape_target(btn3,rQ,th_y,fig_HO,circle_Laser)



dcm_obj = datacursormode(fig_HO);
c_info = getCursorInfo(dcm_obj);
Laser_Target=c_info(2).Position;


% fucked=findobj(fig_HO,'color','g','LineStyle','--');
fucked=findobj(fig_HO,'FaceColor','g');
    delete(fucked);
   
    distance=str2double(get(rQ,'String'));

r=distance;

% tt=linspace(0,2*pi,50);
% xx=r*sin(tt);
% yy=r*cos(tt);
% zz=yy*0;

% Circle_Rotation=-(pi/180)*th_y ;
% 
% 
% new_xyz=Transformation_Rot_y(Laser_Target,Circle_Rotation,[xx;yy;zz]);
% 
% 
% circle_Laser=plot3(new_xyz(1,:),new_xyz(2,:),new_xyz(3,:),'g--');



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

% r=distance;
% 
% tt=linspace(0,2*pi,50);
% xx=r*sin(tt);
% yy=r*cos(tt);
% zz=yy*0;
% 
% Circle_Rotation=-(pi/180)*th_y ;
% 
% 
% new_xyz=Transformation_Rot_y(Laser_Target,Circle_Rotation,[xx;yy;zz]);
% 
% 
% circle_Laser=plot3(new_xyz(1,:),new_xyz(2,:),new_xyz(3,:),'g--');

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



function Put ( btn2,fig_HO,text_orientation)


  Fig_handles = guidata(gcbo);


fucked=findobj(fig_HO,'color','r','LineWidth',2);
    delete(fucked);
    
    
       
    
    dcm_obj = datacursormode(fig_HO);
c_info = getCursorInfo(dcm_obj);

 Laser_Target=c_info(2).Position;

 pos=c_info(1).Position;
 
 delta=Laser_Target-pos;
 
 quiver3(pos(1),pos(2),pos(3),delta(1),delta(2),delta(3),0,'k','MaxHeadSize',0.5,'color','r','LineWidth',2);
 
 norm_delta=norm(delta);
 delta=delta/norm_delta;
 set (text_orientation ,'String',sprintf(' nx=%f \n ny=%f \n nz=%f',delta(1),delta(2),delta(3)),'Fontsize',20);
 
 Fig_handles.Laser_direction=delta;
 Fig_handles.Laser_Head_L_xyz0=pos;
  Fig_handles.Laser_Target=Laser_Target;

%   Fig_handles.nip_point_M
  
  
guidata(gcbo,Fig_handles) ;
 
 
 
 
 
%%

% Accept_Nip_point

  
function  Select_Nip_point( btn4,fig_HO)



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

Substrate_Handle=Fig_handles.Substrate_Handle;
% Fig_handles.numberOfErrors = 0; 
% Save the structure
% guidata(fig_HO,Fig_handles) ;

delete(Substrate_Handle);
   
   
    for ii=1:length(Tape_Roller_Handle)
     
     delete(Tape_Roller_Handle{ii});
       
    end
    
%     distance=str2double(get(rQ,'String'));
% 
% r=distance;
% 
% tt=linspace(0,2*pi,50);
% xx=r*sin(tt);
% yy=r*cos(tt);
% zz=yy*0;
% 
% Circle_Rotation=-(pi/180)*th_y ;
% 
% 
% new_xyz=Transformation_Rot_y(Laser_Target,Circle_Rotation,[xx;yy;zz]);
% 
% 
% circle_Laser=plot3(new_xyz(1,:),new_xyz(2,:),new_xyz(3,:),'g--');


dcm_obj = datacursormode(fig_HO);
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on');
datacursormode on;








function  Accept_Nip_point( btn5,fig_HO,...
         Tape_Roller_Handle,...
       N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,R_cyl,...
        z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,...
         rQ,H_indentation)
    
     
     
       Fig_handles = guidata(gcbo);
     
     
   dcm_obj = datacursormode(fig_HO);
c_info = getCursorInfo(dcm_obj);

% center of the nip-point
nip_point_M=c_info.Position;
tv=nip_point_M;

   

     % plot Tape and Roller
[Tape_Roller_Handle]=Tape_plot_3D_ALL_objects_Free_on_Sub(N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,R_cyl,z_cyl_end,H_indentation);
% assignin('base','Tape_Roller_Handle',Tape_Roller_Handle);
    
% plot Substrate

% [points, Boarders,starting_index]=helical_3D_points_4free_onSub(R_cyl(1),z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,nip_point_M);


if R_cyl(1) ~=0

 if nip_point_M (3) >=0 && nip_point_M (3) <=z_cyl_end
 
[points, Boarders,starting_index]=helical_3D_points_4free_onSub(R_cyl(1),z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,nip_point_M);


 elseif nip_point_M (3) < 0  
     % bottom dome part

[points, Boarders,starting_index]=sub_on_ellips(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,nip_point_M);

 elseif nip_point_M (3) > z_cyl_end
     
     % Upper Dome part
[points, Boarders,starting_index]=sub_on_ellips(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,nip_point_M);


 end
 
 else
  
    [points, Boarders,starting_index]=sub_on_flat(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,nip_point_M);

    
end



xnode=2*No_dev -1;
ynode=No_dev_L-starting_index;

X=zeros(xnode,ynode);
Y=zeros(xnode,ynode);
Z=zeros(xnode,ynode);
Intensity=Z;

for kk=1:xnode
    index=(((1:ynode)*xnode)-kk+1);
    X(kk,1:end)=points(1,index);
    Y(kk,1:end)=points(2,index);
    Z(kk,1:end)=points(3,index);
     Intensity(kk,1:end)=0; % RHS in the thermal model
end
% pause(0.5);

% figure
Substrate_Handle=surf(X,Y,Z,Intensity,'LineStyle','None');

Fig_handles.Substrate_Handle=Substrate_Handle;

% Save the structure
Fig_handles.Tape_Roller_Handle=Tape_Roller_Handle;

Fig_handles.nip_point_M=nip_point_M;

guidata(gcbo,Fig_handles) ;




%%

    distance=str2double(get(rQ,'String'));

r=distance;

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





function  Save_intxt( btn6,sld)

%   Fig_handles.Laser_Target=Laser_Target;

%% Data to be written
% Fig_handles.Laser_direction
%  Fig_handles.Laser_Head_L_xyz0
% Fig_handles.nip_point_M=nip_point_M;


fid10 = fopen('.\Geo_Parameter.txt','r');
  out = textscan(fid10,'%s','delimiter',',');
  
  
  fclose (fid10);
  
out_copy=out;
   

% find index to replace
     Rxyz_ind = contains(out{1},'Rxyz');
%    out{1}{Rxyz_ind};
  
    L_xyz0_ind = contains(out{1},'L_xyz0');
%    out{1}{L_xyz0_ind};
   
     tv_ind = contains(out{1},'Roller_Pos_TV');
     
     
      Laser_head_Rot = contains(out{1},'Laser_head_Rot');
     
%       out{1}{tv_ind};
      
      % remove first index to adjust writing new data
      out{1}(1)=[];
      
      
       Fig_handles = guidata(gcbo);
      
       out{1}{Rxyz_ind}= num2str(Fig_handles.Laser_direction);
         out{1}{L_xyz0_ind}= num2str(Fig_handles.Laser_Head_L_xyz0);
              out{1}{tv_ind}=num2str(Fig_handles.nip_point_M);
              
                rotation_deg=get(sld,'value');
%         set(text_Rot,'String',sprintf('Laser Head Rotation=%d°',rotation_deg));
                     
               out{1}{Laser_head_Rot}=num2str(rotation_deg);
              
              fid11 = fopen('.\Geo_Parameter.txt','w');
              % >> regexprep  % to replace in txt file
              
              for ii=1:2:length(out{1})  % number of lines
      
         fprintf(fid11,'%s\n', strcat(out_copy{1}{ii},', ' ,out{1}{ii}) );

    
              end


