
% This program is able show 3D objects

function counter_ray=optical_3D_objects(th_y,Tape_Sp,...
    R_cyl,z_cyl_end,tv3,W_R,H_indentation)

counter_ray=1;
% Laser_head,L_xyz0


N_tape=Tape_Sp(1); % EVEN number of points, should not be changed !!
W_tape=Tape_Sp(2) ; % width of the tape
R_tape=Tape_Sp(3);
L_flat=Tape_Sp(4);
thick_T=Tape_Sp(5); % thickness of the tape
deg_tape=Tape_Sp(6);

fig=figure('Name','Laser HO finder','NumberTitle','off');

javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

% fig = gcf; % current figure handle
% fig.ToolBar = 'none';
% fig.MenuBar = 'none';
hold on;


Mandrel_plot=plot3D_cylinder4optical_3D_objects(R_cyl,z_cyl_end);



 theta_ind=acosd((R_tape-H_indentation)/R_tape);  %theta_ind in degree
 
 Nip_Mov =R_tape*sind(theta_ind); 



tv=[0;R_cyl(1)-H_indentation;tv3];  % for the tape
tv_R=[0;R_cyl(1)+R_tape+thick_T-H_indentation;tv3];  % for center of Roller



[nip_point_M]=Tape_plot_3D_ALL_objects(N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,R_cyl(1));


plot3(nip_point_M(1),nip_point_M(2),nip_point_M(3),'r*');

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


distance=.5;
r=distance;

tt=linspace(0,2*pi,50);
xx=r*sin(tt);
yy=r*cos(tt);
zz=yy*0;

Circle_Rotation=-(pi/180)*th_y ;


new_xyz=Transformation_Rot_y(nip_point_M,Circle_Rotation,[xx;yy;zz]);


circle_Laser=plot3(new_xyz(1,:),new_xyz(2,:),new_xyz(3,:),'g--');


dcm_obj = datacursormode(fig);
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on');


set(dcm_obj,'UpdateFcn',@NewCallback_datacursor);

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
    'String','5',...
    'Position',[80 55 120 40],'Visible','on');



btn1 = uicontrol('Style', 'pushbutton', 'String', 'Accept',...
    'Position', [80 10 50 40],...
    'Callback', @(btn1,event) Accept (btn1,rQ,th_y,fig,nip_point_M));

btn2 = uicontrol('Style', 'pushbutton', 'String', 'Put Laser',...
    'Position', [300 10 50 40],...
    'Callback', @(btn2,event) Put (btn2,fig,nip_point_M,text_orientation));



function  Accept( btn1,rQ,th_y,fig,nip_point_M)



fucked=findobj(fig,'color','g','LineStyle','--');

    delete(fucked);

    
    distance=str2double(get(rQ,'String'));

r=distance;

tt=linspace(0,2*pi,50);
xx=r*sin(tt);
yy=r*cos(tt);
zz=yy*0;

Circle_Rotation=-(pi/180)*th_y ;


new_xyz=Transformation_Rot_y(nip_point_M,Circle_Rotation,[xx;yy;zz]);


circle_Laser=plot3(new_xyz(1,:),new_xyz(2,:),new_xyz(3,:),'g--');


dcm_obj = datacursormode(fig);
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on');


set(dcm_obj,'UpdateFcn',@NewCallback_datacursor);

hTarget = handle(circle_Laser);
hDatatip = dcm_obj.createDatatip(hTarget);
set(hDatatip,'Host',hTarget);
datacursormode toggle;



function Put ( btn2,fig,nip_point_M,text_orientation)

fucked=findobj(fig,'color','r','LineStyle','--','LineWidth',4);
    delete(fucked);
    
    

dcm_obj = datacursormode(fig);

 s = getCursorInfo(dcm_obj);
 pos=s.Position;
 
 delta=nip_point_M-pos';
 
 quiver3(pos(1),pos(2),pos(3),delta(1),delta(2),delta(3),0,'k--','MaxHeadSize',0.5,'color','r','LineWidth',4);
 
 norm_delta=norm(delta);
 delta=delta/norm_delta;
 set (text_orientation ,'String',sprintf(' nx=%f \n ny=%f \n nz=%f',delta(1),delta(2),delta(3)),'Fontsize',20);
 
 


  
    






