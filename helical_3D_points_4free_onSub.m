% clc;
% clear;

function [thermal_points,Ps,starting_index]=helical_3D_points_4free_onSub(R,z_end,L,w_T,sui,No_dev_L,No_dev,Nip_Mov,nip_point_M)


%% Helical path on cylinder 


% R=10;  % radius of cylinder
% z_end=50;  % long of cylinder

% number_of_section=50;
% [X,Y,Z] = cylinder(R,number_of_section);
% Z=Z*z_end;
% Mandrel_plot=mesh(X,Y,Z);

% h1=figure(1);
%   set(h1,'Visible','off');
% hold on;
% axis equal;
%%


R=R+R/100; % just for representation

% alpha=pi/4; % winding angle
phi=sui*(pi/180); %pi/2-alpha; 
% L=R*pi/2;  % length of the tape
% w_T=R/5;  % width of the tape

% initial points
x0=nip_point_M(1);
y0=nip_point_M(2);
zp=nip_point_M(3);

th=L*cos(phi)/R; %in radian
H=L*sin(phi);

 theta0=atan2(x0,y0);

% No_dev_L=40;
t_L=linspace(theta0,theta0+th,No_dev_L);
x_L=0+R*sin(t_L);
y_L=0+R*cos(t_L);
z_L=zp+linspace(0,H,No_dev_L);

% plot3(x_L,y_L,z_L,'w.-');

% No_dev=10;


%%

% xp=x_L(kk);
% yp=y_L(kk)-R;
% zp=z_L(kk);

th_up=w_T*cos(phi+pi/2)/R; %in radian
H=w_T*sin(phi+pi/2);


t_up=linspace(theta0,theta0+th_up,No_dev);
x_up=0+R*sin(t_up);
y_up=0+R*cos(t_up);
z_up=zp+linspace(0,H,No_dev);

% plot3(x_up,y_up,z_up,'r.-');
%%
% xp=x_L(kk);
% yp=y_L(kk)-R;
% zp=z_L(kk);

th_low=w_T*cos(phi-pi/2)/R; %in radian
H=w_T*sin(phi-pi/2);

t_low=linspace(theta0,theta0+th_low,No_dev);
x_low=0+R*sin(t_low);
y_low=0+R*cos(t_low);
z_low=zp+linspace(0,H,No_dev);

% plot3(x_low,y_low,z_low,'g.-');

%Translation to the new position, Rotation around z

points=[x_up;y_up; z_up];
points2=[x_low;y_low; z_low];

No_row_nodes=(2*No_dev)-1;


thermal_points_Temp=zeros(3,No_dev_L*No_dev);



starting_index=floor(No_dev_L*(Nip_Mov/L)); % after this position , the nip-point will be defined
% starting_index=1;
No=(No_dev_L-starting_index)*No_row_nodes;

t_L2=linspace(0,th,No_dev_L);


for pp=1:No_dev_L
    
    tv=[0, 0, z_L(pp)-zp];
% tv=zeros(3,1);
    Rot=t_L2(pp);

    new_points=Transformation_Rot_z (tv,Rot,points);
      new_points2=Transformation_Rot_z (tv,Rot,points2);

      
thermal_points_Temp(1:3,(1:No_row_nodes) + (pp-1)*No_row_nodes )=[new_points2(1:3,end:-1:1),new_points(1:3,2:end)];

%     plot3(new_points(1,:),new_points(2,:),new_points(3,:),'r.-');
%      plot3(new_points2(1,:),new_points2(2,:),new_points2(3,:),'g.-');
    
end

thermal_points=thermal_points_Temp(1:3,1+((starting_index)*No_row_nodes):end);


% Collect boarder points

P1s=thermal_points(1:2:3,(1:No_row_nodes-1));
P2s=thermal_points(1:2:3,No_row_nodes:No_row_nodes:No);
P3s=thermal_points(1:2:3,No-1:-1:No-No_row_nodes+2);
P4s=thermal_points(1:2:3,No-No_row_nodes+1:-No_row_nodes:No_row_nodes+1);

Ps=[P1s, P2s,P3s, P4s];


%  theta_all=on_out_surface (p,Ps)



