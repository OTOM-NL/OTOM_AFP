
function [h_T]=Tape_plot_3D_ALL_objects_Free_on_Tape(N,W,R,L_flat,tv,th_y,thick_T,deg_T,W_R,theta_ind,R_cyl)


% W_R=1.5*W;  % width of the Roller

%%
%tv  > translation vector



%%
% use th_y to decresse the y-position of the tape related to the helical 
%%

%plot axes system
xo=0;
yo=0;
zo=0;

tv_x =tv(1);
tv_y =tv(2); %3*R;
tv_z=tv(3) ;%10*W;
% TV=[tv_x tv_y tv_z];
%Rotation
Rot_y=[cosd(th_y) 0 sind(th_y) ; 0 1 0; -sind(th_y) 0 cosd(th_y)];
%x_G in general coordinate

xo_points=xo+[0, R, nan, 0, 0, nan, 0, 0];
yo_points=yo+[0, 0, nan, 0, R, nan, 0, 0];
zo_points=zo+[0, 0, nan, 0, 0, nan, 0, W_R];




Temp=Rot_y*[xo_points;yo_points;zo_points];
xo_points=Temp(1,:);
yo_points=Temp(2,:);
zo_points=Temp(3,:);


%x_G in general coordinate
xo_points=xo_points+tv_x ; 
yo_points=yo_points+tv_y;
zo_points=zo_points+tv_z;


xo_text=[xo+R/2, xo, xo];
yo_text=[yo, yo+R/2, yo];
zo_text=[zo, zo, zo+W_R/2];

Temp2=Rot_y*[xo_text;yo_text;zo_text];
xo_text=Temp2(1,:);
yo_text=Temp2(2,:);
zo_text=Temp2(3,:);

xo_text=xo_text+tv_x;
yo_text=yo_text+tv_y;
zo_text=zo_text+tv_z;


plot3(xo_points, ...
       yo_points, ...
      zo_points , 'k--');
text(xo_text, yo_text, zo_text, ['+X';'+Y';'+Z']);
hold on;

%%



R_thick=R+thick_T;   
%plotting tape S1+S2
% N = 360; % EVEN number of points
% W = 12; % width of the tape
% R=10;
[x,y,z] = cylinder(R_thick,N); % for tape
[x_R,y_R,z_R] = cylinder(R,N);  % for Roller
% x=-x; % change the order of data
y=-y+R_thick;  % for tape
y_R=-y_R+R;    % for roller   I have some doubts....  !!!


deg=ceil(deg_T);

if deg_T< theta_ind
    deg=floor(theta_ind);

end
ceil_theta_ind=floor(theta_ind);

% deg=N/8;  % degree for cylinder section of the tape
x = x(:,N/4+ceil_theta_ind:N/4+deg);  % for continuous ends
% x(:,end) = x(:,1);
y = y(:,N/4+ceil_theta_ind:N/4+deg);
% y(:,end) = y(:,1);
z = z(:,N/4+ceil_theta_ind:N/4+deg);
z(:,end) = z(:,1);    
z(2,:) = z(2,:)*W;
z=z+(W_R-W)/2; % to have a same middle   % it should be in the middle


x_R = x_R(:,[1:(90-ceil_theta_ind+1),N/4+ceil_theta_ind-1:N]);  % for continuous ends
% x(:,end) = x(:,1);
y_R = y_R(:,[1:(90-ceil_theta_ind+1),N/4+ceil_theta_ind-1:N]);
% y(:,end) = y(:,1);
z_R = z_R(:,[1:(90-ceil_theta_ind+1),N/4+ceil_theta_ind-1:N]);
% N/4+ceil_theta_ind


%%
phi=-th_y*(pi/180); %pi/2-alpha; 

xp=0;
yp=0;
zp=tv_z;
L=R*theta_ind*(pi/180);

th=L*cos(phi)/R_cyl; %in radian
H=L*sin(phi);

No_dev_L=10;


t_L=linspace(-th,th,No_dev_L);
x_L=xp+R_cyl*sin(t_L);
y_L=yp+R_cyl*cos(t_L);
z_L=zp+linspace(-H,H,No_dev_L);

plot3(x_L,y_L,z_L,'k--');


% th_up=0.5*W_R*cos(phi+pi/2)/0.3; %in radian
% H=0.5*W_R*sin(phi+pi/2);
% 
% No_dev=2;
% 
% t_up=linspace(-th_up,th_up,No_dev);
% x_up=xp+0.3*sin(t_up);
% y_up=yp+0.3*cos(t_up);
% z_up=zp+linspace(-H,H,No_dev);
% 
% plot3(x_up,y_up,z_up,'r.-');
% 
% points=[x_up;y_up; z_up];
% No=(No_dev_L)*No_dev;

% for pp=1:No_dev_L
%     
%     tv=[0, 0 z_L(pp)-tv_z];
% 
%     Rot=t_L(pp);
% 
%     new_points=Transformation_Rot_z (tv,Rot,points);
% %       new_points2=Transformation_Rot_z (tv,Rot,points2);
% 
%   plot3(new_points(1,:),new_points(2,:),new_points(3,:),'r.-');
%       
% % Roller_Def_points(1:3,(1:No_dev) + (pp-1)*No_dev )=new_points(1:3,:);
% 
% Roller_Def_points_X(1:2,pp)= new_points(1,1:2)';
% Roller_Def_points_Y(1:2,pp)= new_points(2,1:2)';
% Roller_Def_points_Z(1:2,pp)= new_points(3,1:2)';
%     
% end


% surf(Roller_Def_points_X,Roller_Def_points_Y,Roller_Def_points_Z)

%%





z_R(2,:) = z_R(2,:)*W_R;

%% to be in the middle
%-W_R/2;   % add -W_R/2
z_R=z_R(1:2,:) -W_R/2;
z=z(1:2,:) -W_R/2;
%%

% x=x+R;

%%
% Rotation matrix and translation seems can be done here
%Y-Rotation
% th_y=-90;


Temp=zeros(3,1);
Temp2=zeros(3,1);
Temp_R=zeros(3,1);
Temp2_R=zeros(3,1);

x_G=x;
y_G=y;
z_G=z;

x_G_R=x_R;
y_G_R=y_R;
z_G_R=z_R;
for i=1:size(x,2)
%     Temp=[x(1,i),y(1,i),z(1,i)];
Temp=Rot_y*[x(1,i);y(1,i);z(1,i)];
Temp2=Rot_y*[x(2,i);y(2,i);z(2,i)];

x_G(1,i)=Temp(1);
y_G(1,i)=Temp(2);
z_G(1,i)=Temp(3);

x_G(2,i)=Temp2(1);
y_G(2,i)=Temp2(2);
z_G(2,i)=Temp2(3);

end

for i=1:size(x_R,2)
   Temp_R=Rot_y*[x_R(1,i);y_R(1,i);z_R(1,i)];
Temp2_R=Rot_y*[x_R(2,i);y_R(2,i);z_R(2,i)];


x_G_R(1,i)=Temp_R(1);
y_G_R(1,i)=Temp_R(2);
z_G_R(1,i)=Temp_R(3);

x_G_R(2,i)=Temp2_R(1);
y_G_R(2,i)=Temp2_R(2);
z_G_R(2,i)=Temp2_R(3); 
end


%%
%trnslation vector
% tv_x =0;
% tv_y =3*R;
% tv_z=10*W;
% TV=[tv_x tv_y tv_z];
%x_G in general coordinate
x_G=x_G+tv_x ; 
y_G=y_G+tv_y;
z_G=z_G+tv_z;   % Tape should be in the middle of the Roller

x_G_R=x_G_R+tv_x ; 
y_G_R=y_G_R+tv_y;
z_G_R=z_G_R+tv_z;
%%
% ax1 = axes;


xnode_Temp=length(z_G_R);
ynode_Temp=(2);

X=zeros(xnode_Temp,ynode_Temp);
Y=zeros(xnode_Temp,ynode_Temp);
Z=zeros(xnode_Temp,ynode_Temp);
X2=zeros(xnode_Temp,ynode_Temp);
Y2=zeros(xnode_Temp,ynode_Temp);
Z2=zeros(xnode_Temp,ynode_Temp);

Temp=zeros(xnode_Temp,ynode_Temp);

for kk=1:xnode_Temp
 
    X(kk,1:end)=x_G_R(1,kk)';
    Y(kk,1:end)=y_G_R(1,kk)';
    Z(kk,1:end)=z_G_R(1,kk)';
    
    X2(kk,1:end)=x_G_R(2,kk)';
    Y2(kk,1:end)=y_G_R(2,kk)';
    Z2(kk,1:end)=z_G_R(2,kk)';
   
end

% X(end,:)=X(1,:);
% Y(end,:)=Y(1,:);
% Z(end,:)=Z(1,:);

% figure(30);

h_R1=fill3(X,Y,Z,Temp,'FaceAlpha',0.8);
h_R2=fill3(X2,Y2,Z2,Temp,'FaceAlpha',0.8);


h_R0=surf(x_G_R,y_G_R,z_G_R,'FaceColor',[0.95 0.7 1]);

h_R=[h_R0;h_R1;h_R2];
% camlight left
camlight(-40,150)
axis equal
hold on
% shading flat
% 'FaceColor','interp','FaceLighting','gouraud'
% shading flat

% view(2)
% ax2 = axes;
% if ~isvector(z_G)
% h2=surf(x_G,y_G,z_G,'FaceColor',[0.278 0.788 0.788]);
% end

nip_point_M=zeros(3,1);
nip_point_M(1)=0.5*(x_G(1,1)+x_G(2,1));
nip_point_M(2)=0.5*(y_G(1,1)+y_G(2,1));
nip_point_M(3)=0.5*(z_G(1,1)+z_G(2,1));

% colormap winter
% shading interp


% linkaxes([ax1,ax2])
% ax2.Visible = 'off';
% ax2.XTick = [];
% ax2.YTick = [];


% colormap(ax1,'hot')
% colormap(ax2,'cool')

% ylim([-1 1])
%%  x(:,end) and y(:,end) should be chosen as start point of surface

 % the long of flat part of the tape
[Z,Y]=meshgrid(linspace(z(1,1),z(2,1),2),linspace(y(1,end)+(L_flat*sind(deg_T)),y(1,end),2));
% [Z,Y]=meshgrid(linspace(W_R/4,W+W_R/4,2),linspace(y(1,end)+(L_flat*cosd(deg)),y(1,end),2));
% Z=2*X-8*Y-18
% surf(X,Y,Z)
n1=- R_thick*sind(deg_T);%   x(1,end); %input('direction of ray in x-direction =');
n2=- R_thick*cosd(deg_T);% y(1,end)-R_thick;% input('direction of ray in y-direction =');
n3= 0.0; %input('direction of ray in z-direction =');
% normalize the Rx,Ry, Rz
normR=norm([n1 n2 n3]);
n1=n1/normR;
n2=n2/normR;
n3=n3/normR;
n=[n1 n2 n3];
% a point on the surface of the plane
x0=x(1,end);
y0=y(1,end);
z0=0;
% plane equation  : n1*(X-x0)+n2(Y-y0)+n3(Z-z0)=0
d=n1*x0+n2*y0+n3*z0;


% x_end_flat=L_flat*cosd(deg_T)+ x(1,end)
X=(n3*Z+n2*Y-d)/(-n1);

%%
% it can be done here for transformation of points
% th_y=0;
% Rot_y=[cosd(th_y) 0 sind(th_y) ; 0 1 0; -sind(th_y) 0 cosd(th_y)];
X_G=X;
Y_G=Y;
Z_G=Z;
for i=1:size(X,2)
%     Temp=[x(1,i),y(1,i),z(1,i)];
Temp=Rot_y*[X(1,i);Y(1,i);Z(1,i)];
Temp2=Rot_y*[X(2,i);Y(2,i);Z(2,i)];
X_G(1,i)=Temp(1);
Y_G(1,i)=Temp(2);
Z_G(1,i)=Temp(3);

X_G(2,i)=Temp2(1);
Y_G(2,i)=Temp2(2);
Z_G(2,i)=Temp2(3);
end

X_G=X_G+tv_x ;
Y_G=Y_G+tv_y;
Z_G=Z_G+tv_z; % Tape should be in the middle of the Roller
%%

h3=surf([x_G X_G],[y_G Y_G],[z_G Z_G],'FaceColor',[0.278 0.788 0.788]);

    h_T= h3; %for the flat and curved part of the tape

% shading flat
% axis([-4*R 3*R -4*R 3*R]) % it should be modified for rotation and translation matrix

