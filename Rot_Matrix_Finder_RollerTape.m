
function [Rot_Roller_axis,n_xyz]=Rot_Matrix_Finder_RollerTape(tv,th_y,R_cyl,z_cyl_end)


%%

% 26 April 2019
if tv(1)==0
    tv(1)=tv(1)+0.000001;
end

if tv(2)==0
    tv(2)=tv(2)+0.000001;
end
if tv(3)==0
    tv(3)=tv(3)+0.000001;
end





%tv  > translation vector

 if length(R_cyl)==3
   
    c1=R_cyl(2);
     c2=R_cyl(3);
else
       c1=R_cyl(1);
     c2=R_cyl(1);
%       R_cyl=R_cyl(1);
    
end
a=R_cyl(1);
b=R_cyl(1);
% c1=3;



x0=tv(1);
y0=tv(2);
z0=tv(3);


% Arbitrary number to discritize the path
No_dev_L=50;
th_y2=th_y*(pi/180); %pi/2-alpha; 
% L=R*pi/2;  % length of the tape
% w_T=R/5;  % width of the tape

% It is assumed L_prim == R_cyl
th_limit=cos(th_y2); %in radian

 if tv(3) >=0 && tv(3) <=z_cyl_end
 
% gradient in cylender
n_xyz=[-2*tv(1); 2*tv(2) ;0];
n_xyz=n_xyz/norm(n_xyz);



theta0=atan2(x0,y0);
% phi=asin((z0-z_cyl_end)/c); % should be modified, Now Assume on the upper dome


theta=linspace(theta0,theta0+th_limit,No_dev_L);

% initial points
% x0=tv(1);
% y0=tv(2);
zp=tv(3);

H=(a)*sin(th_y2);  % to indicate the length of Tape
x_L=0+(a)*sin(theta);
y_L=0+(a)*cos(theta);
z_L=zp+linspace(0,H,No_dev_L);


 elseif tv(3) < 0  
     % bottom dome part

     % gradient in cylender
n_xyz=[-2*tv(1)/a^2; 2*tv(2)/b^2 ;-2*tv(3)/c1^2];
n_xyz=n_xyz/norm(n_xyz);

x0=tv(1);
y0=tv(2);
z0=tv(3)-0*z_cyl_end;

% c2=1;

% make a zero degree winding as a refference to compare for correction 

wind_angle=0*th_y2;

wind_angle=wind_angle- (.001);  % to avoid numerical instability

L=a;
% W=w;

 theta0=atan2(x0,y0);
    phi0=asin(z0/c1);


% wind_angle_lateral=wind_angle;

[th_limit12,status]=ellipsoid_path_integration(a,b,c1,wind_angle,L,theta0,phi0);


number_of_Div=No_dev_L;
% xyz=zeros(3,2*number_of_Div-1);

th_limit=th_limit12(2);   % for other side of Tape

theta=linspace(theta0,theta0+th_limit,number_of_Div);



theta_range=linspace(0,th_limit,number_of_Div);
% theta=linspace(pi/2,1*pi,160);
phi=asin(z0/c1)+ tan(wind_angle)*theta_range;
% phi=pi/6;


% phi0=asin(z0/c1);

% first trajectory
% xp=a.*cos(0).*sin(theta);
% yp=b.*cos(0).*cos(theta);
% zp=c1.*sin(phi); %.*ones(size(theta));


% second trajectory
x_L=a.*cos(phi).*sin(theta);
y_L=b.*cos(phi).*cos(theta);
z_L=c1.*sin(phi); %.*ones(size(theta));


% plot3(x_L,y_L,0+z_L,'w.-');
% hold on;

     
 elseif tv(3) > z_cyl_end
     
     % For upper dome
     
     % gradient in cylender
n_xyz=[-2*tv(1)/a^2; 2*tv(2)/b^2 ;-2*(tv(3)-z_cyl_end)/c2^2];
n_xyz=n_xyz/norm(n_xyz);

x0=tv(1);
y0=tv(2);
z0=tv(3)-z_cyl_end;

% c2=1;

% make a zero degree winding as a refference to compare for correction 

wind_angle=0*th_y2;

wind_angle=wind_angle- (.001);  % to avoid numerical instability

L=a;


 theta0=atan2(x0,y0);
    phi0=asin(z0/c2);


[th_limit12,status]=ellipsoid_path_integration(a,b,c2,wind_angle,L,theta0,phi0);

number_of_Div=No_dev_L;
th_limit=th_limit12(2);   % for other side of Tape

theta=linspace(theta0,theta0+th_limit,number_of_Div);

theta_range=linspace(0,th_limit,number_of_Div);
% theta=linspace(pi/2,1*pi,160);
phi=asin(z0/c2)+ tan(wind_angle)*theta_range;


% second trajectory
x_L=a.*cos(phi).*sin(theta);
y_L=b.*cos(phi).*cos(theta);
z_L=c2.*sin(phi); %.*ones(size(theta));


% plot3(x_L,y_L,z_cyl_end+z_L,'w.-');
% hold on;

  

 end


v1=[0 ;1 ;0];   % the initial laser head on xy plane
%  v2=[0 1 0]';  % normal of the surface
 v2=n_xyz;
 Rot_plane=transpose(fcn_RotationFromTwoVectors(v1, v2)); % find the rotation matrix
%  



%%
%% This part fix the roller to the desired direction
% CAN be modified
% a=R_cyl(1);
% b=R_cyl(1);
% c=R_cyl(1);

% WE do correction because the transformed axis is rotated around initial
% y-axis, This correction make sure that the rotated 3D-axis conform the
% final axis system

% before the correction 
% plot3(x_L,y_L,z_L,'b.-');

% vector of the tangent 
Vec_tang=[x_L(2)-x_L(1),y_L(2)-y_L(1),z_L(2)-z_L(1)];

% roller axis vector showing the cylinder path in general coordinate-
% X-axis of local coordinate
xo_points=[0, a/100];  % just to indicate the direction
yo_points=[0, 0];
zo_points=[0, 0];


% execute transformation on Local coordinate axis to check difference
% between transformed vector and helical path

Temp=Rot_plane*[xo_points;yo_points;zo_points];
xo_points=Temp(1,:)+tv(1);
yo_points=Temp(2,:)+tv(2);
zo_points=Temp(3,:)+tv(3);

Vec_Xaxis_roller=[xo_points(2)-xo_points(1),yo_points(2)-yo_points(1),zo_points(2)-zo_points(1)];

% plot3(xo_points,yo_points,zo_points,'r-')


 theta_2Refconfig=acosd((Vec_Xaxis_roller*Vec_tang')/(norm(Vec_tang)*norm(Vec_Xaxis_roller)));

 vector_cross=cross(Vec_Xaxis_roller,Vec_tang);
%  Vec_tang is the refference curve 

%% 16April 2019
% because we assume a stupid normal of surface vector
n_xyz([1 3])=-n_xyz([1 3]);
%%


% the value between normal surface and the cross product to see CCW or CW
% should be Rotated 
 theta_normal_cross_product=acosd((vector_cross*n_xyz)/(norm(vector_cross)*norm(n_xyz)));


 
%  if  tv(3)<0 || tv(3) > z_cyl_end 
%  
%      % in upper or bottom Dome parts
%      
% th_y=sign(theta_normal_cross_product-90)* theta_2Refconfig+th_y;
% %  elseif  
% %      % in cylinderical parts
% %      th_y=-theta_2Refconfig+th_y;
%  else
%      
%     th_y=-theta_2Refconfig; 
% 
%  end


% 16th April 2019
Dir=180;

% In Uot model, we make the roller direction zero to the path direction.
% The path direction in UOT is not Zero like here!!


 if tv(3) >=0 && tv(3) <=z_cyl_end
     th_y=-sign(cosd(theta_normal_cross_product))* theta_2Refconfig + Dir;%+ th_y;
 else
 
th_y=sign(cosd(theta_normal_cross_product))* theta_2Refconfig + Dir+ th_y;
 end

% Make a Rotation that transform axis to final axis
Rot_y=[cosd(th_y) 0 sind(th_y) ; 0 1 0; -sind(th_y) 0 cosd(th_y)];

% The rotation matrix from initial Config to Final Config of coordinate
% axis
% Rot_z=([cos(-theta0)  -sin(-theta0) 0 ;  sin(-theta0)  cos(-theta0) 0 ; 0 0 1]);


Rot_Roller_axis=Rot_plane*Rot_y;

