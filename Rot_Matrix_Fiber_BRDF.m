

function [Rot_Roller_axis,phi_Fiber]=Rot_Matrix_Fiber_BRDF(tv,Fib_or,R_cyl,z_cyl_end,Pr_v,n_xyz)
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


No_dev_L=10;
Hel_Aux=0; %pi/2-alpha; 
% L=R*pi/2;  % length of the tape
% w_T=R/5;  % width of the tape

% It is assumed L_prim == R_cyl, make auxiliary helical path
th_limit=pi/6; %in radian ?!!

 if tv(3) >=0 && tv(3) <=z_cyl_end
 
% gradient in cylender
% n_xyz=[2*tv(1); 2*tv(2) ;0];
% n_xyz=n_xyz/norm(n_xyz);


% starting point
theta0=atan2(x0,y0);
% phi=asin((z0-z_cyl_end)/c); % should be modified, Now Assume on the upper dome


theta=linspace(theta0,theta0+th_limit,No_dev_L);

% initial points
% x0=tv(1);
% y0=tv(2);
zp=tv(3);

H=(a)*sin(Hel_Aux);  % to indicate the length of Tape
x_L=0+(a)*sin(theta);
y_L=0+(a)*cos(theta);
z_L=zp+linspace(0,H,No_dev_L);


 elseif tv(3) < 0  
     % bottom dome part

     % gradient in cylender
%  n_xyz=[-2*tv(1)/a^2; 2*tv(2)/b^2 ;-2*tv(3)/c1^2];
%  n_xyz=n_xyz/norm(n_xyz);


%starting point
x0=tv(1);
y0=tv(2);
z0=tv(3)-0*z_cyl_end;

% c2=1;

% make a zero degree winding as a refference to compare for correction 

wind_angle=Hel_Aux;

wind_angle=wind_angle- (.001);  % to avoid numerical instability

% L=a;
% W=w;

 theta0=atan2(x0,y0);
    phi0=asin(z0/c1);


% wind_angle_lateral=wind_angle;

% [th_limit12,status]=ellipsoid_path_integration(a,b,c1,wind_angle,L,theta0,phi0);


number_of_Div=No_dev_L;
% xyz=zeros(3,2*number_of_Div-1);

% th_limit=th_limit12(2);   % for other side of Tape

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
%  n_xyz=[-2*tv(1)/a^2; 2*tv(2)/b^2 ;-2*(tv(3)-z_cyl_end)/c2^2];
%  n_xyz=n_xyz/norm(n_xyz);

x0=tv(1);
y0=tv(2);
z0=tv(3)-z_cyl_end;

% c2=1;

% make a zero degree winding as a refference to compare for correction 

% For axuilary helical path
wind_angle=Hel_Aux;

wind_angle=wind_angle- (.001);  % to avoid numerical instability

% L=a;


 theta0=atan2(x0,y0);
    phi0=asin(z0/c2);


% [th_limit12,status]=ellipsoid_path_integration(a,b,c2,wind_angle,L,theta0,phi0);

number_of_Div=No_dev_L;
% th_limit= pi/6;
% th_limit12(2);   % for other side of Tape -Axuilary helical path

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


v1=[0 ;0 ;1];   % the initial laser head on xy plane
%  v2=[0 1 0]';  % normal of the surface
 v2=n_xyz;
 Rot_plane=(fcn_RotationFromTwoVectors(v1, v2)); % find the rotation matrix
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

  vector_cross_Pr_v=cross(Pr_v,Vec_tang);

% the value between normal surface and the cross product to see CCW or CW
% should be Rotated 
 theta_normal_cross_product_Pr_v=acosd((vector_cross_Pr_v*n_xyz)/(norm(vector_cross_Pr_v)*norm(n_xyz)));



 phi_Fiber=0-(sign(cosd(theta_normal_cross_product_Pr_v))*acosd((Pr_v*Vec_tang')/(norm(Vec_tang)*norm(Pr_v))));
 



% roller axis vector showing the cylinder path in general coordinate-
% X-axis of local coordinate
xo_points=[0,   a/10];  % just to indicate the direction
yo_points=[0,0];
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


% the value between normal surface and the cross product to see CCW or CW
% should be Rotated 
 theta_normal_cross_product=acosd((vector_cross*n_xyz)/(norm(vector_cross)*norm(n_xyz)));

%  if theta_2Refconfig >90
%      theta_2Refconfig=theta_2Refconfig-180;
%  end

 
%  if  tv(3)<0 || tv(3) > z_cyl_end 
 

%     Fib_or=-sign(cosd(theta_normal_cross_product))* theta_2Refconfig+Fib_or+180;

%  else
%      
    Fib_or=-sign(cosd(theta_normal_cross_product))* theta_2Refconfig+Fib_or+180;

%  end


% Make a Rotation that transform axis to final axis
% Rot_y=[cosd(Fib_or) 0 sind(Fib_or) ; 0 1 0; -sind(Fib_or) 0 cosd(Fib_or)];

Rot_z=([cosd(-Fib_or)  -sind(-Fib_or) 0 ;  sind(-Fib_or)  cosd(-Fib_or) 0 ; 0 0 1]);


Rot_Roller_axis=Rot_plane *Rot_z;

