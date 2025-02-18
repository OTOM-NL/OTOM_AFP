
% Laser Head distribution management




function Divergence_factor=laser_Divergence_management (Ax,Ay,nx,ny,Max_Div_Angle)


% Vector which makes maximum divergence >> at the edges
Vec_C_point_actual=[Ax,Ay,0]; %-0

% normal of the surface
v2=[0;0;1];
% plot3(laser_point_Actual(1),laser_point_Actual(2),laser_point_Actual(3),'r*');

% b2=v2'+ (Vec_C_point_actual*Divergence_factor);

% Find the Divergence factor
F_Div=@ (Divergence_factor) acosd((((Divergence_factor*Vec_C_point_actual)+v2')*v2)/norm((Divergence_factor*Vec_C_point_actual)+v2'))-Max_Div_Angle;

x0 = [0 500]; % initial interval
Divergence_factor = fzero(F_Div,x0);



area=2*Ax*2*Ay;


% Should be changed by the user by scrolling



% v1=[0 0 1]';   % the initial laser head on xy plane
%  v2=[Rx Ry Rz]';
%  Rot_plane=transpose(fcn_RotationFromTwoVectors(v1, v2)); % find the rotation matrix


xp_laser=linspace(-Ax,Ax,nx);
yp_laser=linspace(-Ay,Ay,ny);
[X,Y]=meshgrid(xp_laser,yp_laser);  % ny * nx = size
Z=zeros(size(X));

% laser_point_Actual=zeros(1,3);
% laser_point_Actual_all=zeros(nx*ny,4);




% Vec_C_point_actual=laser_point_Actual; %-0
% plot3(laser_point_Actual(1),laser_point_Actual(2),laser_point_Actual(3),'r*');



length_line=max([Ax,Ay]);


No=nx*ny;
Points=zeros(No*3,3);

counter=1;
for ii=1:nx
    for jj=1:ny
    
%         counter=counter+1;
%         laser_point_Actual_all(counter,4)=amp*fun(xp_laser(ii),yp_laser(jj));
        
%         Power_Actual(jj,ii)= laser_point_Actual_all(counter,4);
%         index=(ii-1)*ny+jj
        
        Vec_C_point_actual=[xp_laser(ii),yp_laser(jj),0];
        b2=v2'+ (Vec_C_point_actual*Divergence_factor);
        % to be in harmony with size of the laser head
b2=(b2/norm(b2))*length_line;
% make the line discrete and show only one-time for graphical performance
Points(counter:counter+2,:)=[Vec_C_point_actual;Vec_C_point_actual+b2;nan nan nan];
       
counter=counter+3;
    end
end


% figure;

%
surf(X,Y,Z,Z,'LineStyle','None');
hold on;
plot3(Points(:,1),Points(:,2),Points(:,3),'c--');
axis equal;
 colorbar off;
view([20 10]);






end


