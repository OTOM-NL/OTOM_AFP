
% Should included some transformation from global axis to local axis
function new_points=Transformation_G2L (tv,Rot_Roller_axis,points)

%tv >> tranformation vector
new_points=[];
if ~isempty(points) 

tv_x =-tv(1);
tv_y =-tv(2); %3*R;
tv_z=-tv(3) ;%10*W;
% TV=-1*[tv_x tv_y tv_z];
%Rotation
% Rot_y=([cosd(th_y) 0 sind(th_y) ; 0 1 0; -sind(th_y) 0 cosd(th_y)]);
% Rot_Roller_axis


x_p=points(1,:);
y_p=points(2,:);
z_p=points(3,:);


%x_G in general coordinate
x_p=x_p+tv_x ; 
y_p=y_p+tv_y;
z_p=z_p+tv_z;


Temp=Rot_Roller_axis\[x_p;y_p;z_p];
x_p=Temp(1,:);
y_p=Temp(2,:);
z_p=Temp(3,:);


new_points=[x_p;y_p;z_p];

end