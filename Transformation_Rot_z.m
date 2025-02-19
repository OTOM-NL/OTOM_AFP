
% Should included some transformation from global axis to local axis
function new_points=Transformation_Rot_z(tv,th_y,points)

%tv >> tranformation vector
new_points=[];
if ~isempty(points) 

tv_x =tv(1);
tv_y =tv(2); %3*R;
tv_z=tv(3) ;%10*W;
% TV=-1*[tv_x tv_y tv_z];
%Rotation
Rot_z=([cos(th_y)  -sin(th_y) 0 ;  sin(th_y)  cos(th_y) 0 ; 0 0 1]);


x_p=points(1,:);
y_p=points(2,:);
z_p=points(3,:);




Temp=Rot_z\[x_p;y_p;z_p];
x_p=Temp(1,:);
y_p=Temp(2,:);
z_p=Temp(3,:);


%x_G in general coordinate
x_p=x_p+tv_x ; 
y_p=y_p+tv_y;
z_p=z_p+tv_z;


new_points=[x_p;y_p;z_p];

end