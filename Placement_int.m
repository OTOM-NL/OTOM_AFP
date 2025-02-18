
function xyz_int=Placement_int (laser_source_P,laser_direction)

Lx=laser_source_P(1);
Ly=laser_source_P(2);
Lz=laser_source_P(3);

Rx=laser_direction(1);
Ry=laser_direction(2);
Rz= laser_direction(3);

norm_Rxyz=norm([Rx,Ry,Rz]);

Rx=Rx/norm_Rxyz;
Ry=Ry/norm_Rxyz;
Rz= Rz/norm_Rxyz;

% get data in the general axis system

t=-Lz/Rz;
x_int=(Rx*t)+Lx;
y_int=(Ry*t)+Ly;
z_int=0;



if t==Inf
    x_int=[];
z_int=[];
y_int=[];
    
end

   xyz_int=[x_int ;y_int ;z_int ];
            


%   V_AB=xyz_int- laser_source_P;
%             V_AB=V_AB/norm(V_AB);
%             
% %             we assume the normal of surface is [0 0 1]
%             
%             theta_V_AB_Line3D=acosd(V_AB'*laser_direction );
            
%             if abs(theta_V_AB_Line3D)<=1e-1  % parallel
%                 
%                 xyz_int=[x_int ;y_int ;z_int ];
%                 return;
%                 
%             else
%                 xyz_int=[];
%                 
%             end
            
            