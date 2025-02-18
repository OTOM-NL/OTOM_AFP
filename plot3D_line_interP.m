function plot3D_line_interP(R,z_cylinder_end,laser_direction_G,Lx,Ly,Lz,xyz_int,fileID_Int,fileID_Ray)


%plot line and intersection point
%coordinate of 3D line

Rx=laser_direction_G(1);
Ry=laser_direction_G(2);
Rz=laser_direction_G(3);

if ~isempty (xyz_int)
xc=linspace(Lx,xyz_int(1),2);
yc=linspace(Ly,xyz_int(2),2);
zc=linspace(Lz,xyz_int(3),2);
fprintf(fileID_Int,'   %12.8f   %12.8f  %12.8f \r\n',xyz_int(1),xyz_int(2),xyz_int(3));
else
    t=0:1:3;

xc=Rx*t+Lx;
yc=Ry*t+Ly;
zc=Rz*t+Lz;
end

% plot3(xc,yc,zc,'g:','Linewidth',2)


   fprintf(fileID_Ray,' %12.8f    %12.8f   %12.8f  \r\n ',xc(1),yc(1),zc(1) );
  fprintf(fileID_Ray,' %12.8f    %12.8f   %12.8f  \r\n ',xc(2),yc(2),zc(2) );
    fprintf(fileID_Ray,' nan nan nan \r\n');
  




%intersection points
% plot3( x_int, y_int, z_int ,'k.'  )

% text([x_int],[y_int],[z_int], ['P1']);

%   fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3) ,laser_point_Actual_all(ii,4)*First_hit_energy);
%                  fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3) ,laser_point_Actual_all(ii,4)*First_hit_energy);
               

% % % text([x_int],[y_int],[z_int], ['P1';'P2']);
