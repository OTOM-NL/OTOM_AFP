
function [kr1,beta]=refl_Roller_plot3D(R_Roller,xyz_int,laser_direction,tv,Rot_Roller_axis,fileID_Normals,fileID_Refs,W_R,varargin)


 x_int=xyz_int(1);
 y_int=xyz_int(2);
 z_int=xyz_int(3);
 
  Rx=laser_direction(1);
  Ry=laser_direction(2);
  Rz= laser_direction(3);
 
 

t=linspace(0,R_Roller/5,2);


S_normal=zeros(1,3);
Ref_line=zeros(1,3);


Tol=1e-3;


if x_int
if (abs(x_int^2+y_int^2 - R_Roller^2)/(R_Roller^2) < Tol)  % check if it is on the surface

    norm_n=norm([x_int y_int]);
    
   %normal in x-y plane 
S_normal(1,1:2)=([x_int(1) y_int ])/norm_n;

if isempty(varargin)
    
normal_line_X(1,:)=S_normal(1,1)*t+x_int;
normal_line_Y (1,:)=S_normal(1,2)*t+y_int;

normal_G=Transformation_L2G (tv,Rot_Roller_axis,[normal_line_X(1,:);normal_line_Y(1,:);[z_int z_int]]);
% normal_G(:,2)=normal_G(:,1)*5;
% plot3(normal_G(1,:),normal_G(2,:),normal_G(3,:),'y:');


   fprintf(fileID_Normals,' %12.8f    %12.8f   %12.8f  \r\n ',normal_G(1,1),normal_G(2,1),normal_G(3,1));
  fprintf(fileID_Normals,' %12.8f    %12.8f   %12.8f  \r\n ',normal_G(1,2),normal_G(2,2),normal_G(3,2));
    fprintf(fileID_Normals,' nan nan nan \r\n');

end
% text(normal_line_X(1,2),normal_line_Y(1,2),z_int , sprintf('n'));

% S_normal(1,1:3)=[S_normal(1,1:2) 0 ];
else
   
 S_normal(1,1:3)=[0 0 sign(z_int-W_R/2)]  ;  % to determine which side of the cylinder is hited by ray
 
 if isempty(varargin)
     
 normal_line_Z (1,:)=S_normal(1,3)*t+z_int;
 normal_G=Transformation_L2G (tv,Rot_Roller_axis,[[x_int x_int];[y_int y_int];normal_line_Z]);
%  plot3(normal_G(1,:),normal_G(2,:),normal_G(3,:),'y:');
 
  fprintf(fileID_Normals,' %12.8f    %12.8f   %12.8f  \r\n ',normal_G(1,1),normal_G(2,1),normal_G(3,1));
  fprintf(fileID_Normals,' %12.8f    %12.8f   %12.8f  \r\n ',normal_G(1,2),normal_G(2,2),normal_G(3,2));
    fprintf(fileID_Normals,' nan nan nan \r\n');
 end
%     hold on
end








Ray=([Rx Ry Rz]);

kr1=zeros(1,3);


% if S_normal(1,1)

kr1(1,1:3)=2*(-Ray*S_normal(1,:)')*S_normal(1,:)-(-Ray);
%kr1 direction should be outward of the surface !

if isempty(varargin)
    
Ref_line_X(1,:)=kr1(1,1)*t+x_int(1);
Ref_line_Y (1,:)=kr1(1,2)*t+y_int(1);
Ref_line_Z (1,:)=kr1(1,3)*t+z_int(1);
% plot3(Ref_line_X(1,:),Ref_line_Y(1,:),Ref_line_Z(1,:),'r--')


Ref_G=Transformation_L2G (tv,Rot_Roller_axis,[Ref_line_X(1,:);Ref_line_Y(1,:);Ref_line_Z(1,:)]);
% Ref_G(:,2)=Ref_G(:,1)*5;
% plot3(Ref_G(1,:),Ref_G(2,:),Ref_G(3,:),'b--');

 fprintf(fileID_Refs,' %12.8f    %12.8f   %12.8f  \r\n ',Ref_G(1,1),Ref_G(2,1),Ref_G(3,1));
  fprintf(fileID_Refs,' %12.8f    %12.8f   %12.8f  \r\n ',Ref_G(1,2),Ref_G(2,2),Ref_G(3,2));
    fprintf(fileID_Refs,' nan nan nan \r\n');
end
% text(Ref_line_X(1,2),Ref_line_Y(1,2),Ref_line_Z(1,2) , sprintf('Ref'));

% Ref_line=[Ref_line_X;Ref_line_Y;Ref_line_Z];  % reflection line from the surface

end

% if S_normal(2,1)
% kr2 (1,1:3)=2*(Ray*S_normal(2,:)') *S_normal(2,:)-Ray;
% 
% 
% Ref_line_X(2,:)=kr2(1,1)*(-t)+x_int(2);
% Ref_line_Y (2,:)=kr2(1,2)*(-t)+y_int(2);
% Ref_line_Z (2,:)=kr2(1,3)*(-t)+z_int(2);
% plot3(Ref_line_X(2,:),Ref_line_Y(2,:),Ref_line_Z(2,:),'b--')
% text(Ref_line_X(2,2),Ref_line_Y(2,2),Ref_line_Z(2,2) , sprintf('Ref%d',2));
% 
% 
% 
% 
% end



kr1=kr1';  % just for output in general program

beta=0.5*acosd(dot(laser_direction,kr1)/(norm(laser_direction)*norm(kr1)));
% found in 16 Jan 2020
% \beta is the angle between reflection and laser for measuring
% %                 angle of incident
% 0.5* is because the angle between incident ray and reflection is twice as
% angle between normal and incident ray

kr1=Transformation_L2G (zeros(3,1),Rot_Roller_axis,kr1);
