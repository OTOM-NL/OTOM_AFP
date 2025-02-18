
function [kr1,beta]=refl_Mandrel_plot3D(R_cyl,xyz_int,laser_direction,z_cyl_end,fileID_Normals,fileID_Refs,varargin)




if length(R_cyl)==3
   
    c1=R_cyl(2);
     c2=R_cyl(3);
     
else
       c1=R_cyl(1);
     c2=R_cyl(1);
%       R_cyl=R_cyl(1);
    
end



 x_int=xyz_int(1);
 y_int=xyz_int(2);
 z_int=xyz_int(3);
 
  Rx=laser_direction(1);
  Ry=laser_direction(2);
  Rz= laser_direction(3);
 
 

% t=linspace(0,R_cyl(1)/10,2);


S_normal=zeros(1,3);
% Ref_line=zeros(1,3);



if ~isempty(x_int)



if R_cyl(1) ~=0

    t=linspace(0,R_cyl(1)/10,2);




    
    if z_int >= 0 & z_int <=z_cyl_end
if (abs(x_int^2+y_int^2 - R_cyl(1)^2) < 1e-3)   % Seems not necessary

    norm_n=norm([x_int y_int]);
    
   %normal in x-y plane 
S_normal(1,1:2)=([x_int(1) y_int ])/norm_n;


if isempty(varargin)
normal_line_X(1,:)=S_normal(1,1)*t+x_int;
normal_line_Y (1,:)=S_normal(1,2)*t+y_int;

% plot3(normal_line_X(1,:),normal_line_Y(1,:),[z_int z_int],'y:');


   fprintf(fileID_Normals,' %12.8f    %12.8f   %12.8f  \r\n ',normal_line_X(1,1),normal_line_Y(1,1),z_int);
  fprintf(fileID_Normals,' %12.8f    %12.8f   %12.8f  \r\n ',normal_line_X(1,2),normal_line_Y(1,2),z_int);
    fprintf(fileID_Normals,' nan nan nan \r\n');
end

% text(normal_line_X(1,2),normal_line_Y(1,2),z_int , sprintf('n'));
% hold on
% else
%  S_normal(1,1:3)=[0 0 sign(z_int-R_cyl(1)/2)]  ;  % to determine which side of the cylinder is hited by ray
%  normal_line_Z (1,:)=S_normal(1,3)*t+z_int;
% 
%  plot3([x_int x_int],[y_int y_int],normal_line_Z,'b:');
%     hold on   

end


    elseif z_int >= z_cyl_end
        %for upper ellipsoid  2nd Dome
          %gradiant calculations
          
           a=R_cyl(1);
            b=R_cyl(1);
%             c1=R_cyl;
            
          
            x0_E1=0;
            y0_E1=0;
            z0_E1=z_cyl_end;
            
            
                    nx1=(x_int-x0_E1)/a^2;
                    ny1=(y_int-y0_E1)/b^2;
                    nz1=(z_int-z0_E1)/c2^2;
                    
                    Norm_nxyz= norm([nx1,ny1,nz1]);
                    nx1=nx1/Norm_nxyz;
                    ny1=ny1/Norm_nxyz;
                    nz1=nz1/Norm_nxyz;
                    n_xyz1=[nx1 ny1 nz1];
                    
                    S_normal(1,1:3)=n_xyz1;
if isempty(varargin)
normal_line_X(1,:)=S_normal(1,1)*t+x_int;
normal_line_Y (1,:)=S_normal(1,2)*t+y_int;
normal_line_Z (1,:)=S_normal(1,3)*t+z_int;

% plot3(normal_line_X(1,:),normal_line_Y(1,:),normal_line_Z(1,:),'y:');


   fprintf(fileID_Normals,' %12.8f    %12.8f   %12.8f  \r\n ',normal_line_X(1,1),normal_line_Y(1,1),normal_line_Z(1,1));
  fprintf(fileID_Normals,' %12.8f    %12.8f   %12.8f  \r\n ',normal_line_X(1,2),normal_line_Y(1,2),normal_line_Z(1,2));
    fprintf(fileID_Normals,' nan nan nan \r\n');
        
end

   elseif z_int <= 0
       % for lower ellipsoid  1st Dome
       
    %gradiant calculations
            a=R_cyl(1);
            b=R_cyl(1);
%             c=R_cyl(1);
            
          
            x0_E2=0;
            y0_E2=0;
            z0_E2=0;
            
            
                    nx1=(x_int-x0_E2)/a^2;
                    ny1=(y_int-y0_E2)/b^2;
                    nz1=(z_int-z0_E2)/c1^2;
                    
                    Norm_nxyz= norm([nx1,ny1,nz1]);
                    nx1=nx1/Norm_nxyz;
                    ny1=ny1/Norm_nxyz;
                    nz1=nz1/Norm_nxyz;
                    n_xyz1=[nx1 ny1 nz1];
                    
                    S_normal(1,1:3)=n_xyz1;

                    if isempty(varargin)
                        
normal_line_X(1,:)=S_normal(1,1)*t+x_int;
normal_line_Y (1,:)=S_normal(1,2)*t+y_int;
normal_line_Z (1,:)=S_normal(1,3)*t+z_int;

% plot3(normal_line_X(1,:),normal_line_Y(1,:),normal_line_Z(1,:),'y:');

 fprintf(fileID_Normals,' %12.8f    %12.8f   %12.8f  \r\n ',normal_line_X(1,1),normal_line_Y(1,1),normal_line_Z(1,1));
  fprintf(fileID_Normals,' %12.8f    %12.8f   %12.8f  \r\n ',normal_line_X(1,2),normal_line_Y(1,2),normal_line_Z(1,2));
    fprintf(fileID_Normals,' nan nan nan \r\n');
                    end

    end


else
%     R_cyl(1) ==0
    % reflection from flat surface - Placement

    t=linspace(0,z_cyl_end/10,2);

       S_normal(1,1:3)=[0 0 1];
       
if isempty(varargin)
    
normal_line_X(1,:)=S_normal(1,1)*t+x_int;
normal_line_Y (1,:)=S_normal(1,2)*t+y_int;
normal_line_Z (1,:)=S_normal(1,3)*t+z_int;

% plot3(normal_line_X(1,:),normal_line_Y(1,:),[z_int z_int],'y:');

 fprintf(fileID_Normals,' %12.8f    %12.8f   %12.8f  \r\n ',normal_line_X(1,1),normal_line_Y(1,1),normal_line_Z(1,1));
  fprintf(fileID_Normals,' %12.8f    %12.8f   %12.8f  \r\n ',normal_line_X(1,2),normal_line_Y(1,2),normal_line_Z(1,2));
    fprintf(fileID_Normals,' nan nan nan \r\n');
end


end





%%  calculating reflection from normal and incoming ray!

Ray=([Rx Ry Rz]);

kr1=zeros(1,3);


% if S_normal(1,1)

kr1(1,1:3)=2*(-Ray*S_normal(1,:)')*S_normal(1,:)-(-Ray);
%kr1 direction should be outward of the surface !

if isempty(varargin)
Ref_line_X(1,:)=kr1(1,1)*t+x_int(1);
Ref_line_Y (1,:)=kr1(1,2)*t+y_int(1);
Ref_line_Z (1,:)=kr1(1,3)*t+z_int(1);
% plot3(Ref_line_X(1,:),Ref_line_Y(1,:),Ref_line_Z(1,:),'r--');



 fprintf(fileID_Refs,' %12.8f    %12.8f   %12.8f  \r\n ',Ref_line_X(1,1),Ref_line_Y(1,1),Ref_line_Z(1,1));
  fprintf(fileID_Refs,' %12.8f    %12.8f   %12.8f  \r\n ',Ref_line_X(1,2),Ref_line_Y(1,2),Ref_line_Z(1,2));
    fprintf(fileID_Refs,' nan nan nan \r\n');


% text(Ref_line_X(1,2),Ref_line_Y(1,2),Ref_line_Z(1,2) , sprintf('Ref'));

% Ref_line=[Ref_line_X;Ref_line_Y;Ref_line_Z];  % reflection line from the surface
end

end



kr1=kr1';  % just for output in general program



beta=0.5*acosd(dot(laser_direction,kr1)/(norm(laser_direction)*norm(kr1)));

% \beta is the angle between reflection and laser for measuring
%                 angle of incident
% 0.5* is because the angle between incident ray and reflection is twice as
% angle between normal and incident ray
