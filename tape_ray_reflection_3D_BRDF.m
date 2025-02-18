
%normal and reflection for Tape -Ray 3D
function [kr1_active_G,beta_macro]=tape_ray_reflection_3D_BRDF (R,S1,S2,xyz_int_L,x0,y0,laser_direction,tv,Rot_Roller_axis,fileID_Normals,fileID_Refs,xyz_int_G,BRDF_Tape)


% BRDF_Tape=[Fib_Rot,div,sig_T,sig_F,Amp,threshold]


% Fib_Rot= 90; % in degree



% Rx=laser_direction(1);
% Ry=laser_direction(2);
% Rz= laser_direction(3);

x_int=xyz_int_L(1);
y_int=xyz_int_L(2);
% z_int=xyz_int_L(3);



t=linspace(0,R/5,2);
% % kr1=zeros(1,3);

if S1
    
    norm_n=norm([x_int y_int-R]); % the circle is not located at center, at (x-0)^2 +(y-R)^2=R^2
    
    S_normal=([x_int y_int-R ])/norm_n;
    
    normal_line_X=S_normal(1)*t+x_int;
    normal_line_Y =S_normal(2)*t+y_int;
    
    %     plot3(normal_line_X,normal_line_Y,[z_int z_int],'k:')
    %
    %     text(normal_line_X(2),normal_line_Y(2),z_int , sprintf('n'));
    %     hold on
    
    
    
    
elseif S2
    norm_n=norm([x0 y0-R]);
    
    S_normal=([x0 y0-R ])/norm_n;
    
    normal_line_X=S_normal(1)*t+x_int;
    normal_line_Y =S_normal(2)*t+y_int;
    
    %     plot3(normal_line_X,normal_line_Y,[z_int z_int],'k:')
    %
    %     text(normal_line_X(2),normal_line_Y(2),z_int , sprintf('n'));
    %     hold on
    
    % else
    %     disp('no intersection, no normal vector')
    
end

%% reflection from BRDF

% All in local axis

% Fib_Rot= 90; % in degree

n_macro_L=[S_normal 0];
% n_macro_G=Transformation_L2G (zeros(3,1),Rot_Roller_axis,n_macro_L');

% figure(1);
%
% % t=linspace(0,0.5,2);
%   Pr_N(1,:)=n_macro_G(1)*t+xyz_int_G(1);
%                 Pr_N (2,:)=n_macro_G(2)*t+xyz_int_G(2);
%                      Pr_N (3,:)=n_macro_G(3)*t+xyz_int_G(3);
%                                   plot3(Pr_N(1,:),Pr_N(2,:),Pr_N(3,:),'k--');




[~,beta_macro,S_normal_active,kr1_active]=Rot_Matrix_Finder_local_BRDF_Tape(BRDF_Tape,n_macro_L,laser_direction);

%  v1=[0 ;0; 1];
%  n_normal_T=Rot_Roller_axis*Rot_Roller_axis_L*v1;


[~,n]=size(S_normal_active);

for ii=1:n
    
    Dir_G= Rot_Roller_axis*S_normal_active(1:3,ii);
    
    normal_line_X(1,:)=Dir_G(1)*t+xyz_int_G(1);
    normal_line_Y (1,:)=Dir_G(2)*t+xyz_int_G(2);
    normal_line_Z (1,:)=Dir_G(3)*t+xyz_int_G(3);
    
    
    
    
    
    fprintf(fileID_Normals,' %12.8f    %12.8f   %12.8f  \r\n ',normal_line_X(1,1),normal_line_Y(1,1),normal_line_Z(1,1));
    fprintf(fileID_Normals,' %12.8f    %12.8f   %12.8f  \r\n ',normal_line_X(1,2),normal_line_Y(1,2),normal_line_Z(1,2));
    fprintf(fileID_Normals,' nan nan nan \r\n');
    
    
    Dir_ref_G= Rot_Roller_axis*kr1_active(1:3,ii);
    
    Ref_line_X(1,:)=Dir_ref_G(1)*t+xyz_int_G(1);
    Ref_line_Y (1,:)=Dir_ref_G(2)*t+xyz_int_G(2);
    Ref_line_Z (1,:)=Dir_ref_G(3)*t+xyz_int_G(3);
    
    
    
    fprintf(fileID_Refs,' %12.8f    %12.8f   %12.8f  \r\n ',Ref_line_X(1,1),Ref_line_Y(1,1),Ref_line_Z(1,1));
    fprintf(fileID_Refs,' %12.8f    %12.8f   %12.8f  \r\n ',Ref_line_X(1,2),Ref_line_Y(1,2),Ref_line_Z(1,2));
    fprintf(fileID_Refs,' nan nan nan \r\n');
    
    
    
end






%
%            if isempty(varargin)






kr1_active_G= Rot_Roller_axis*kr1_active(1:3,:);

kr1_active_G(4,:)= kr1_active(4,:);















%%
% if S1 | S2
%      tR=linspace(0,R,2);
%     Ray=([Rx Ry Rz]);
%
%        S_normal(1,1:3)=[S_normal(1,:) 0 ];
%         Ray=-Ray;
%        kr1(1,1:3)=2*(Ray*S_normal(1,:)')*S_normal(1,:)-Ray;
% Ref_line_X(1,:)=kr1(1,1)*tR+x_int;
% Ref_line_Y (1,:)=kr1(1,2)*tR+y_int;
% Ref_line_Z (1,:)=kr1(1,3)*tR+z_int;
%
% %%
% % transform to global coordinate
% %%
%
% %%
% points(1,:)=normal_line_X;
% points(2,:)=normal_line_Y;
% points(3,:)=([z_int z_int]);
%
% new_points=Transformation_L2G (tv,Rot_Roller_axis,points);
% normal_line_X=new_points(1,:);
% normal_line_Y=new_points(2,:);
% z_int=new_points(3,1);
%
%
%
% %  plot3(normal_line_X,normal_line_Y,[z_int z_int],'y:');
%
% if isempty(varargin)
%
%
%
%    fprintf(fileID_Normals,' %12.8f    %12.8f   %12.8f  \r\n ',normal_line_X(1,1),normal_line_Y(1,1),z_int);
%   fprintf(fileID_Normals,' %12.8f    %12.8f   %12.8f  \r\n ',normal_line_X(1,2),normal_line_Y(1,2),z_int);
%     fprintf(fileID_Normals,' nan nan nan \r\n');
%
%
% %       text(normal_line_X(2),normal_line_Y(2),z_int , sprintf('n'));
%     hold on
% %%
%
%
% points(1,:)=Ref_line_X;
% points(2,:)=Ref_line_Y;
% points(3,:)=Ref_line_Z;
%
% new_points=Transformation_L2G (tv,Rot_Roller_axis,points);
% Ref_line_X=new_points(1,:);
% Ref_line_Y=new_points(2,:);
% Ref_line_Z=new_points(3,:);
%
% % plot3(Ref_line_X(1,:),Ref_line_Y(1,:),Ref_line_Z(1,:),'b--')
%
%  fprintf(fileID_Refs,' %12.8f    %12.8f   %12.8f  \r\n ',Ref_line_X(1,1),Ref_line_Y(1,1),Ref_line_Z(1,1));
%   fprintf(fileID_Refs,' %12.8f    %12.8f   %12.8f  \r\n ',Ref_line_X(1,2),Ref_line_Y(1,2),Ref_line_Z(1,2));
%     fprintf(fileID_Refs,' nan nan nan \r\n');
%
% end
% text(Ref_line_X(1,2),Ref_line_Y(1,2),Ref_line_Z(1,2) , sprintf('Ref'));


% kr1=kr1';  % just for match in general program !

% beta=0.5*acosd(dot(laser_direction,kr1)/(norm(laser_direction)*norm(kr1)));
% \beta is the angle between normal and laser for measuring
%                 angle of incident
% 0.5* is because the angle between incident ray and reflection is twice as
% angle between normal and incident ray

end