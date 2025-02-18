
%% NOT USING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


function laser_point_Actual_all=Laser_ray_generator (Rx,Ry,Rz,Ax,Ay,nx,ny,L_xyz0,ID)

% ID=laser energy distribution definition
% Ray distribution is uniform
% an ID which define the laser distribution
% ID=0 uniform, 1 for tophat, 2 for Gaussian, 3 for linear

% SAme input energy !
% compute the amplitude
area=2*Ax*2*Ay;
switch ID   
    
    case 0  % uniform
         fun=@(x,y) 1;
        amp=1;
    case 1  % top-hat
        
    case  2 % Gaussian
        fun=@(x,y) exp(-(x.^2+y.^2));
        amp=area/integral2(fun,-Ax,Ax,-Ay,Ay);
end







v1=[0 0 1]';   % the initial laser head on xy plane
 v2=[Rx Ry Rz]';
 Rot_plane=transpose(fcn_RotationFromTwoVectors(v1, v2)); % find the rotation matrix
 
  
%  Ax=5;
%  Ay=2;
%  nx=100;
%  ny=40;
 xp_laser=linspace(-Ax,Ax,nx);
 yp_laser=linspace(-Ay,Ay,ny);
 
 
 [X,Y]=meshgrid(xp_laser,yp_laser);  % ny * nx = size
Z=zeros(size(X));
figure(20);
 surf(X,Y,Z);
 hold on;
 
 
 
 
X_Actual=zeros(size(X));
Y_Actual=zeros(size(X));
Z_Actual=zeros(size(X));
 Power_Actual=zeros(size(X));

laser_point_Actual=zeros(1,3);
laser_point_Actual_all=zeros(nx*ny,4);


%%
 
laser_point_Actual_Temp=zeros(1,3);
X_Actual_Temp=zeros(length(Y));
Y_Actual_Temp=zeros(length(Y));
Z_Actual_Temp=zeros(length(Y));
  
%   plot3(point_x(:,1),point_x(:,2),point_x(:,3),'b--');
  

% for ii=1:length(xp_laser)
for jj=1:length(yp_laser)
 
    laser_point_Actual_Temp= [0,yp_laser(jj), 0]* (Rot_plane);
    X_Actual_Temp(jj) =laser_point_Actual_Temp(1);
    Y_Actual_Temp(jj) =laser_point_Actual_Temp(2);
    Z_Actual_Temp(jj) =laser_point_Actual_Temp(3);
    %         laser_point_Actual_all(counter,1:3)=laser_point_Actual;
    
    
end
  
  plot3(X_Actual_Temp,Y_Actual_Temp,Z_Actual_Temp,'b--');


%    laser_point_Actual_Temp= [0 ,yp_laser(jj), 0]* (Rot_plane);  % to define direction
    
      C = cross( laser_point_Actual_Temp,v2);  % find normal of two vectors
  C=C/norm(C);   % to have the unit size 
  

    counter=0;
 for jj=1:length(yp_laser)
  for ii=1:length(xp_laser)
      counter=counter+1;
  laser_point_Actual=C*xp_laser(ii) + [ X_Actual_Temp(jj), Y_Actual_Temp(jj), Z_Actual_Temp(jj)];
 
  plot3(laser_point_Actual(1),laser_point_Actual(2),laser_point_Actual(3),'r*');
  
  
     X_Actual(jj,ii) =laser_point_Actual(1);
        Y_Actual(jj,ii) =laser_point_Actual(2);
        Z_Actual(jj,ii) =laser_point_Actual(3);
  
             laser_point_Actual_all(counter,1:3)=laser_point_Actual;
  
  
  %energy will be specified here
          laser_point_Actual_all(counter,4)=amp*fun(xp_laser(ii),yp_laser(jj));
                Power_Actual(jj,ii)= laser_point_Actual_all(counter,4);
  
  end
 end


  
    
    %%
    
    


% counter=0;
% for ii=1:length(xp_laser)
%     for jj=1:length(yp_laser)
%         counter=counter+1;
%         laser_point_Actual= [xp_laser(ii),yp_laser(jj), 0]* (Rot_plane);
%         X_Actual(jj,ii) =laser_point_Actual(1);
%         Y_Actual(jj,ii) =laser_point_Actual(2);
%         Z_Actual(jj,ii) =laser_point_Actual(3);
%         laser_point_Actual_all(counter,1:3)=laser_point_Actual;
%         
%         %energy will be specified here
%         laser_point_Actual_all(counter,4)=amp*fun(xp_laser(ii),yp_laser(jj));
%         
%         Power_Actual(jj,ii)= laser_point_Actual_all(counter,4);
%     end
% end

 


 %%
 xp=[0,Rx];
  yp=[0,Ry];
   zp=[0,Rz];
   plot3(xp,yp,zp);
   
    xp=[0,0];
  yp=[0,0];
   zp=[0,1];
   plot3(xp,yp,zp);
  surf(X_Actual,Y_Actual,Z_Actual,Power_Actual); % for general figure (1)
  
  
  
  
  %%
  

counter=0;
for ii=1:length(xp_laser)
    for jj=1:length(yp_laser)
        counter=counter+1;
        laser_point_Actual= [xp_laser(ii),yp_laser(jj), 0]* (Rot_plane);
        X_Actual(jj,ii) =laser_point_Actual(1);
        Y_Actual(jj,ii) =laser_point_Actual(2);
        Z_Actual(jj,ii) =laser_point_Actual(3);
        laser_point_Actual_all(counter,1:3)=laser_point_Actual;
        
        %energy will be specified here
        laser_point_Actual_all(counter,4)=amp*fun(xp_laser(ii),yp_laser(jj));
        
        Power_Actual(jj,ii)= laser_point_Actual_all(counter,4);
    end
end
  %%
  

   
%   C = cross( laser_point_Actual,v2);
%   C=C/norm(C);
%   
%   t=linspace(0,Ay/2,10);
%   for ii=1:10
%   point_x(ii,1:3)=C*t(ii);
%  
%   
%   end
%   
%   plot3(point_x(:,1),point_x(:,2),point_x(:,3),'b--');
%%
 
 
 X_Actual=X_Actual+L_xyz0(1);
  Y_Actual=Y_Actual+L_xyz0(2);
 Z_Actual=Z_Actual+L_xyz0(3);
 
   hold on;
   figure(1)
 surf(X_Actual,Y_Actual,Z_Actual,Power_Actual); 
 
%   hold on;
%  surf(X_Actual,Y_Actual,Z_Actual,Power_Actual); % for general figure (1)

 
%  figure (10);
% %  surf(X_Actual,Y_Actual,Z_Actual,Power_Actual);
% 
%  mesh(X,Y,Power_Actual);
%  axis equal;
%  title('Laser energy distribution');

laser_point_Actual_all(:,1:3)=laser_point_Actual_all(:,1:3)+(ones(nx*ny,3)*diag(L_xyz0));

% laser_point_Actual_all(:,1:3)=laser_point_Actual_all(:,1:3)+(L_xyz0);
%   
 % add all points to [x0 y0 z0]
 
