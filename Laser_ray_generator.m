
function laser_point_Actual_all=Laser_ray_generator (Rx,Ry,Rz,Ax,Ay,nx,ny,L_xyz0,ID,Laser_head_Rot,Divergence_factor)


% Divergence_factor=2.4;

%  fid23 = fopen('.\Supp_files\Laser_characteristic.txt','r');
%
%
%  out = textscan(fid23,'%s ','delimiter',',');
%  fclose(fid23);
%
% Divergence_factor=str2double(out{1}{2});

% ID=laser energy distribution definition
% Ray distribution is uniform
% an ID which define the laser distribution
% ID=0 uniform, 1 for tophat, 2 for Gaussian, 3 for linear

% SAme input energy !
% compute the amplitude
area=2*Ax*2*Ay;

Gauss_Par_X=ID(2);
Gauss_Par_Y=ID(3);




% ID=[2,Gauss_Par_X,Gauss_Par_Y]

switch ID(1)
    
    case 0  % uniform
        fun=@(x,y) 1;
        amp=1;
    case 1  % linear
         x0=ID(6);
        y0=ID(7);
        
          fun=@(x,y)  1+(Gauss_Par_X.*x+x0 )+(Gauss_Par_Y.*y+y0 );
        amp=1;
        
        
    case  2 % Gaussian
        fun=@(x,y) exp(-((x*Gauss_Par_X).^2+(y*Gauss_Par_Y).^2));
        amp=area/integral2(fun,-Ax,Ax,-Ay,Ay);
    case 3
        
        Gauss_Par_X2=ID(4);
        Gauss_Par_Y2=ID(5);
        
        x0=ID(6);
        y0=ID(7);
        
        
        fun1=@(x,y) exp(-(( (x-x0).*Gauss_Par_X).^2+(  (y-y0).*Gauss_Par_Y).^2)) ;
        
        fun2=@(x,y) exp(-(((x-x0).*Gauss_Par_X2).^2+((y-y0).*Gauss_Par_Y).^2));
        
        fun3=@(x,y) exp(-(((x-x0)*Gauss_Par_X2).^2+((y-y0).*Gauss_Par_Y2).^2));
        
        fun4=@(x,y) exp(-(((x-x0)*Gauss_Par_X).^2+((y-y0).*Gauss_Par_Y2).^2));
        
end




v1=[0 0 1]';   % the initial laser head on xy plane
v2=[Rx Ry Rz]';
Rot_plane=transpose(fcn_RotationFromTwoVectors(v1, v2)); % find the rotation matrix

%  r=vrrotvec(v1,v2);
% m = vrrotvec2mat(r);

%  Ax=5;
%  Ay=2;
%  nx=100;
%  ny=40;


%%
%  theta_ray=acosd(v2(3)/norm(v2));
% phi_ray=atan2d(v2(2),v2(1));

% delta_phi=-10;%*sign(phi_ray); %deg minus for this value
% delta_theta=+0;%*sign(theta_ray); %deg + positive for this value

%%


xp_laser=linspace(-Ax,Ax,nx);
yp_laser=linspace(-Ay,Ay,ny);

[X,Y]=meshgrid(xp_laser,yp_laser);  % ny * nx = size
% Z=zeros(size(X));
% figure(20);
%  surf(X,Y,Z);
%  hold on;




X_Actual=zeros(size(X));
Y_Actual=zeros(size(X));
Z_Actual=zeros(size(X));
Power_Actual=zeros(size(X));

% laser_point_Actual=zeros(1,3);
laser_point_Actual_all=zeros(nx*ny,7);


% theta=atan(Ry/Rx);
% theta= atan(Rz/Rx)+atan(Ry/Rx);
% theta= atan( -0.504360 /-0.863494 )+pi/6;
theta=Laser_head_Rot*(pi/180);

% rot_error=1;

% in each loop we modify the rotate laser head
% while abs(rot_error) > 1e-3

counter=0;

if ID(1) ==3
    
    
    amp1=area/integral2(fun1,x0,Ax,y0,Ay);
    amp2=area/integral2(fun2,-Ax,x0,y0,Ay);
    amp3=area/integral2(fun3,-Ax,x0,-Ay,y0);
    amp4=area/integral2(fun4,x0,Ax,-Ay,y0);
    
    amp=(amp1+amp2+amp3+amp4)/16;
    
    
    for ii=1:length(xp_laser)
        for jj=1:length(yp_laser)
            counter=counter+1;
            
            new_points=Transformation_Rot_z([0 0 0],-theta,[xp_laser(ii);yp_laser(jj); 0]);
            laser_point_Actual= [new_points(1),new_points(2), new_points(3)]* (Rot_plane);
            %         laser_point_Actual= [xp_laser(ii),yp_laser(jj), 0]* (Rot_plane);
            X_Actual(jj,ii) =laser_point_Actual(1);
            Y_Actual(jj,ii) =laser_point_Actual(2);
            Z_Actual(jj,ii) =laser_point_Actual(3);
            laser_point_Actual_all(counter,1:3)=laser_point_Actual;
            
            
            if xp_laser(ii)>=x0 && yp_laser(jj) >=y0
                
                
                laser_point_Actual_all(counter,4)=amp* fun1(xp_laser(ii),yp_laser(jj));
                
                
            elseif xp_laser(ii)<=x0 && yp_laser(jj) >=y0
                laser_point_Actual_all(counter,4)=amp* fun2(xp_laser(ii),yp_laser(jj));
                
            elseif xp_laser(ii)<=x0 && yp_laser(jj) <=y0
                
                laser_point_Actual_all(counter,4)=amp* fun3(xp_laser(ii),yp_laser(jj));
            elseif xp_laser(ii)>=x0 && yp_laser(jj) <=y0
                laser_point_Actual_all(counter,4)=amp* fun4(xp_laser(ii),yp_laser(jj));
            end
            
            Power_Actual(jj,ii)= laser_point_Actual_all(counter,4);
            
            Vec_C_point_actual=laser_point_Actual; %-0
             b2=v2'+ (Vec_C_point_actual*Divergence_factor);
             b2=b2/norm(b2);
             laser_point_Actual_all(counter,5:7)=b2;
            
            
            
        end
    end
    
    
    
    
else
    
    
    
    for ii=1:length(xp_laser)
        for jj=1:length(yp_laser)
            counter=counter+1;
            
            new_points=Transformation_Rot_z([0 0 0],-theta,[xp_laser(ii);yp_laser(jj); 0]);
            laser_point_Actual= [new_points(1),new_points(2), new_points(3)]* (Rot_plane);
            %         laser_point_Actual= [xp_laser(ii),yp_laser(jj), 0]* (Rot_plane);
            X_Actual(jj,ii) =laser_point_Actual(1);
            Y_Actual(jj,ii) =laser_point_Actual(2);
            Z_Actual(jj,ii) =laser_point_Actual(3);
            laser_point_Actual_all(counter,1:3)=laser_point_Actual;
            
            
            
            
            
            
            %energy will be specified here
            laser_point_Actual_all(counter,4)=amp*fun(xp_laser(ii),yp_laser(jj));
            Power_Actual(jj,ii)= laser_point_Actual_all(counter,4);
            
            %% Calculate \theta and \phi for each location
            %                 sign(xp_laser(ii))
            %                 sign(yp_laser(jj))clc
            % (xp_laser(ii)/Ax)
            %                 (yp_laser(jj)/Ay)
            
            Vec_C_point_actual=laser_point_Actual; %-0
            % plot3(laser_point_Actual(1),laser_point_Actual(2),laser_point_Actual(3),'r*');
            
            b2=v2'+ (Vec_C_point_actual*Divergence_factor);
            
            %                 Rx2=cosd(phi_ray +  (xp_laser(ii)/Ax)* delta_phi)*sind(theta_ray+ (yp_laser(jj)/Ay)*   delta_theta);
            %                 Ry2=sind(phi_ray + (xp_laser(ii)/Ax)*   delta_phi)*sind(theta_ray+ (yp_laser(jj)/Ay)* delta_theta);
            %                 Rz2=cosd(theta_ray+ (yp_laser(jj)/Ay)* delta_theta);
            
            %                 b2=[Rx2 Ry2 Rz2];
            b2=b2/norm(b2);
            
            laser_point_Actual_all(counter,5:7)=b2;
            %                 laser_point_Actual_all(counter,5:7)=[Rx2 Ry2 Rz2]/norm([Rx2 Ry2 Rz2]);
            %%
            %                 if ii==1 & jj==1
            %                     max_divergence=acosd(b2*v2);
            %                 end
            
            
            
        end
    end
    
    
end
% this formulation has the error and should be compensated after a loop

% norm_tri=norm([X_Actual(1,1)-X_Actual(1,end),Y_Actual(1,1)-Y_Actual(1,end),Z_Actual(1,1)-Z_Actual(1,end)]);
% rot_error=asin(  (Y_Actual(1,1)-Y_Actual(1,end)) / norm_tri);
%
%
% theta=theta+rot_error;
%
%
% end


%% Find the factor based on the Given Divergence!
% theta divergence in degree
% theta_div=5;
% F_Div=@ (Divergence_factor) acosd((((Divergence_factor*Vec_C_point_actual)+v2')*v2)/norm((Divergence_factor*Vec_C_point_actual)+v2'))-theta_div;
%
%  x0 = [0 100]; % initial interval
% Divergence_factor = fzero(F_Div,x0);



%  xp=[0,Rx];
%   yp=[0,Ry];
%    zp=[0,Rz];
%    plot3(xp,yp,zp);



X_Actual=X_Actual+L_xyz0(1);
Y_Actual=Y_Actual+L_xyz0(2);
Z_Actual=Z_Actual+L_xyz0(3);

laser_point_Actual_all(:,4)=(laser_point_Actual_all(:,4)/ sum(Power_Actual(:)))*nx*ny;

Power_Actual=(Power_Actual/ sum(Power_Actual(:)))*nx*ny;



hold on;
surf(X_Actual,Y_Actual,Z_Actual,Power_Actual,'LineStyle','--', 'FaceColor',[0.3 0.75 0.9]);

%   hold on;
% figure;
%  surf(X_Actual,Y_Actual,Z_Actual,Power_Actual); % for general figure (1)
%  title('Actual laser head distribution');
%  colorbar;


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

