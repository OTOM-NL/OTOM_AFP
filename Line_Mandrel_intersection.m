

function [xyz_int]=Line_Mandrel_intersection(laser_source_P,laser_direction,R_cyl,z_cyl_end)



if R_cyl(1) ~=0


%first with cylinderical part
%second ellipsoid parts

if length(R_cyl)==3
    
    c1=R_cyl(2);
    c2=R_cyl(3);
    
else
    c1=R_cyl(1);
    c2=R_cyl(1);
    %       R_cyl=R_cyl(1);
    
end


xyz_int=[];

% kr1 reflection direction !

Lx=laser_source_P(1);
Ly=laser_source_P(2);
Lz=laser_source_P(3);

Rx=laser_direction(1);
Ry=laser_direction(2);
Rz= laser_direction(3);



%  Rxyz was normalised in General_3D_Optical

if  Rx==0
    Rx=Rx+ 0.00001;
end

if  Ry==0
    Ry=Ry+ 0.00001;
end

if  Rz==0
    Rz=Rz+ 0.00001;
end

norm_Rxyz=norm([Rx,Ry,Rz]);

Rx=Rx/norm_Rxyz;
Ry=Ry/norm_Rxyz;
Rz= Rz/norm_Rxyz;



Rxyz=[Rx;Ry;Rz];



%% Find intersection with cylinder

%assume that the center of cylinder in axis center

y0=Ly -(Ry/Rx)*Lx;
x_int=[];
y_int=[];
z_int=[];

% f=@(x) x^2 -R^2 + ((Ry/Rx)*x + y0 )^2;  % for Matlab 2015 or later
%  syms x ;% for old Matlab version 2013
% f= x^2 -R_cyl^2 + ((Ry/Rx)*x + y0 )^2;    % for old Matlab version 2013
% x_int=double(solve(f,'Real', true));



x_int=zeros(2,1);
x_int(1)= -(Rx*(Ry*y0 + (R_cyl(1)^2*Rx^2 + R_cyl(1)^2*Ry^2 - Rx^2*y0^2)^(1/2)))/(Rx^2 + Ry^2);
x_int(2)= -(Rx*(Ry*y0 - (R_cyl(1)^2*Rx^2 + R_cyl(1)^2*Ry^2 - Rx^2*y0^2)^(1/2)))/(Rx^2 + Ry^2);

if imag(x_int)
    x_int=[];
end
if isnan(x_int)
    x_int=[];
end


%if we have correct real results in 2D x-y plane, we continue the
%computations, otherwise we stop, because there is no intersection.

for ii=1:length(x_int)
    
    %     ii
    %     length(x_int)
    %     x_int
    % if x_int(ii)
    
    y_int(ii)=(Ry/Rx)*(x_int(ii))+ y0;
    z_int(ii)=((x_int(ii)-Lx)/Rx)*Rz + Lz;
    % z_int=((y_int-Ly)/Ry)*Rz + Lz;
    
    % now check whether z_int in the region of z_cylinder or not?
    % if z_cylinder cover the z_int space is fine, otherwise
    
    
    if z_int(ii) <0 | z_int(ii) > z_cyl_end
        %         warning('There is no intersection between line and surface!')
        %         x_int=[];
        %         y_int=[];
        %         z_int=[];
        %        xyz_int=[x_int ;y_int ;z_int ];
        
    elseif z_int(ii) >= 0 & z_int(ii) <=z_cyl_end
        
        
        
        n_xyz=[2*x_int(ii) 2*y_int(ii) 0];
        n_xyz=n_xyz/norm(n_xyz);
        
        
        theta_Fn=acosd(n_xyz*Rxyz);
        
        % check which intersection are the right by the help of gradient
        % vector
        
        if abs(theta_Fn)>=90
            
            % for sure only there is one point with these
            % conditions because the surface is convex
            
            %                      xyz_int=[x_int(ii) ;y_int(ii) ;z_int(ii) ];
            %                         plot3( nx1_P, ny1_P, nz1_P,'k--');
            % return
            
            
            xyz_int=[x_int(ii) ;y_int(ii) ;z_int(ii) ];
            
            V_AB=xyz_int- laser_source_P;
            V_AB=V_AB/norm(V_AB);
            
            
            
            theta_V_AB_Line3D=acosd(V_AB'*laser_direction );
            
            if abs(theta_V_AB_Line3D)<=1e-1  % parallel
                
                xyz_int=[x_int(ii) ;y_int(ii) ;z_int(ii) ];
                return;
                
            else
                xyz_int=[];
                
            end
            
            
            
            
            
            
            
            
        end
        
        
    end
    
    
    
    %     end
    % else
    %     warning('There is no intersection between line and surface!')
    
end

%%%%%%%%%%%%%%%% if the program did not "return", so there is no intersection !!!!
%      x_int=[];
%       y_int=[];
%         z_int=[];


%%

% Find intersection with 3D ellipsoidic parts -1

% Somw Part  tv > z_cyl_end  , c2

a=R_cyl(1);
b=R_cyl(1);
%             c1=R_cyl(1);

%             if length(R_cyl)==3
%                 c1=R_cyl(2);
%             end
%




x0_L=Lx;
y0_L=Ly;
z0_L=Lz;

x0_E1=0;
y0_E1=0;
z0_E1=z_cyl_end;

x0=x0_L-x0_E1;
y0=y0_L-y0_E1;
z0=z0_L-z0_E1;


t1 = -(a*b*c2*(Rx^2*b^2*c2^2 - Rx^2*b^2*z0^2 - Rx^2*c2^2*y0^2 + 2*Rx*Ry*c2^2*x0*y0 + 2*Rx*Rz*b^2*x0*z0 + Ry^2*a^2*c2^2 - Ry^2*a^2*z0^2 - Ry^2*c2^2*x0^2 + 2*Ry*Rz*a^2*y0*z0 + Rz^2*a^2*b^2 - Rz^2*a^2*y0^2 - Rz^2*b^2*x0^2)^(1/2) + Rx*b^2*c2^2*x0 + Ry*a^2*c2^2*y0 + Rz*a^2*b^2*z0)/(Rx^2*b^2*c2^2 + Ry^2*a^2*c2^2 + Rz^2*a^2*b^2);
t2=-(Rx*b^2*c2^2*x0 - a*b*c2*(Rx^2*b^2*c2^2 - Rx^2*b^2*z0^2 - Rx^2*c2^2*y0^2 + 2*Rx*Ry*c2^2*x0*y0 + 2*Rx*Rz*b^2*x0*z0 + Ry^2*a^2*c2^2 - Ry^2*a^2*z0^2 - Ry^2*c2^2*x0^2 + 2*Ry*Rz*a^2*y0*z0 + Rz^2*a^2*b^2 - Rz^2*a^2*y0^2 - Rz^2*b^2*x0^2)^(1/2) + Ry*a^2*c2^2*y0 + Rz*a^2*b^2*z0)/(Rx^2*b^2*c2^2 + Ry^2*a^2*c2^2 + Rz^2*a^2*b^2);



%%

if ~ imag(t1)
    xp1=t1*Rx+x0_L;
    yp1=t1*Ry+y0_L;
    zp1=t1*Rz+z0_L;
    
    % for upper side
    %                 zp1 >=z0_E &&    zp1 <= z0_E+c2
    
    % for bottom side
    % zp1 <=z0_E &&    zp1 >= z0_E-c;
    
    if zp1 >=z0_E1 &&    zp1 <= z0_E1+c2;
        
        %gradiant calculations
        nx1=(xp1-x0_E1)/a^2;
        ny1=(yp1-y0_E1)/b^2;
        nz1=(zp1-z0_E1)/c2^2;
        
        Norm_nxyz= norm([nx1,ny1,nz1]);
        nx1=nx1/Norm_nxyz;
        ny1=ny1/Norm_nxyz;
        nz1=nz1/Norm_nxyz;
        n_xyz1=[nx1 ny1 nz1];
        
        
        theta_Fn1=acosd(n_xyz1*Rxyz);
        
        if abs(theta_Fn1)>=90
            
            %                         plot3( nx1_P, ny1_P, nz1_P,'k--');
            xyz_int=[xp1 ;yp1 ;zp1 ];
            
            V_AB=xyz_int- laser_source_P;
            V_AB=V_AB/norm(V_AB);
            
            
            
            theta_V_AB_Line3D=acosd(V_AB'*laser_direction );
            
            if abs(theta_V_AB_Line3D)<=1e-1  % parallel
                
                xyz_int=[xp1 ;yp1 ;zp1 ];
                return;
                
            else
                xyz_int=[];
                
            end
            
        end
        
        
    end
    
end

% for second point


if ~ imag(t2)
    xp2=t2*Rx+x0_L;
    yp2=t2*Ry+y0_L;
    zp2=t2*Rz+z0_L;
    
    if zp2 >=z0_E1 &&    zp2 <= z0_E1+c2;
        
        %gradiant calculations
        nx2=(xp2-x0_E1)/a^2;
        ny2=(yp2-y0_E1)/b^2;
        nz2=(zp2-z0_E1)/c2^2;
        
        Norm_nxyz= norm([nx2,ny2,nz2]);
        nx2=nx2/Norm_nxyz;
        ny2=ny2/Norm_nxyz;
        nz2=nz2/Norm_nxyz;
        n_xyz2=[nx2 ny2 nz2];
        
        
        theta_Fn2=acosd(n_xyz2*Rxyz);
        
        if abs(theta_Fn2)>=90
            
            
            %                          plot3( nx1_P, ny1_P, nz1_P,'k--');
            xyz_int=[xp2 ;yp2 ;zp2 ];
            
            V_AB=xyz_int- laser_source_P;
            V_AB=V_AB/norm(V_AB);
            
            
            
            theta_V_AB_Line3D=acosd(V_AB'*laser_direction);
            
            if abs(theta_V_AB_Line3D)<=1e-1  % parallel
                
                xyz_int=[xp2 ;yp2 ;zp2 ];
                return;
                
            else
                xyz_int=[];
                
            end
            
            
            %                         plot3( nx1_P, ny1_P, nz1_P,'k--');
            %                            xyz_int=[xp2 ;yp2 ;zp2 ];
            %                             return;
            
        end
        
        
    end
    
    
    %                           text(xp2,yp2,zp2,sprintf('P2=%f, %f, %f',xp1,yp1,zp1))
    %                     plot3(xp2,yp2,zp2,'r*');
end




%% Lower Dome 1st Dome


% Find intersection with 3D ellipsoidic parts -2



a=R_cyl(1);
b=R_cyl(1);
%             c2=R_cyl(1);

%              if length(R_cyl)==3
%                 c2=R_cyl(2);
%             end



x0_E2=0;
y0_E2=0;
z0_E2=0;

x0=x0_L-x0_E2;
y0=y0_L-y0_E2;
z0=z0_L-z0_E2;


t1 = -(a*b*c1*(Rx^2*b^2*c1^2 - Rx^2*b^2*z0^2 - Rx^2*c1^2*y0^2 + 2*Rx*Ry*c1^2*x0*y0 + 2*Rx*Rz*b^2*x0*z0 + Ry^2*a^2*c1^2 - Ry^2*a^2*z0^2 - Ry^2*c1^2*x0^2 + 2*Ry*Rz*a^2*y0*z0 + Rz^2*a^2*b^2 - Rz^2*a^2*y0^2 - Rz^2*b^2*x0^2)^(1/2) + Rx*b^2*c1^2*x0 + Ry*a^2*c1^2*y0 + Rz*a^2*b^2*z0)/(Rx^2*b^2*c1^2 + Ry^2*a^2*c1^2 + Rz^2*a^2*b^2);
t2=-(Rx*b^2*c1^2*x0 - a*b*c1*(Rx^2*b^2*c1^2 - Rx^2*b^2*z0^2 - Rx^2*c1^2*y0^2 + 2*Rx*Ry*c1^2*x0*y0 + 2*Rx*Rz*b^2*x0*z0 + Ry^2*a^2*c1^2 - Ry^2*a^2*z0^2 - Ry^2*c1^2*x0^2 + 2*Ry*Rz*a^2*y0*z0 + Rz^2*a^2*b^2 - Rz^2*a^2*y0^2 - Rz^2*b^2*x0^2)^(1/2) + Ry*a^2*c1^2*y0 + Rz*a^2*b^2*z0)/(Rx^2*b^2*c1^2 + Ry^2*a^2*c1^2 + Rz^2*a^2*b^2);



%%

if ~ imag(t1)
    xp1=t1*Rx+x0_L;
    yp1=t1*Ry+y0_L;
    zp1=t1*Rz+z0_L;
    
    % for upper side
    %                 zp1 >=z0_E &&    zp1 <= z0_E+c1
    
    % for bottom side
    % zp1 <=z0_E &&    zp1 >= z0_E-c1;
    
    if zp1 <=z0_E2 &&    zp1 >= z0_E2-c1;
        
        %gradiant calculations
        nx1=(xp1-x0_E2)/a^2;
        ny1=(yp1-y0_E2)/b^2;
        nz1=(zp1-z0_E2)/c1^2;
        
        Norm_nxyz= norm([nx1,ny1,nz1]);
        nx1=nx1/Norm_nxyz;
        ny1=ny1/Norm_nxyz;
        nz1=nz1/Norm_nxyz;
        n_xyz1=[nx1 ny1 nz1];
        
        
        
        %%
        
        
        theta_Fn1=acosd(n_xyz1*Rxyz);
        
        %                     if abs(theta_Fn1)>=90
        %
        % %                         plot3( nx1_P, ny1_P, nz1_P,'k--');
        %                            xyz_int=[xp1 ;yp1 ;zp1 ];
        %                             return;
        %
        %                     end
        
        
        
        if abs(theta_Fn1)>=90
            
            
            %                          plot3( nx1_P, ny1_P, nz1_P,'k--');
            xyz_int=[xp1 ;yp1 ;zp1 ];
            
            V_AB=xyz_int- laser_source_P;
            V_AB=V_AB/norm(V_AB);
            
            
            
            theta_V_AB_Line3D=acosd(V_AB'*laser_direction);
            
            if abs(theta_V_AB_Line3D)<=1e-1  % parallel
                
                xyz_int=[xp1 ;yp1 ;zp1 ];
                return;
                
            else
                xyz_int=[];
                
            end
            
            
            %                         plot3( nx1_P, ny1_P, nz1_P,'k--');
            %                            xyz_int=[xp2 ;yp2 ;zp2 ];
            %                             return;
            
        end
        
        
        
        
        %%
        
    end
    
end

% for second point


if ~ imag(t2)
    xp2=t2*Rx+x0_L;
    yp2=t2*Ry+y0_L;
    zp2=t2*Rz+z0_L;
    
    if zp2 <=z0_E2 &&    zp2 >= z0_E2-c1;
        
        %gradiant calculations
        nx2=(xp2-x0_E2)/a^2;
        ny2=(yp2-y0_E2)/b^2;
        nz2=(zp2-z0_E2)/c1^2;
        
        Norm_nxyz= norm([nx2,ny2,nz2]);
        nx2=nx2/Norm_nxyz;
        ny2=ny2/Norm_nxyz;
        nz2=nz2/Norm_nxyz;
        n_xyz2=[nx2 ny2 nz2];
        
        
        theta_Fn2=acosd(n_xyz2*Rxyz);
        
        %                     if abs(theta_Fn2)>=90
        %
        % %                         plot3( nx1_P, ny1_P, nz1_P,'k--');
        %                            xyz_int=[xp2 ;yp2 ;zp2 ];
        %                             return;
        %
        %                     end
        
        
        if abs(theta_Fn2)>=90
            
            
            %                          plot3( nx1_P, ny1_P, nz1_P,'k--');
            xyz_int=[xp2 ;yp2 ;zp2 ];
            
            V_AB=xyz_int- laser_source_P;
            V_AB=V_AB/norm(V_AB);
            
            
            
            theta_V_AB_Line3D=acosd(V_AB'*laser_direction);
            
            if abs(theta_V_AB_Line3D)<=1e-1  % parallel
                
                xyz_int=[xp2 ;yp2 ;zp2 ];
                return;
                
            else
                xyz_int=[];
                
            end
            
            
            %                         plot3( nx1_P, ny1_P, nz1_P,'k--');
            %                            xyz_int=[xp2 ;yp2 ;zp2 ];
            %                             return;
            
        end
        
        
        
        
        
        
    end
    
    
    %                           text(xp2,yp2,zp2,sprintf('P2=%f, %f, %f',xp1,yp1,zp1))
    %                     plot3(xp2,yp2,zp2,'r*');
end









%% mandrel
% [x_E,y_E,z_E]=sphere(30);
% x_E=(x_E*a)+x0_E;
% y_E=y_E*b+y0_E;
% z_E=z_E*c1+z0_E;
%
% % consider the semi-ellipsoid- Upper bound
% % Has been devided based on the z-direction!
%
% [m,n]=size(z_E);
%
% start_ind=ceil(m/2);
%
% surf_val=zeros(size(x_E));
%
%
% % upper side
% % surf(x_E(start_ind:end,:),y_E(start_ind:end,:),z_E(start_ind:end,:), surf_val(start_ind:end,:));
%
% % bottom side
% surf(x_E(1:start_ind,:),y_E(1:start_ind,:),z_E(1:start_ind,:), surf_val(1:start_ind,:));
%
% axis equal;
% camlight right;
%          axis fill;


else 
%    R_cyl(1) ==0 ;
   
   
   xyz_int=Placement_int (laser_source_P,laser_direction);
    
end



