
% This function only for cylinder and NOT functioning anymore!!


function [xyz_int]=Line_cylinder_intersection(laser_source_P,laser_direction,R_cyl,z_cyl_end)

% kr1 reflection direction !

Lx=laser_source_P(1);
Ly=laser_source_P(2);
Lz=laser_source_P(3);

 Rx=laser_direction(1);
 Ry=laser_direction(2);
 Rz= laser_direction(3);
 
 if  Rx==0
     Rx=Rx+ 0.00001;
 end
 
  if  Ry==0
     Ry=Ry+ 0.00001;
  end
 
 if  Rz==0
     Rz=Rz+ 0.00001;
 end
 
 
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
x_int(1)= -(Rx*(Ry*y0 + (R_cyl^2*Rx^2 + R_cyl^2*Ry^2 - Rx^2*y0^2)^(1/2)))/(Rx^2 + Ry^2);
x_int(2)= -(Rx*(Ry*y0 - (R_cyl^2*Rx^2 + R_cyl^2*Ry^2 - Rx^2*y0^2)^(1/2)))/(Rx^2 + Ry^2);

if imag(x_int)
x_int=[];
end
if isnan(x_int)
x_int=[];
end


%if we have correct real results in 2D x-y plane, we continue the
%computations, otherwise we stop, because there is no intersection.

if x_int
    
    y_int=(Ry/Rx)*(x_int)+ y0;
    z_int=((x_int-Lx)/Rx)*Rz + Lz;
    % z_int=((y_int-Ly)/Ry)*Rz + Lz;
    
    % now check whether z_int in the region of z_cylinder or not?
    % if z_cylinder cover the z_int space is fine, otherwise
    
    
    if z_int<0 | z_int > z_cyl_end
%         warning('There is no intersection between line and surface!')
        x_int=[];
        y_int=[];
        z_int=[];
        
    elseif z_int>0 & z_int < z_cyl_end
        
%         disp('The intersection point has been calculated with cylinder!')
    else
        
        %compute the intersection from the line
        %decide from which ends of the cylinder the ray leave the body
        
        x_int_new=([0 z_cyl_end] -Lz)*(Rx/Rz)+Lx;
        y_int_new=([0 z_cyl_end] -Lz)*(Ry/Rz)+Ly;
        
%         disp(z_int)
          Temp=x_int_new.^2+y_int_new.^2; % distance to the center
        
        if z_int(1) >0 & z_int(1) < z_cyl_end
         
            
            
%             x_int (2)=   x_int_new(abs(x_int_new)<R_cyl);
%             y_int (2)= y_int_new(abs(y_int_new)<R_cyl);
              x_int (2)=   x_int_new(Temp<R_cyl^2);
              y_int (2)=   y_int_new(Temp<R_cyl^2);
            z_int(2)=((Rz/Rx)*(x_int (2)-Lx))+Lz;
            
        elseif z_int(2) >0 & z_int(2) < z_cyl_end
%           x_int (1)=   x_int_new(abs(x_int_new)<R_cyl);
%             y_int (1)= y_int_new(abs(y_int_new)<R_cyl);  % should be
%             modified in other functions too
             
             x_int (1)=   x_int_new(Temp<R_cyl^2);
              y_int (1)=   y_int_new(Temp<R_cyl^2);
            
                z_int(1)=((Rz/Rx)*(x_int (1)-Lx))+Lz;
        end
        
        
    end
else
%     warning('There is no intersection between line and surface!')
  
end


    if ~(isempty(z_int) )
% [temp,index_x]=min(abs(x_int-Lx));
% [temp,index_y]=min(abs(y_int-Ly));
% [temp,index_z]=min(abs(z_int-Lz));



norm1=norm([(x_int(1)-Lx),(y_int(1)-Ly),(z_int(1)-Lz)]);
norm2=norm([(x_int(2)-Lx),(y_int(2)-Ly),(z_int(2)-Lz)]);

if norm1 <norm2
    index=1;
else
    index=2;
end


x_int=x_int(index) ;
y_int=y_int(index) ;
z_int=z_int(index);

% x_int=x_int(index_x) ;
% y_int=y_int(index_y) ;
% z_int=z_int(index_z);


    end
    
%     x_int
%     y_int
%     z_int
    

        
    if ~isempty(x_int)
        
        Dir_P2_to_P1=sign([x_int-Lx   y_int-Ly   z_int-Lz ]);
        
        if Dir_P2_to_P1==sign([Rx,Ry,Rz])
            
            xyz_int=[x_int ;y_int ;z_int ];
            
             else
        xyz_int=[];
        end
        
    else
        xyz_int=[];
        
        
    end

%%
%plot 3D line and intersection points
%  plot3D_cylinder_line(R_cyl,z_cyl_end,Rx,Lx,Ry,Ly,Rz,Lz,x_int,y_int,z_int)
%  plot3D_line_interP(R_cyl,z_cyl_end,Rx,Lx,Ry,Ly,Rz,Lz,x_int,y_int,z_int)
%  
 
 
 

%find the normal vectors from the surface
%%
%gradient of the circle in x-y plane
% g=x i+y j 
%to ensure that the inteersection point is on the circle


