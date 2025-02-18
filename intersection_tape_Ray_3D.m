function [intersection_points,S1,S2]=intersection_tape_Ray_3D(laser_source_P,laser_direction,R,W,deg_T,n1,n2,n3,x0,y0,z0,L_flat,W_R,theta_ind)


% W_R=1.5*W;   % we suppose that Roller width is twith than Width of tape


% Ray
x_int=[];
y_int=[];
z_int=[];


Lx=laser_source_P(1);
Ly=laser_source_P(2);
Lz=laser_source_P(3);

 Rx=laser_direction(1);
 Ry=laser_direction(2);
 Rz= laser_direction(3);



% Lx= -20; % input('x position of Laser = ');    %corrdinate of laser
% Ly=-20;% input('y position of Laser =');
% Lz=5;%  input('z position of Laser =');
% 
% Rx=1.1; %input('direction of ray in x-direction =');
% Ry= 1.99;% input('direction of ray in y-direction =');
% Rz= 0.01; %input('direction of ray in z-direction =');

% normalize the Rx,Ry, Rz
% normR=norm([Rx Ry Rz]);
% Rx=Rx/normR;
% Ry=Ry/normR;
% Rz=Rz/normR;

%line equation 
%y=mx+y0  , which m=Ry/Rx
m=Ry/Rx;
y0_L=Ly -(Ry/Rx)*Lx;

% syms x    % for Matlab 2013a
% f=@(x) x^2 -R^2 + ((Ry/Rx)*x + y0_L -R )^2;  % for matlab 2015 or later
% f= x^2 -R^2 + ((Ry/Rx)*x + y0_L -R )^2;     % for Matlab 2013a
% x_int_trial=double(solve(f,'Real', true));
% x_int=x_int_trial(x_int_trial<0);


% ans =
%  
x_int_trial=zeros(2,1);
x_int_trial(1)= (Rx*(R*Ry - Ry*y0_L + (R^2*Ry^2 + 2*R*Rx^2*y0_L - Rx^2*y0_L^2)^(1/2)))/(Rx^2 + Ry^2);
x_int_trial(2)= -(Rx*(Ry*y0_L - R*Ry + (R^2*Ry^2 + 2*R*Rx^2*y0_L - Rx^2*y0_L^2)^(1/2)))/(Rx^2 + Ry^2);
 
if imag(x_int_trial)
x_int_trial=[];
end

[temp,index_x]=min(abs(x_int_trial-Lx));
x_int=x_int_trial(index_x) ;




S1=false(1);
S2=false(1);



deg=(deg_T);

if deg_T< theta_ind
    deg=(theta_ind);

end





if   x_int<=0 &  x_int>=-R*sind(deg)
    y_int=(Ry/Rx)*(x_int)+ y0_L;   %use the line equation
    
    if y_int>=0 &  y_int<=R*(1-cosd(deg))
        z_int=((x_int-Lx)/Rx)*Rz + Lz;
         if z_int>=-W/2 & z_int <=W/2
%         if z_int>=-W_R/4 & z_int <=W-W_R/4
%             disp('intersection between line and tape in curved rigion')
%     if z_int>=W_R/4 & z_int <=W+W_R/4
            S1=true(1);
        else
            z_int=[];
            x_int=[];
            y_int=[];
        end
    else
        x_int=[];
        y_int=[];
        
    end
else
    x_int=[];
    % go for calculation for plane intersection
end
%%
%calculation for intersection with part plane- flat
if isempty(x_int)
    x_int=((n1/n2)*(x0)-y0_L+y0)/((n1/n2)+m);
    
    if   x_int<=-R*sind(deg) &  x_int>=-R*sind(deg)-L_flat*cosd(deg_T)
        y_int=(Ry/Rx)*x_int+y0_L;
        
        if y_int>=R*(1-cosd(deg)) &  y_int<=R*(1-cosd(deg))+L_flat*sind(deg_T)
            z_int=(Rz/Rx)*(x_int-Lx)+Lz;
           if  z_int>=-W/2 & z_int <=W/2
%             if z_int>=-W_R/4 & z_int <=W-W_R/4   % the relation between tape and roller width
%                 disp('intersection between line and tape in curved rigion');
%     if  z_int>=W_R/4 & z_int <=W+W_R/4
                S2=true(1);
                
                
            else
                z_int=[];
                x_int=[];
                y_int=[];
            end
        else
            x_int=[];
            y_int=[];
            
        end
    else
        x_int=[];
    end
    
end

%%
% if x_int
% %     disp('There is a intersection with tape ! ')
% else
% %     disp('There is no intersection with tape !')
% end
%  intersection_points=([x_int;y_int;z_int]);
 
 
 
 if ~isempty(x_int)
        
        Dir_P2_to_P1=sign([x_int-Lx   y_int-Ly   z_int-Lz ]);
        
        if Dir_P2_to_P1==sign([Rx,Ry,Rz])
            
            intersection_points=[x_int ;y_int ;z_int ];
            
             else
        intersection_points=[];
        end
        
    else
        intersection_points=[];
        
        
    end
 

