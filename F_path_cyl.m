

function F=F_path_cyl(L,tv,delta_s0,sui,R,No_dev_L)



if length(R)==3
   
    c1=R(2);
     c2=R(3);
else
       c1=R(1);
     c2=R(1);
%       R_cyl=R_cyl(1);
    
end

a=R(1);
b=R(1);



x0=tv(1);
y0=tv(2);
z0=tv(3);


phi=sui*(pi/180); %pi/2-alpha; 

%%
% for cylinder

th=L*cos(phi)/a; %in radian
H=L*sin(phi);

 theta0=atan2(x0,y0);

% No_dev_L=40;
t_L=linspace(theta0,theta0+th,No_dev_L);
x_L=0+a*sin(t_L);
y_L=0+a*cos(t_L);
z_L=z0+linspace(0,H,No_dev_L);





%   phi=atan((a/(c2*m))*(sin(th)));%+asin(1*z0/c1);
%         
% 
%  xp=b.*cos(phi).*cos(th+rot);
%     yp=a.*cos(phi).*sin(th+rot);       
%         zp=c2.*sin(phi) +z0_E2 ; %.*ones(size(theta));

       F=norm([x_L(2)-x0,y_L(2)-y0,z_L(2)-z0])-delta_s0;
        
        

      
    