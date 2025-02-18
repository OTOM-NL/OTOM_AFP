

function F=F_path(x0,y0,z0,a,c2,m,b,z0_E2,delta_s0,rot,th)



  phi=atan((a/(c2*m))*(sin(th)));%+asin(1*z0/c1);
        

 xp=b.*cos(phi).*cos(th+rot);
    yp=a.*cos(phi).*sin(th+rot);       
        zp=c2.*sin(phi) +z0_E2 ; %.*ones(size(theta));

       F=norm([xp-x0,yp-y0,zp-z0])-delta_s0;
        
        

      
    