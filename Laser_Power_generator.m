
function Power_Actual=Laser_Power_generator (Ax,Ay,nx,ny,ID)

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
        
        % %         if x>0 && y >0
        fun1=@(x,y) exp(-(( (x-x0).*Gauss_Par_X).^2+(  (y-y0).*Gauss_Par_Y).^2)) ;
        
        fun2=@(x,y) exp(-(((x-x0).*Gauss_Par_X2).^2+((y-y0).*Gauss_Par_Y).^2));
        %             amp=area/integral2(fun,-Ax,Ax,-Ay,Ay);
        
        fun3=@(x,y) exp(-(((x-x0)*Gauss_Par_X2).^2+((y-y0).*Gauss_Par_Y2).^2));
        %             amp=area/integral2(fun,-Ax,Ax,-Ay,Ay);
        
        fun4=@(x,y) exp(-(((x-x0)*Gauss_Par_X).^2+((y-y0).*Gauss_Par_Y2).^2));
        %             amp=area/integral2(fun,-Ax,Ax,-Ay,Ay);
        
        
end


%%
%  theta_ray=acosd(v2(3)/norm(v2));
% phi_ray=atan2d(v2(2),v2(1));
%
% delta_phi=+0;%*sign(phi_ray); %deg minus for this value
% delta_theta=+0;%*sign(theta_ray); %deg + positive for this value

%%


xp_laser=linspace(-Ax,Ax,nx);
yp_laser=linspace(-Ay,Ay,ny);

[X,~]=meshgrid(xp_laser,yp_laser);  % ny * nx = size
% Z=zeros(size(X));



Power_Actual=zeros(size(X));

% laser_point_Actual=zeros(1,3);
% laser_point_Actual_all=zeros(nx*ny,7);

if ID(1) ==3
    
    
    amp1=area/integral2(fun1,x0,Ax,y0,Ay);
    amp2=area/integral2(fun2,-Ax,x0,y0,Ay);
    amp3=area/integral2(fun3,-Ax,x0,-Ay,y0);
    amp4=area/integral2(fun4,x0,Ax,-Ay,y0);
    
    amp=(amp1+amp2+amp3+amp4)/16;
    
    
    
    for ii=1:length(xp_laser)
        for jj=1:length(yp_laser)
%             counter=counter+1;
            
            
            
            if xp_laser(ii)>=x0 && yp_laser(jj) >=y0
                
                
                Power_Actual(jj,ii)=amp* fun1(xp_laser(ii),yp_laser(jj));
                
                
            elseif xp_laser(ii)<=x0 && yp_laser(jj) >=y0
                Power_Actual(jj,ii)=amp* fun2(xp_laser(ii),yp_laser(jj));
                
            elseif xp_laser(ii)<=x0 && yp_laser(jj) <=y0
                
                Power_Actual(jj,ii)=amp* fun3(xp_laser(ii),yp_laser(jj));
            elseif xp_laser(ii)>=x0 && yp_laser(jj) <=y0
                Power_Actual(jj,ii)=amp* fun4(xp_laser(ii),yp_laser(jj));
            end
            
            %         Power_Actual(jj,ii)= laser_point_Actual_all(counter,4);
        end
    end
    
    
    
    
else
    
    
    
    for ii=1:length(xp_laser)
        for jj=1:length(yp_laser)
            
            
            %energy will be specified here
            Power_Actual(jj,ii)=amp*fun(xp_laser(ii),yp_laser(jj));
            
            
            
        end
    end
    
end

Power_Actual=(Power_Actual/ sum(Power_Actual(:)))*nx*ny;

