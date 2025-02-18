
% syms a b  c lam

function [th_limit,counter]=ellipsoid_path_integration(a,b,c,path_angle,L_sub,theta0,phi0)

% syms th;

% a=5; along x-axis
% b=2; along y-axis
% c=7; along z-axis
% path_angle=0;
% path_angle=pi/2;
% path_angle=path_angle+ (.001);  % to avoid numerical instability

% f=a*cos(phi0+tan(path_angle)*th)*sin(th+theta0);
% g=b*cos(phi0+tan(path_angle)*th)*cos(th+theta0);
% h=c*sin(phi0+tan(path_angle)*th);

% f_prim=diff(f,th);
% g_prim=diff(g,th);
% h_prim=diff(h,th);

% f_prim = a*cos(path_angle*th)*cos(th) - a*path_angle*sin(path_angle*th)*sin(th)
% g_prim =- b*cos(path_angle*th)*sin(th) - b*path_angle*sin(path_angle*th)*cos(th)
%   h_prim = c*path_angle*cos(path_angle*th)

% f_prim = a*cos(th + theta0)*cos(phi0 + th*tan(path_angle)) - a*sin(th + theta0)*tan(path_angle)*sin(phi0 + th*tan(path_angle));
% g_prim =- b*sin(th + theta0)*cos(phi0 + th*tan(path_angle)) - b*cos(th + theta0)*tan(path_angle)*sin(phi0 + th*tan(path_angle));
%   h_prim = c*tan(path_angle)*cos(phi0 + th*tan(path_angle));


% Tot_fgh_prim=f_prim^2+g_prim^2+h_prim^2;


% diff(sqrt(Tot_fgh_prim),th)
% initial th can be pi/2-path_angle-pi/200

% counter=0;
% counter2=0; 
% Length1=L_sub*2;
% Length2=L_sub*2;
% 
% th_limit1=pi/2-path_angle; % to avoid Nan in the results
% th_limit2=pi/2-path_angle; % to avoid Nan in the results
% 
% Tol=1e-3;
% 
% while abs(Length1-L_sub) >Tol
%     th_limit1=th_limit1*(1-((Length1-L_sub)/Length1 ));
% 
% Length1 = double(vpaintegral(sqrt(Tot_fgh_prim),th,[ 0 th_limit1]));
%  
% counter=counter+1;
% end

% th_limit1

%%
% based on the numerical approach not symbolic

counter=0;
counter2=0; 
Length1=L_sub*2;
Length2=L_sub*2;

fun = @(xx) F_ellipsoid_path(a,b,xx,c,path_angle,theta0,phi0);    % function of x alone

th_limit1=pi/2-path_angle; % to avoid Nan in the results
th_limit2=pi/2-path_angle; % to avoid Nan in the results

Tol=1e-3;

while abs(Length1-L_sub) >Tol
    th_limit1=th_limit1*(1-((Length1-L_sub)/Length1 ));
% Length = double(vpaintegral(sqrt(Tot_fgh_prim),th,[0 th_limit]));
 % Length1 = double(vpaintegral(sqrt(Tot_fgh_prim),th,[ 0 th_limit1]));
  Length1= integral(fun,0,th_limit1);
counter=counter+1;
end

% end_limit = fzero(fun,[0 pi]);


% th_limit1

%%







while abs(Length2-L_sub) >Tol
    th_limit2=th_limit2*(1-((Length2-L_sub)/Length2 ));

 % Length2 = double(vpaintegral(sqrt(Tot_fgh_prim),th,[ -th_limit2 0]));
 
Length2 =integral(fun,-th_limit2, 0);

counter2=counter2+1;
end

th_limit=[th_limit1 th_limit2];
 