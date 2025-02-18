
function Tot_fgh_prim=F_ellipsoid_path(a,b,th,...
    c,path_angle,theta0,phi0)


%%
% f=a*cos(phi0+tan(path_angle)*th)*sin(th+theta0);
% g=b*cos(phi0+tan(path_angle)*th)*cos(th+theta0);
% h=c*sin(phi0+tan(path_angle)*th);

% f_prim=diff(f,th);
% g_prim=diff(g,th);
% h_prim=diff(h,th);

f_prim = a.*cos(th + theta0).*cos(phi0 + th.*tan(path_angle)) - a.*sin(th + theta0).*tan(path_angle).*sin(phi0 + th.*tan(path_angle));
g_prim =- b.*sin(th + theta0).*cos(phi0 + th.*tan(path_angle)) - b.*cos(th + theta0).*tan(path_angle).*sin(phi0 + th.*tan(path_angle));
h_prim = c.*tan(path_angle).*cos(phi0 + th.*tan(path_angle));


Tot_fgh_prim=sqrt(f_prim.^2+g_prim.^2+h_prim.^2);


%        q = integral(Tot_fgh_prim,0,th);
%
%           F=q-L_sub;

% end
