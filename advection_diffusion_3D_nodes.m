
function advection_diffusion_3D_nodes(Lx,Ly,Lz,eleNx,eleNy,eleNz)

% clc;
% clear all;
% Lx=20e-3;  % 20 mm
% Ly=10e-3;  % 10 mm
% Lz=10e-3;  % 5 mm 


% eleNx=30-1;
% eleNy=15-1;

% eleNx=2+1;  % number of element eleNx-1
% eleNy=1+1;  % number of element eleNy-1
% eleNz=2+1;  % number of element eleNz-1

% dx=Lx/eleNx;  
% dy=Ly/eleNy;
% xnode=0:dx:Lx;
% ynode=0:dy:Ly;

counter=0;


xnode=linspace(0,Lx,eleNx);
ynode=linspace(0,Ly,eleNy);
znode=linspace(0,Lz,eleNz);


fileID1 = fopen('nodes.txt','w');
%  fileID2 = fopen('elements.txt','w');
 
 
 n=(eleNx)*(eleNy)*eleNz;
 
% fileID5 = fopen('Mesh_3D.plt','w');
% 
% fprintf(fileID5,' "TITLE = 3D-Temperature  by Amin zaami" \n');
% fprintf(fileID5,'VARIABLES = "X", "Y", "Z", "Temp" \n');
% fprintf(fileID5,'zone  N=    %d E=   %d DATAPACKING = POINT, ZONETYPE = FEBRICK  \n',n, (eleNx-1)*(eleNy-1)*(eleNz-1) );
% 

 
 
 

 for kk=1:eleNz
for ii=1: eleNx
    
    for jj=1: eleNy
   counter=counter+1;
    fprintf(fileID1,'%d, %f, %f  %f \r\n',counter, xnode(ii), ynode(jj), znode(kk));     
        
%     fprintf(fileID5,' %f %f  %f %f \r\n', xnode(ii), ynode(jj), znode(kk), rand(1,1));   
    end
    
end
 end

%% element creation
% % 
% shift=(eleNx)*(eleNy);
% num_ele_plane=(eleNx-1)*(eleNy-1);  %number of element in a plane
% 
% counter_e=0;
% n_th=0;
% for ii=1: (eleNx-1)*(eleNy-1)*(eleNz-1)    % number of element
%     
%    
%   
%   
%    if rem(ii,eleNx-1)==1
%      counter_e=counter_e+1;
%    end
%    
%    
%    
%    node_plane=[counter_e, counter_e+1,eleNx+counter_e+1,eleNx+counter_e];
%     node_plane_up=shift+[counter_e, counter_e+1,eleNx+counter_e+1,eleNx+counter_e];
%    
% %    fprintf(fileID2,'%d, %d, %d, %d, %d, %d, %d, %d, %d \r\n',ii, [node_plane node_plane_up]); 
%       fprintf(fileID5,' %d %d %d %d %d %d %d %d \r\n', [node_plane node_plane_up]); 
%    
%    if mod(ii,num_ele_plane)==0
%       n_th=n_th+1;
%      counter_e=n_th*shift-1;
%       
%    end
%     
%     counter_e=counter_e+1;
%   
% end
% 
% 








   fclose(fileID1);
%      fclose(fileID2);
%        fclose(fileID5);