
% clc;
% clear all;
function Node_elem_creeation(Lx,Ly,eleNx,eleNy)

% Lx=20e-3;
% Ly=10e-3;

eleNx=eleNx-1;
eleNy=eleNy-1;

dx=Lx/eleNx;  
dy=Ly/eleNy;

xnode=0:dx:Lx;
ynode=0:dy:Ly;

counter=0;




fileID1 = fopen('nodes.txt','w');
 fileID2 = fopen('elements.txt','w');

for ii=1: eleNy+1
    
    for jj=1: eleNx+1
   counter=counter+1;
    fprintf(fileID1,'%d, %f, %f \r\n',counter, xnode(jj), ynode(ii));     
        
    end
    
end

%% element creation


counter_e=0;
for ii=1: eleNx*eleNy
    
   
  
  
   if rem(ii,eleNx)==1
     counter_e=counter_e+1;
   end
   fprintf(fileID2,'%d, %d, %d, %d, %d \r\n',ii, counter_e, counter_e+1,eleNx+counter_e+2,eleNx+counter_e+1);     
    
    counter_e=counter_e+1;
  
end

   fclose(fileID1);
     fclose(fileID2);
