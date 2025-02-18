
    
    
function  Tecplot_Matlab_creator_2D (tt,N,xnode,ynode,v_original,Lx,Ly,T,fileID6)
    
    
              
    
        %      fileID6 = fopen(sprintf('Temp_3D-FDM_V_0to%f.plt',v_original),'a');
    
    
    
    
    
            fprintf(fileID6,' "TITLE = 2D-Temperature  by Amin zaami" \n');
            fprintf(fileID6,'VARIABLES = "X", "Y", "Temp" \n');
            fprintf(fileID6,'zone  N=    %d E=   %d F=FEPOINT  ET=QUADRILATERAL  \n',N, (xnode-1)*(ynode-1) );
            fprintf(fileID6,'STRANDID=1, SOLUTIONTIME=  %d \n',tt );
    
            eleNx=xnode;
            eleNy=ynode;
%             eleNz=znode;
    
    
            xnode_P=linspace(0,Lx,eleNx);
            ynode_P=linspace(0,Ly,eleNy);
%             znode_P=linspace(0,Lz,eleNz);
            counter=0;
    
    
%             for kk=1:znode
                for ii=1: xnode
    
                    for jj=1: ynode
                        counter=counter+1;
    
                        fprintf(fileID6,' %f %f  %f \r\n', xnode_P(ii), ynode_P(jj), T(counter));
                    end
    
                end
%             end
    
    
    
    % element creation
    
%             shift=(eleNx)*(eleNy);
%             num_ele_plane=(eleNx-1)*(eleNy-1);  %number of element in a plane
    
%             node_plane=zeros(4,1);
%             node_plane_up=zeros(4,1);
    
    
            counter_e=0;
%             n_th=0;
            for ii=1: (eleNx-1)*(eleNy-1)    % number of element
    
                if rem(ii,eleNy-1)==1
                    counter_e=counter_e+1;
                end
    
                node_plane=[counter_e, counter_e+1,eleNy+counter_e+1,eleNy+counter_e];
%                 node_plane_up=shift+[counter_e, counter_e+1,eleNy+counter_e+1,eleNy+counter_e];
    
                fprintf(fileID6,' %d %d %d %d \r\n',node_plane );
    
%                 if mod(ii,num_ele_plane)==0
%                     n_th=n_th+1;
%                     counter_e=n_th*shift-1;
%     
%                 end
    
                counter_e=counter_e+1;
    
            end
    
    
    
    
    
%         fclose(fileID6);