
% function to calculate neighboring nodes and corresponding weight function, and some additional data for Thermal analysis 

function Pre_Thermal_UOT(W_tape,thick_T,R_tape,L_flat,deg_tape,...
    R_cyl, W_R,node_space,Angular_space,L_flat_space,...
    H_indentation,nip_point_M_all,Rot_Roller_axis_all,...
    UOT_pathfile,CV_mesh,text_status)



mkdir(UOT_pathfile);

steps=length(nip_point_M_all);

% n_dis_sub_all=zeros(steps,ii,1:4)
Points_in_domain_all_sub=cell(steps,1);
Points_in_domain_all_Tape=cell(steps,1);
 Tape_points_Data=cell(steps,4);

for ss=1:steps
    
   set(text_status,'String',[num2str((ss/steps)*100) '%']);
    drawnow;
    % should be modified
    
%     fileID21 = fopen(strcat(UOT_pathfile, sprintf('Ponits_in_domain%d.txt',ss)),'w');
%     
%     fprintf(fileID21,' X         ,Y         ,Z       \r\n');
    
    
    
    fileID1 = fopen(strcat(UOT_pathfile, sprintf('Cylinder_ints%d.txt',ss)),'r');
    fileID2 = fopen(strcat(UOT_pathfile, sprintf('Tape_ints%d.txt',ss)),'r');
    
    
    % >> Pause to let the computer write data into txt file
    %     pause(0.01); %% was 0.05
    
    % for the substrate
    
    % S_ind >> source index
    % get XYZ coordiante 3, Energy 1, and incoming ray index 1
    node = textscan(fileID1,' %f %f %f %f %f','Delimiter',',','HeaderLines',1) ;
    Nodes_optical=cell2mat(node);
    
    %      S_ind_Nodes_optical=cell2mat(node(5));
    
    %% for the tape
       % get XYZ coordiante 3, Energy 1, and incoming ray index 1
    node_T = textscan(fileID2,' %f %f %f %f %f','Delimiter',',','HeaderLines',1) ;
    Nodes_optical_T=cell2mat(node_T);
    
    
    
    theta_ind=acosd((R_tape-H_indentation)/R_tape);  %theta_ind in degree
    
    %     Nip_Mov =R_tape*sind(theta_ind);
    
    
    %% Generating thermal points and decide which optical point is in the domain
    %we already know that the domain for the tape, do not need to calculate
    %Now the program should decide which part is exposed for laying substrate which depends on Nip-point
    
    
    mat_size=size(CV_mesh);
    
    X=reshape(CV_mesh(ss,1,:,:),mat_size(3:4));
    Y= reshape(CV_mesh(ss,2,:,:),mat_size(3:4));
    Z=reshape(CV_mesh(ss,3,:,:),mat_size(3:4));
    

    
    % make boarder points > CCW > to know whether a point is inside or not!
    Ps1=[X(1:end,1),Y(1:end,1),Z(1:end,1) ];
    Ps2=[X(end,2:end)' ,Y(end,2:end)' ,Z(end,2:end)' ];
    Ps3=[X(end-1:-1:1,end),Y(end-1:-1:1,end),Z(end-1:-1:1,end)];
    Ps4=[X(1,end-1:-1:2)',Y(1,end-1:-1:2)',Z(1,end-1:-1:2)' ];
    
    Boarders=[Ps1;Ps2;Ps3;Ps4]';
    
    
    
    %    plot3(Boarders(1,:),Boarders(2,:),Boarders(3,:),'c--')
    %    axis equal;
    %    hold on;
    
    % Transpose because numbering in Thermal problem is like this
    X=X';
    Y=Y';
    Z=Z';
    % all Eulerian points
    points=[X(:),Y(:),Z(:)];
    
    
    
    
    
%     
%     if R_cyl(1) ~= 0
%         
%         UpperBoarder_index=Boarders(2,:)>0;
%         BottomBoarder_index=Boarders(2,:)<0;
%         
% %         UpperBoarder= Boarders(:,UpperBoarder_index);
% %         BottomBoarder= Boarders(:,BottomBoarder_index);
%         
%         
%     end
    
    
    [Len_Nodes_optical,~]=size(Nodes_optical);
%     Points_in_domain=zeros(Len_Nodes_optical,13);
     Points_in_domain=zeros(Len_Nodes_optical,9);
     
    counter=0;
    
    dist_opt_sub=length(points);
    
    
    
   [ynode,xnode]=size(X);
%%

% No_dev_L_Cut= ;%length(E_points)/node_num_Z;
% xnode=node_num_Z;
% ynode=No_dev_L_Cut;


% X=zeros(xnode,ynode);
% Y=zeros(xnode,ynode);
% Z=zeros(xnode,ynode);
% Intensity=zeros(xnode,ynode);

% for kk=1:xnode
%     index=((1:ynode)*xnode)-kk+1;
% %     X(kk,1:end)=points(index,1)';
% %     Y(kk,1:end)=points(index,2)';
% %     Z(kk,1:end)=points(index,3)';
% %     Intensity(kk,1:end)=points(index,4)'; % RHS in the thermal model
% end



S(:,:,1) = X;
S(:,:,2) = Y;
S(:,:,3) = Z;

T = Nodes_optical(:,1:3); %[Points.XData',Points.YData',Points.ZData'];

d = 1.1*norm(squeeze(S(1,1,:))'-squeeze(S(2,2,:))');

% tic
inSurf =zeros(1,length(T(:,1)));
[s1,s2,~]=size(S);    
Gi = 1;
Gj = 1;
%Main code
for i=1:(s1-Gi)
    for j = 1:(s2-Gj)
        p1 =squeeze(S(i,j,:))';
        p2 =squeeze(S(i,j+Gj,:))';
        p3 =squeeze(S(i+Gi,j,:))';
        p4 =squeeze(S(i+Gi,j+Gj,:))';
        for z=1:size(T,1);
            if inSurf(z) == 0 && norm(T(z,:)-(0.5*(p1+p3)))< d;
                if SideCheck(p1,p2,p3,T(z,:))==1 || SideCheck(p4,p3,p2,T(z,:))==1; 
                    inSurf(z)= z;
                end;
            end
        end;
    end;
end;

inSurf(inSurf==0)=[]; 
    
    
 for ii=1:length(inSurf)
    
    ind=inSurf(ii);   
        
      
        %     plot3(UpperBoarder([1],:),UpperBoarder([2],:),UpperBoarder([3],:),'r*')
        
        
%         if (THETA_all >=359 & THETA_all<=365)
            
            
            counter=counter+1;
            Points_in_domain(counter,1:5)=Nodes_optical(ind,1:5);
            %                      plot3(Nodes_optical(ind,1),Nodes_optical(ind,2),Nodes_optical(ind,3),'cp');
            
            
%             fprintf(fileID21,' %12.8f %12.8f %12.8f \r\n',Nodes_optical(ind,1),Nodes_optical(ind,2),Nodes_optical(ind,3));
%             
            
            
            % find four neighbouring points
            for dd=1:length(points)
                dist_opt_sub(dd)=norm(abs(Points_in_domain(counter,1:3)- points(dd,1:3)  ));
                
            end
            
%             n_dis >>  is a vector  > data are saved for each point
            [~,n]=sort(dist_opt_sub);
            % find nearest Z in specified theta angle
           
%             n_dis_sub_all(ss,ind,1:4)=n_dis(1:4);

            P_temp(1,1:3)=points(n(1),1:3)  ;
            P_temp(2,1:3)= points(n(2),1:3) ;
%             P_temp(3,1:3)= points(n(3),1:3) ;
%             P_temp(4,1:3)=points(n(4),1:3) ;
            
            
            
            w1_temp=1/norm([Nodes_optical(ind,1)-P_temp(1,1),...
                Nodes_optical(ind,2)-P_temp(1,2),...
                Nodes_optical(ind,3)-P_temp(1,3)]);
            
            w2_temp=1/norm([Nodes_optical(ind,1)-P_temp(2,1),...
                Nodes_optical(ind,2)-P_temp(2,2),...
                Nodes_optical(ind,3)-P_temp(2,3)]);
            
%             w3_temp=1/norm([Nodes_optical(ind,1)-P_temp(3,1),...
%                 Nodes_optical(ind,2)-P_temp(3,2),...
%                 Nodes_optical(ind,3)-P_temp(3,3)]);
%             
%             w4_temp=1/norm([Nodes_optical(ind,1)-P_temp(4,1),...
%                 Nodes_optical(ind,2)-P_temp(4,2),...
%                 Nodes_optical(ind,3)-P_temp(4,3)]);
            
            
            cte=w1_temp+w2_temp;%+w3_temp+w4_temp;
            
            % normalized weighting values for each point based on the distance from the
            % optical point
            w1=w1_temp/cte;
            w2=w2_temp/cte;
%             w3=w3_temp/cte;
%             w4=w4_temp/cte;

            
%             Points_in_domain(counter,6:13)=[n(1:4),w1 w2 w3 w4] ;
 Points_in_domain(counter,6:9)=[n(1:2),w1 w2 ] ;
            % Save weight variables w1... w4 and neighbouring point index
            % n1...n4
%                    fprintf(fileID21,' %12.8f %12.8f %12.8f %d %d %d %d  \r\n',Nodes_optical(ind,1),Nodes_optical(ind,2),Nodes_optical(ind,3),n(1:4));
            


            
        
        
        
    end
    
    
    % To prevent zero elements of the matrix
    
    if counter <Len_Nodes_optical
        Points_in_domain(counter+1:end,:)=[];
        
    end
    
 Points_in_domain_all_sub{ss}= Points_in_domain  ;
    
    
    
    %% Now calculation for Tape
    %define the thermal points for the tape
    % all the optical points of the tape are in the doimain
    
    R_tape_th=R_tape+thick_T;
    
    
    % th_v indicate which degree is the tangent for the Tape part
    %Get from Optical_UOT.mat
    
    Rot_Roller_axis= reshape(Rot_Roller_axis_all(ss,:,:),[3 3]);
    Roller_Pos_TV=nip_point_M_all(ss,:,:);
    
    
    [points_T_G,z,dl,th_v]=Tape_points (deg_tape,R_tape_th,W_tape,L_flat,Roller_Pos_TV,Rot_Roller_axis,W_R,...
        node_space,Angular_space,L_flat_space,theta_ind);
    
    % cell matrix to fill
    Tape_points_Data{ss,4}=th_v;
     Tape_points_Data{ss,1}=points_T_G;
      Tape_points_Data{ss,2}=z;
       Tape_points_Data{ss,3}=dl;
    
    
    
    points_T_G=points_T_G';
    Temp_norm=zeros(length(points_T_G),1);
    
%       counter=0;
    
    [row_T,col]=size(Nodes_optical_T);
    
%     Points_in_domain_T=zeros(row_T,13);
       Points_in_domain_T=zeros(row_T,9);
       
    for ii=1:row_T
        
        % limit the scope search, search in these limits
        
        %     points_T_G(points_T_G(:,1)>)
        %  X_limit= X_op-delta_X  <X_op<X_op+delta_X;
        %  Y_limit=   Y_op-delta_Y  <Y_op<Y_op+delta_Y;
        % Z_limit=     Z_op-delta_Z  <Z_op<Z_op+delta_Z;
        
        %X_limit,Y_limit, and Z_limit should have a same size
        
        % XYZ_limits=[X_limit,Y_limit,Z_limit]
        
%      counter=counter+1;
            Points_in_domain_T(ii,1:5)=Nodes_optical_T(ii,1:5);
            
            
        
        for jj=1:length(points_T_G)
            Temp_norm(jj)=norm(points_T_G(jj,1:3)-Nodes_optical_T(ii,1:3));
        end
        
        [Temp,n_T]=sort(Temp_norm');
        % store n_T(1),n_T(2),n_T(3),n_T(4)
        
        
        
         
        P_temp(1,1:3)=points_T_G(n_T(1),1:3);
        P_temp(2,1:3)= points_T_G(n_T(2),1:3);
%         P_temp(3,1:3)= points_T_G(n_T(3),1:3);
%         P_temp(4,1:3)=points_T_G(n_T(4),1:3);
        
        
        %%
        
        w1_temp=1/norm([Nodes_optical_T(ii,1)-P_temp(1,1),...
            Nodes_optical_T(ii,2)-P_temp(1,2),...
            Nodes_optical_T(ii,3)-P_temp(1,3)]);
        
        w2_temp=1/norm([Nodes_optical_T(ii,1)-P_temp(2,1),...
            Nodes_optical_T(ii,2)-P_temp(2,2),...
            Nodes_optical_T(ii,3)-P_temp(2,3)]);
        
%         w3_temp=1/norm([Nodes_optical_T(ii,1)-P_temp(3,1),...
%             Nodes_optical_T(ii,2)-P_temp(3,2),...
%             Nodes_optical_T(ii,3)-P_temp(3,3)]);
%         
%         w4_temp=1/norm([Nodes_optical_T(ii,1)-P_temp(4,1),...
%             Nodes_optical_T(ii,2)-P_temp(4,2),...
%             Nodes_optical_T(ii,3)-P_temp(4,3)]);
        
        
        cte=w1_temp+w2_temp;%+w3_temp+w4_temp;
        
        % normalized weighting values for each point based on the distance from the
        % optical point
        w1_T=w1_temp/cte;
        w2_T=w2_temp/cte;
%         w3_T=w3_temp/cte;
%         w4_T=w4_temp/cte;
        
         % Save weight variables w1... w4 and neighbouring point index
            % n1...n4
            % since all points from tape are inside
%              Points_in_domain_T(ii,6:13)=[n_T(1:4),w1_T w2_T w3_T w4_T] ;
         Points_in_domain_T(ii,6:9)=[n_T(1:2),w1_T w2_T ] ;
    end
    
     Points_in_domain_all_Tape{ss}= Points_in_domain_T  ;
    
    
    
    fclose(fileID1);
    fclose(fileID2);
    
end





%% Data to Save

  dir2save=UOT_pathfile;
     save(strcat(dir2save,'\Pre_Thermal_UOT.mat'),'Points_in_domain_all_sub',...
         'Points_in_domain_all_Tape','Tape_points_Data');

     
% Ly=2*w;
% Lx=L_prim-Nip_Mov;
% Lz=thick_sub;
% 
% mat_size=size(CV_mesh);
% ynode=mat_size(4);
% xnode=mat_size(3);
% 
% znode=Node_z_3D_thermal;
% 
% if znode==0
%     znode=4;
% end
%      
%      advection_diffusion_3D_nodes(Lx,Ly,Lz,xnode,ynode,znode);

