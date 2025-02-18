

%% PLOT section of ALL: Rays, intersections > for reading the Data

function Optical_post_processing(jobname,R_cyl,z_cyl_end,nip_point_M_all,Rot_Roller_axis_all,...
    N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,H_indentation,tv_def_all,CV_mesh,UOT_pathfile,...
      Laser_head,Total_energy)

ss=1;


%     fileID1 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Cylinder_ints%d.txt'),jobname,ss),'r');   % file includes the xyz + ID
%     fileID2 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Tape_ints%d.txt'),jobname,ss),'r');     % file includes the xyz + ID
%     fileID3 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Roller_ints%d.txt'),jobname,ss),'r');
%
%
%     %%  instead of plotting in each loop, plot all afterwards
%     fileID4 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\int_xyz%d.txt'),jobname,ss),'r');
%     fileID5 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Laser_Rays%d.txt'),jobname,ss),'r');
%     fileID6 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Normal_vectors_Mandrel%d.txt'),jobname,ss),'r');
%     fileID7 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Reflection_vector_Mandrel%d.txt'),jobname,ss),'r');
%     fileID8 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Normal_vectors_Tape%d.txt'),jobname,ss),'r');
%     fileID9 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Reflection_vector_Tape%d.txt'),jobname,ss),'r');
%
%     fileID11 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Ref_int_xyz%d.txt'),jobname,ss),'r');
%




fig_HO=figure('Name','Optical UOT Analysis','NumberTitle','off');
cameratoolbar(fig_HO, 'Show');
javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

set(gcf,'color','w')

% create structure of handles
% Fig_handles = guihandles(fig_HO);
% % Add some additional data as a new field called numberOfErrors
% % Fig_handles.numberOfErrors = 0;
% % Save the structure
% guidata(fig_HO,Fig_handles)




Fig_handles = guidata(gcbo);


% plot Mandrel
% First is cylinder , second is bottom dome part, third is upper Dome part
Mandrel_plot=plot3D_cylinder(R_cyl,z_cyl_end);
assignin('base','Mandrel_plot',Mandrel_plot);


set(fig_HO,'Visible','on');

% plot3(tv(1),tv(2),tv(3),'g*','Markersize',5)


%    Roller_Pos_TV=nip_point_M_all(ss,:,:);


%     tv(1)=Roller_Pos_TV(1)-(-tv_def(1));
%     tv(2)=Roller_Pos_TV(2)-tv_def(2);
%     tv(3)=Roller_Pos_TV(3)+tv_def(3);
%
%      Rot_Roller_axis= reshape(Rot_Roller_axis_all(ss,:,:),[3 3]);
% output is real center of roller without deformation, end of arrow of
% the local axis
%   [n1,n2,n3,x0,y0,z0,Fig_handles.h_R,Fig_handles.h_T,Center_roller]=Tape_plot_3D_Tshow(N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,Rot_Roller_axis,H_indentation,tv_def);

% assignin('base','h_R',Fig_handles.h_R);
% assignin('base','h_T',Fig_handles.h_T);




%     XYZ_int = textscan(fileID4,'%f %f %f','Delimiter',',','HeaderLines',0) ;
%     XYZ_int=cell2mat(XYZ_int);
%     Fig_handles.Laser_int_points=plot3( XYZ_int(:,1), XYZ_int(:,2), XYZ_int(:,3) ,'k.');
%
%
%
%     XYZ_Rays = textscan(fileID5,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     XYZ_Rays=cell2mat(XYZ_Rays);
%    Fig_handles.Graphics_laser= plot3( XYZ_Rays(:,1), XYZ_Rays(:,2), XYZ_Rays(:,3) ,'g:','Linewidth',2 );
%
%
%
%     Normal_vectors_Mandrel = textscan(fileID6,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     Normal_vectors_Mandrel=cell2mat(Normal_vectors_Mandrel);
%     Fig_handles.Normal_lines1=plot3( Normal_vectors_Mandrel(:,1), Normal_vectors_Mandrel(:,2), Normal_vectors_Mandrel(:,3) ,'y:');
%
%
%
%
%     Reflection_vector_Mandrel = textscan(fileID7,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     Reflection_vector_Mandrel=cell2mat(Reflection_vector_Mandrel);
%     Fig_handles.Graphics_reflection1=plot3( Reflection_vector_Mandrel(:,1), Reflection_vector_Mandrel(:,2), Reflection_vector_Mandrel(:,3) ,'r--');
%




%     Normal_vectors_Tape = textscan(fileID8,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     Normal_vectors_Tape=cell2mat(Normal_vectors_Tape);
%    Fig_handles.Normal_lines2= plot3( Normal_vectors_Tape(:,1), Normal_vectors_Tape(:,2), Normal_vectors_Tape(:,3) ,'y:');
%
%
%
%
%     Reflection_vector_Tape = textscan(fileID9,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     Reflection_vector_Tape=cell2mat(Reflection_vector_Tape);
%     Fig_handles.Graphics_reflection2=plot3( Reflection_vector_Tape(:,1), Reflection_vector_Tape(:,2), Reflection_vector_Tape(:,3) ,'b--');
%
%
%
%     Ref_int_xyz = textscan(fileID11,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     Ref_int_xyz=cell2mat(Ref_int_xyz);
%     Fig_handles.reflection_int=plot3( Ref_int_xyz(:,1), Ref_int_xyz(:,2), Ref_int_xyz(:,3) ,'y.');


%% Changing the graphical object, Visibility or deleting

%     Graphics_laser=findobj(h,'color','g','LineStyle',':');
%     assignin('base','Graphics_laser',Fig_handles.Graphics_laser);
%     % set(Graphics_laser,'Visible','off');
% %     Graphics_reflection1=findobj(h,'color','r');
% %     Graphics_reflection2=findobj(h,'color','b');
%     Fig_handles.Graphics_reflection=[Fig_handles.Graphics_reflection1;Fig_handles.Graphics_reflection2];
%     assignin('base','Graphics_reflection',Fig_handles.Graphics_reflection);
%
% %     reflection_int=findobj(h,'color','y','Marker','.'); % intersection from the reflection
%     assignin('base','reflection_int',Fig_handles.reflection_int);
%
% %     Laser_int_points=findobj(h,'color','k','Marker','.'); % intersection from the reflection
%     assignin('base','Laser_int_points',Fig_handles.Laser_int_points);
%
% %     Normal_lines=findobj(h,'color','y','LineStyle',':'); % intersection from the reflection
% Fig_handles.Normal_lines=[Fig_handles.Normal_lines1;Fig_handles.Normal_lines2];
%     assignin('base','Normal_lines',Fig_handles.Normal_lines);

% Thermal_points=findobj(h,'color','m','Marker','d'); % intersection from the reflection
% assignin('base','Thermal_points',Thermal_points);




% set(Graphics_reflection,'Visible','off');

Laser_head_Ax=Laser_head(1);
Laser_head_Ay=Laser_head(2);

Laser_head_nx=Laser_head(3);
Laser_head_ny=Laser_head(4);
ID=[0 0 0];

Power_Actual=Laser_Power_generator (Laser_head_Ax,Laser_head_Ay,Laser_head_nx,Laser_head_ny,ID);

%% Create a push button

Fig_handles.Graphics_laser=[];
Fig_handles.Graphics_reflection=[];
Fig_handles.h_R=[];
Fig_handles.h_T=[];
Fig_handles.reflection_int=[];
Fig_handles.Laser_int_points=[];
Fig_handles.Normal_lines=[];
Fig_handles.Sub=[];

guidata(gcbo,Fig_handles);

rb1 = uicontrol('Style', 'radiobutton','Position',[20 50 90 40],'Value',true,...
    'String','Rays','Callback', @(rb1,event) Visible_plotButtonPushed(rb1,Fig_handles.Graphics_laser),'Backgroundcolor','w');

btn1 = uicontrol('Style', 'pushbutton', 'String', 'Clear Laser rays',...
    'Position', [20 20 90 20],...
    'Callback', 'delete(Graphics_laser)','Backgroundcolor','w');


rb2 = uicontrol('Style', 'radiobutton','Position',[150 50 90 40],'Value',true,...
    'String','Reflections','Callback', @(rb2,event) Visible_plotButtonPushed(rb2,Fig_handles.Graphics_reflection),'Backgroundcolor','w');

btn2 = uicontrol('Style', 'pushbutton', 'String', 'Clear Reflections',...
    'Position', [150 20 90 20],...
    'Callback', 'delete(Graphics_reflection)','Backgroundcolor','w');



rb3 = uicontrol('Style', 'radiobutton','Position',[280 50 90 40],'Value',true,...
    'String','Mandrel','Callback', @(rb3,event) Visible_plotButtonPushed(rb3,Mandrel_plot),'Backgroundcolor','w');

btn3 = uicontrol('Style', 'pushbutton', 'String', 'Clear Mandrel',...
    'Position', [280 20 90 20],...
    'Callback', 'delete(Mandrel_plot)','Backgroundcolor','w');


rb4 = uicontrol('Style', 'radiobutton','Position',[400 50 90 40],'Value',true,...
    'String','Roller','Callback', @(rb4,event) Visible_plotButtonPushed(rb4,Fig_handles.h_R),'Backgroundcolor','w');

btn4 = uicontrol('Style', 'pushbutton', 'String', 'Clear Roller',...
    'Position', [400 20 90 20],...
    'Callback', 'delete(h_R)','Backgroundcolor','w');


rb5 = uicontrol('Style', 'radiobutton','Position',[520 50 90 40],'Value',true,...
    'String','Tape-Substrate','Callback', @(rb5,event) Visible_plotButtonPushed(rb5,[Fig_handles.Sub,Fig_handles.h_T]),'Backgroundcolor','w');

btn5 = uicontrol('Style', 'pushbutton', 'String', 'Clear Tape-Sub',...
    'Position', [520 20 90 20],...
    'Callback', 'delete(Tape_sub)','Backgroundcolor','w');




rb6 = uicontrol('Style', 'radiobutton','Position',[640 50 90 40],'Value',true,...
    'String','1st reflections','Callback', @(rb6,event) Visible_plotButtonPushed(rb6,Fig_handles.reflection_int),'Backgroundcolor','w');
btn6 = uicontrol('Style', 'pushbutton', 'String', 'Clear reflection_int',...
    'Position', [640 20 90 20],...
    'Callback', 'delete(reflection_int)','Backgroundcolor','w');


rb7 = uicontrol('Style', 'radiobutton','Position',[760 50 90 40],'Value',true,...
    'String','Laser_int','Callback', @(rb7,event) Visible_plotButtonPushed(rb7,Fig_handles.Laser_int_points),'Backgroundcolor','w');

btn7 = uicontrol('Style', 'pushbutton', 'String', 'Clear Laser_int_points',...
    'Position', [760 20 90 20],...
    'Callback', 'delete(Laser_int_points)','Backgroundcolor','w');


%
% rb8 = uicontrol('Style', 'radiobutton','Position',[900 50 90 40],'Value',true,...
%     'String','Thermal_P','Callback', @(rb8,event) Visible_plotButtonPushed(rb8,Thermal_points));

rb8 = uicontrol('Style', 'radiobutton','Position',[900 50 90 40],'Value',true,...
    'String','Normal lines','Callback', @(rb8,event) Visible_plotButtonPushed(rb8,Fig_handles.Normal_lines),'Backgroundcolor','w');


%     string=[]
btn8 = uicontrol('Style', 'pushbutton', 'String', 'Clear Normal_lines',...
    'Position', [900 20 90 20],...
    'Callback', 'delete(Normal_lines)','Backgroundcolor','w');




text_orientation = uicontrol('Style','text','Units','normalized',...
    'String',sprintf('Total Absorbed Energy: \n Tape=%f \n Mandrel=%f \n Roller=%f',0 ,0,0),'Fontsize',10,...
    'Position',[0.8 0.9 0.2 0.08],'Visible','on','Backgroundcolor','w','HorizontalAlignment','right');

% jFigPeer = get(handle(gcf),'JavaFrame');
% jWindow = jFigPeer.fFigureClient.getWindow;
% com.sun.awt.AWTUtilities.setWindowOpacity(jWindow,0.7)

sld_transparency= uicontrol('Style', 'slider',...
        'Min',0,'Max',1,'Value',1,...
       'Units','normalized'...
    ,'Position', [0.09 0.98 0.08 0.02],'Callback',  @(sld_transparency,event) change_transp(sld_transparency)); 


pop_items = uicontrol('Style', 'popupmenu','String', {1:length(nip_point_M_all)},...
    'Units','normalized'...
    ,'Position', [0 0.98 0.08 0.02],...
    'Callback',  @(pop_items,event) show_step_outputs(text_orientation,pop_items,jobname,nip_point_M_all,...
    N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,Rot_Roller_axis_all,H_indentation,tv_def_all,...
    rb1,rb2,rb3,rb4,rb6,rb7,rb8,...
    sld_transparency));


% pop_LaserID = uicontrol('Style', 'popupmenu','String', {'ID-1','ID-2','ID-3'},...
%     'Units','normalized'...
%     ,'Position', [0 0.94 0.08 0.02]    );

text_box_ID = uicontrol('Style', 'edit','String', '0 0 0',...
    'Units','normalized'...
    ,'Position', [0 0.94 0.08 0.02]    );



% 'Callback', @surfzlim


[cdata,map] = imread('complete.jpg');

h2=msgbox('Please select a step number from Pop-up menu at top left corner!',...
    'Notice','custom',cdata,map);
javaFrame    = get(h2,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
set(h2,'color','w');




    function change_transp (sld_transparency)
        data=get(gca,'ch');
        
        for i=1:7
            data(i).Color(4) = get(sld_transparency,'value');
            
        end
        
    end


    function Visible_plotButtonPushed(rb1,Graphics_laser)
        %         x = linspace(0,2*pi,100);
        %         y = sin(x);
        %         plot(ax,x,y)
        if rb1.Value ==false
            set(Graphics_laser,'Visible','off');
        else
            set(Graphics_laser,'Visible','on');
        end
    end

    function show_step_outputs(text_orientation,src,jobname,nip_point_M_all,...
            N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,Rot_Roller_axis_all,H_indentation,tv_def,...
            rb1,rb2,rb3,rb4,rb6,rb7,rb8,sld_transparency)
        val = src.Value;
        str = src.String;
        %         str{val};
        %         disp(['Selection: ' str{val}]);
        
        % Fig_handles = guidata(gcbo);
        
       ID=str2num( get(text_box_ID,'String'));
        Power_Actual=Laser_Power_generator (Laser_head_Ax,Laser_head_Ay,Laser_head_nx,Laser_head_ny,ID);
        
        ss=val;
        tv_def=tv_def_all(ss,:,:);
        
        % jobname='Example0';
        %     fileID1 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Cylinder_ints%d.txt'),jobname,ss),'r');   % file includes the xyz + ID
        %     fileID2 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Tape_ints%d.txt'),jobname,ss),'r');     % file includes the xyz + ID
        %     fileID3 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Roller_ints%d.txt'),jobname,ss),'r');
        
        fileID1 = fopen(strcat(UOT_pathfile, sprintf('Cylinder_ints%d.txt',ss)),'r');
        fileID2 = fopen(strcat(UOT_pathfile, sprintf('Tape_ints%d.txt',ss)),'r');
        fileID3 = fopen(strcat(UOT_pathfile, sprintf('Roller_ints%d.txt',ss)),'r');
        
        
        %%  instead of plotting in each loop, plot all afterwards
        %     fileID4 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\int_xyz%d.txt'),jobname,ss),'r');
        %     fileID5 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Laser_Rays%d.txt'),jobname,ss),'r');
        %     fileID6 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Normal_vectors_Mandrel%d.txt'),jobname,ss),'r');
        %     fileID7 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Reflection_vector_Mandrel%d.txt'),jobname,ss),'r');
        %     fileID8 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Normal_vectors_Tape%d.txt'),jobname,ss),'r');
        %     fileID9 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Reflection_vector_Tape%d.txt'),jobname,ss),'r');
        %
        %     fileID11 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Ref_int_xyz%d.txt'),jobname,ss),'r');
        fileID4 = fopen(strcat(UOT_pathfile, sprintf('int_xyz%d.txt',ss)),'r');
        fileID5 = fopen(strcat(UOT_pathfile, sprintf('Laser_Rays%d.txt',ss)),'r');
        fileID6 = fopen(strcat(UOT_pathfile, sprintf('Normal_vectors_Mandrel%d.txt',ss)),'r');
        fileID7 = fopen(strcat(UOT_pathfile, sprintf('Reflection_vector_Mandrel%d.txt',ss)),'r');
        fileID8 = fopen(strcat(UOT_pathfile, sprintf('Normal_vectors_Tape%d.txt',ss)),'r');
        fileID9 = fopen(strcat(UOT_pathfile, sprintf('Reflection_vector_Tape%d.txt',ss)),'r');
        fileID11 = fopen(strcat(UOT_pathfile, sprintf('Ref_int_xyz%d.txt',ss)),'r');
        
        
        
        Mandrel_E= textscan(fileID1,' %*f %*f %*f %f %f ','Delimiter',',','HeaderLines',1) ;
        Mandrel_E=cell2mat(Mandrel_E);
        
        Mandrel_energy=(Total_energy/(Laser_head_nx*Laser_head_ny))*Power_Actual(Mandrel_E(:,2)).*Mandrel_E(:,1);
        
        
        Tape_E = textscan(fileID2,' %*f %*f %*f %f %f','Delimiter',',','HeaderLines',1) ;
        Tape_E=cell2mat(Tape_E);
        
          Tape_energy=(Total_energy/(Laser_head_nx*Laser_head_ny))*Power_Actual(Tape_E(:,2)).*Tape_E(:,1);
        
        
%         Rolelr_energy = textscan(fileID3,' %*f %*f %*f %f %f','Delimiter',',','HeaderLines',1) ;
%         Rolelr_energy=cell2mat(Rolelr_energy);
        
            Roller_E = textscan(fileID3,' %*f %*f %*f %f %f','Delimiter',',','HeaderLines',1) ;
        Roller_E=cell2mat(Roller_E);
        
          Roller_energy=(Total_energy/(Laser_head_nx*Laser_head_ny))*Power_Actual(Roller_E(:,2)).*Roller_E(:,1);
          
        
        set (text_orientation ,'String',sprintf('Total Absorbed Energy: \n Tape=%f \n Mandrel=%f \n Roller=%f',sum(Tape_energy),sum(Mandrel_energy),sum(Roller_energy)));
        
        
        
        
        
        % h=figure(1);
        
        
        
        % delete([h_R])
        delete([Fig_handles.Sub,Fig_handles.h_R',Fig_handles.h_T,Fig_handles.Graphics_laser,Fig_handles.Graphics_reflection',Fig_handles.reflection_int,Fig_handles.Laser_int_points,Fig_handles.Normal_lines']);
        % plot Mandrel
        % First is cylinder , second is bottom dome part, third is upper Dome part
        
        
        
        Roller_Pos_TV=nip_point_M_all(ss,:,:);
        
        
        tv(1)=Roller_Pos_TV(1)+(-tv_def(1));
        tv(2)=Roller_Pos_TV(2)-tv_def(2);
        tv(3)=Roller_Pos_TV(3)-tv_def(3);
        
        Rot_Roller_axis= reshape(Rot_Roller_axis_all(ss,:,:),[3 3]);
        % output is real center of roller without deformation, end of arrow of
        % the local axis
        [n1,n2,n3,x0,y0,z0,Fig_handles.h_R,Fig_handles.h_T,Center_roller]=Tape_plot_3D_Tshow(N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,Rot_Roller_axis,H_indentation,tv_def);
        
        % assignin('base','h_R',Fig_handles.h_R);
        % assignin('base','h_T',Fig_handles.h_T);
        
        
        mat_size=size(CV_mesh);
        
        X=reshape(CV_mesh(ss,1,:,:),mat_size(3:4));
        Y= reshape(CV_mesh(ss,2,:,:),mat_size(3:4));
        Z=reshape(CV_mesh(ss,3,:,:),mat_size(3:4));
        %         zeroval=zeros(size(Z))
        Fig_handles.Sub=mesh(X,Y,Z,'Facecolor','c','FaceAlpha',0.4);
        
        
        
        XYZ_int = textscan(fileID4,'%f %f %f','Delimiter',',','HeaderLines',0) ;
        XYZ_int=cell2mat(XYZ_int);
        Fig_handles.Laser_int_points=plot3( XYZ_int(:,1), XYZ_int(:,2), XYZ_int(:,3) ,'k.');
        
        XYZ_Rays = textscan(fileID5,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
        XYZ_Rays=cell2mat(XYZ_Rays);
        Fig_handles.Graphics_laser= plot3( XYZ_Rays(:,1), XYZ_Rays(:,2), XYZ_Rays(:,3) ,'g:','Linewidth',0.001 );
        
        Normal_vectors_Mandrel = textscan(fileID6,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
        Normal_vectors_Mandrel=cell2mat(Normal_vectors_Mandrel);
        Fig_handles.Normal_lines1=plot3( Normal_vectors_Mandrel(:,1), Normal_vectors_Mandrel(:,2), Normal_vectors_Mandrel(:,3) ,'y:');
        
        
        Reflection_vector_Mandrel = textscan(fileID7,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
        Reflection_vector_Mandrel=cell2mat(Reflection_vector_Mandrel);
        Fig_handles.Graphics_reflection1=plot3( Reflection_vector_Mandrel(:,1), Reflection_vector_Mandrel(:,2), Reflection_vector_Mandrel(:,3) ,'r--','Linewidth',0.01 );
        
        
        Normal_vectors_Tape = textscan(fileID8,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
        Normal_vectors_Tape=cell2mat(Normal_vectors_Tape);
        Fig_handles.Normal_lines2= plot3( Normal_vectors_Tape(:,1), Normal_vectors_Tape(:,2), Normal_vectors_Tape(:,3) ,'y:');
        
        
        Reflection_vector_Tape = textscan(fileID9,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
        Reflection_vector_Tape=cell2mat(Reflection_vector_Tape);
        Fig_handles.Graphics_reflection2=plot3( Reflection_vector_Tape(:,1), Reflection_vector_Tape(:,2), Reflection_vector_Tape(:,3) ,'b--','Linewidth',0.01 );
        
        Ref_int_xyz = textscan(fileID11,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
        Ref_int_xyz=cell2mat(Ref_int_xyz);
        Fig_handles.reflection_int=plot3( Ref_int_xyz(:,1), Ref_int_xyz(:,2), Ref_int_xyz(:,3) ,'y.');
        
        
        
        
        
        %% Changing the graphical object, Visibility or deleting
        
        %     Graphics_laser=findobj(h,'color','g','LineStyle',':');
        assignin('base','Graphics_laser',Fig_handles.Graphics_laser);
        % set(Graphics_laser,'Visible','off');
        %     Graphics_reflection1=findobj(h,'color','r');
        %     Graphics_reflection2=findobj(h,'color','b');
        Fig_handles.Graphics_reflection=[Fig_handles.Graphics_reflection1;Fig_handles.Graphics_reflection2];
        assignin('base','Graphics_reflection',Fig_handles.Graphics_reflection);
        
        %     reflection_int=findobj(h,'color','y','Marker','.'); % intersection from the reflection
        assignin('base','reflection_int',Fig_handles.reflection_int);
        
        %     Laser_int_points=findobj(h,'color','k','Marker','.'); % intersection from the reflection
        assignin('base','Laser_int_points',Fig_handles.Laser_int_points);
        assignin('base','h_R',Fig_handles.h_R);
        
        
        %     Normal_lines=findobj(h,'color','y','LineStyle',':'); % intersection from the reflection
        Fig_handles.Normal_lines=[Fig_handles.Normal_lines1;Fig_handles.Normal_lines2];
        assignin('base','Normal_lines',Fig_handles.Normal_lines);
        
        % Thermal_points=findobj(h,'color','m','Marker','d'); % intersection from the reflection
        % assignin('base','Thermal_points',Thermal_points);
        
        Fig_handles.Tape_sub=[Fig_handles.Sub,Fig_handles.h_T];
        assignin('base','Tape_sub',Fig_handles.Tape_sub);
        
        
        
        
        % set(Graphics_reflection,'Visible','off');
        
        
        %        rb5 = uicontrol('Style', 'radiobutton','Position',[520 50 90 40],'Value',true,...
        %         'String','Tape-Substrate','Callback', @(rb5,event) Visible_plotButtonPushed(rb5,[Fig_handles.Sub,Fig_handles.h_T]),'Backgroundcolor','w');
        
        %     btn5 = uicontrol('Style', 'pushbutton', 'String', 'Clear Tape-Sub',...
        %         'Position', [520 20 90 20],...
        %         'Callback', 'delete(Tape_sub)','Backgroundcolor','w');
        
        
        
        set(rb1,'Callback', @(rb1,event) Visible_plotButtonPushed(rb1,Fig_handles.Graphics_laser));
        set(rb2,'Callback', @(rb2,event) Visible_plotButtonPushed(rb2,Fig_handles.Graphics_reflection));
        % mandrel is not updating
        %set(rb3,'Callback', @(rb3,event) Visible_plotButtonPushed(rb1,Graphics_laser));
        set(rb4,'Callback', @(rb4,event) Visible_plotButtonPushed(rb4,Fig_handles.h_R));
        
        set(rb5,'Callback', @(rb5,event) Visible_plotButtonPushed(rb5,[Fig_handles.Sub,Fig_handles.h_T]));
        
        set(rb6,'Callback', @(rb6,event) Visible_plotButtonPushed(rb6,Fig_handles.reflection_int));
        set(rb7,'Callback', @(rb7,event) Visible_plotButtonPushed(rb7,Fig_handles.Laser_int_points));
        set(rb8,'Callback', @(rb8,event) Visible_plotButtonPushed(rb8,Fig_handles.Normal_lines));
        
        
        
        fclose(fileID1);
        fclose(fileID2);
        fclose(fileID3);
        fclose(fileID4);
        fclose(fileID5);
        fclose(fileID6);
        fclose(fileID7);
        fclose(fileID8);
        fclose(fileID9);
        fclose(fileID11);
        
        
        
        
        
        % 4 transparency of Setup.fig
data=get(gca,'ch');

for i=1:7
data(i).Color(4) = get(sld_transparency,'value');

end
        
        
        guidata(gcbo,Fig_handles) ;
        
    end
end
