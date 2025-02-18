


function save_fitting (Poly_fit_degree,uit,fig,Material_name)
% fit_degree=str2double(Poly_fit_degree.String);

 ft = fittype( strcat('poly',Poly_fit_degree) );
 

    %%
%     fitresult_Kx,fitresult_Ky,fitresult_Cp,fitresult_rho
    
    
    
    Density=uit(:,3);
    specificheat=uit(:,2);
    Temperature= uit(:,1);
%     ThermalconductivityAxial
ThermalconductivityAxial= uit(:,4);
Thermalconductivitytransverse=uit(:,5);
Thermalconductivitytransverse_z=uit(:,6);

%%

% Density

   % Fit model to data.
    % [fitresult, gof] = fit( force, normal_dis, ft )
    [fitresult_rho] = fit( Temperature, Density, ft );
    
    % plot(fitresult,force, normal_dis,'predfunc');
    figure(66);
    javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

    plot(fitresult_rho,Temperature, Density,'fit','residuals');
    legend Location SouthWest;
    title('residuals');
    subplot(2,1,1);
    legend Location NorthWest;
    title(strcat('\rho: ',formula(fitresult_rho)));
    hold on;
    
    
    %%
    
    
% specificheat

   % Fit model to data.
    % [fitresult, gof] = fit( force, normal_dis, ft )
    [fitresult_Cp] = fit( Temperature, specificheat, ft );
    
    % plot(fitresult,force, normal_dis,'predfunc');
    figure(62);
     javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    plot(fitresult_Cp,Temperature, specificheat,'fit','residuals');
    legend Location SouthWest;
    title('residuals');
    subplot(2,1,1);
    legend Location NorthWest;
    title(strcat('Cp: ',formula(fitresult_Cp)));
    hold on;
    
    
    %%
    
        
% ThermalconductivityAxial

   % Fit model to data.
    % [fitresult, gof] = fit( force, normal_dis, ft )
    [fitresult_Kx] = fit( Temperature, ThermalconductivityAxial, ft );
    
    % plot(fitresult,force, normal_dis,'predfunc');
   figure(63);
   javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

    plot(fitresult_Kx,Temperature, ThermalconductivityAxial,'fit','residuals');
    legend Location SouthWest;
    title('residuals');
    subplot(2,1,1);
    legend Location NorthWest;
    title(strcat('Kx: ',formula(fitresult_Kx)));
    hold on;
    
    
     %%
%     Thermalconductivitytransverse


    % Fit model to data.
    % [fitresult, gof] = fit( force, normal_dis, ft )
    [fitresult_Ky] = fit( Temperature, Thermalconductivitytransverse, ft );
    
    % plot(fitresult,force, normal_dis,'predfunc');
     figure(64);
      javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    plot(fitresult_Ky,Temperature, Thermalconductivitytransverse,'fit','residuals');
    legend Location SouthWest;
    title('residuals');
    subplot(2,1,1);
    legend Location NorthWest;
    title(strcat('Ky: ',formula(fitresult_Ky)));
    hold on;
    
    [fitresult_Kz] = fit( Temperature, Thermalconductivitytransverse_z, ft );
    
    % plot(fitresult,force, normal_dis,'predfunc');
     figure(65);
      javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    plot(fitresult_Kz,Temperature, Thermalconductivitytransverse_z,'fit','residuals');
    legend Location SouthWest;
    title('residuals');
    subplot(2,1,1);
    legend Location NorthWest;
    title(strcat('Kz: ',formula(fitresult_Kz)));
    hold on;
    
    % Save fitting data
%         dir2save=UOT_pathfile;
     save(strcat('.\Supp_files\', Material_name,'.mat'),'fitresult_Kz','fitresult_Ky','fitresult_Kx',...
     'fitresult_Cp','fitresult_rho'    );
close(fig);

