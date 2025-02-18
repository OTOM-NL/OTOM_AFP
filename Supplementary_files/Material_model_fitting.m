
fit_degree=4;

 ft = fittype( strcat('poly','4') );
 

    %%
%     fitresult_Kx,fitresult_Ky,fitresult_Cp,fitresult_rho
    
    
    
    Density=Density(2:end);
    specificheat=specificheat(2:end);
    Temperature= Temperature(2:end);
ThermalconductivityAxial= ThermalconductivityAxial(2:end);
Thermalconductivitytransverse=Thermalconductivitytransverse(2:end);


%%

% Density

   % Fit model to data.
    % [fitresult, gof] = fit( force, normal_dis, ft )
    [fitresult_rho] = fit( Temperature, Density, ft );
    
    % plot(fitresult,force, normal_dis,'predfunc');
    figure(61);
    plot(fitresult_rho,Temperature, Density,'fit','residuals');
    legend Location SouthWest;
    title('residuals');
    subplot(2,1,1);
    legend Location NorthWest;
    title(strcat('rho: ',formula(fitresult_rho)));
    hold on;
    
    
    %%
    
    
% specificheat

   % Fit model to data.
    % [fitresult, gof] = fit( force, normal_dis, ft )
    [fitresult_Cp] = fit( Temperature, specificheat, ft );
    
    % plot(fitresult,force, normal_dis,'predfunc');
    figure(62);
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
    plot(fitresult_Ky,Temperature, Thermalconductivitytransverse,'fit','residuals');
    legend Location SouthWest;
    title('residuals');
    subplot(2,1,1);
    legend Location NorthWest;
    title(strcat('Ky: ',formula(fitresult_Ky)));
    hold on;
    
