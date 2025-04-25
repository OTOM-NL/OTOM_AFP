function err=ADS_com_Sensors (ADS_Vars,ADS_fileName,ADS_PathName,TwinCat_file_name)

%% Input : Define variables, should be got from TwinCAT
% data recieve from TwinCAt , process in Matlab and send to TwinCAt again !
asm = NET.addAssembly(strcat(ADS_PathName,ADS_fileName));
import TwinCAT.Ads.*;

% Create instance of class TcAdsClient
adsClt=TwinCAT.Ads.TcAdsClient;

% tcClient.Connect(851)
adsClt.Connect(851);

% Defining variables: for sensors

% Var3=adsClt.ReadSymbolInfo(strcat(TwinCat_file_name,'.',ADS_Vars{ii,2}))




% Var3.ReadDouble
