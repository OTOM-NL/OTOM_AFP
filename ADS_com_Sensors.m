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

err=[];

[m,n]=size(ADS_Vars);

for ii=1:m
    
    if ADS_Vars{ii}
        
        
        % for getting data name from TwinCat
        
        
        try
            
            hVar = adsClt.CreateVariableHandle( strcat(TwinCat_file_name,'.',ADS_Vars{ii,2})); %name of the variable
        catch err
            adsClt.Dispose();
%             disp(err.message);
             msgbox ('Error');
        end
        
              
        try
            %% Create AdsStream instance
            
            % for Int Datatype
            
            Var=adsClt.ReadSymbolInfo(strcat(TwinCat_file_name,'.',ADS_Vars{ii,2}));
            Data_type=char(Var.Datatype);
            
            
            switch  Data_type(6)
                case 'I'
            
               dataStreamRead = AdsStream(2); %number of y
            binRead = AdsBinaryReader(dataStreamRead);
            adsClt.Read(hVar,dataStreamRead);
                  readValue(ii) = binRead.ReadInt16();
                  
                msgbox(sprintf('Value= %f',readValue(ii)));
                        
            valueToWrite = 3+ii;
            
            dataStreamWrite= AdsStream(2);
            binWrite = AdsBinaryWriter(dataStreamWrite);
            binWrite.Write(int16(valueToWrite));
            
            adsClt.Write(hVar,dataStreamWrite);
%             adsClt.DeleteVariableHandle(hVar);
        
                case 'R'
                  
                       % for Real Datatype
            
            dataStreamRead = AdsStream(8); %number of y
            binRead = AdsBinaryReader(dataStreamRead);
            adsClt.Read(hVar,dataStreamRead);
                  readValue(ii) = binRead.ReadDouble();
                  
                  
                     msgbox(sprintf('Value= %f',readValue(ii)));
                  
                  
                    valueToWrite = 1.025+ii;
            
            dataStreamWrite= AdsStream(8);
            binWrite = AdsBinaryWriter(dataStreamWrite);
            binWrite.Write(double(valueToWrite));
            
            adsClt.Write(hVar,dataStreamWrite);
            
%             adsClt.DeleteVariableHandle(hVar);
            end


            
         
            
            %%
            % program should be run here ??!
            
       
      
            
            
        catch err
            adsClt.Dispose();
%             disp(['Error: ' err.message]);
             msgbox ('Error');
        end
        
        
        
        
%         ADS_Vars{ii,1};
    end
    
end



%% Dispose ADS client
adsClt.Dispose();









