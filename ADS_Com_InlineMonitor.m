function err=ADS_Com_InlineMonitor (ADS_Vars,ADS_fileName,ADS_PathName,TwinCat_file_name)

% This function
% 1- recieve variable from TwinCat  >> input collection
% 2-Compute the temperature
% 3- Send back variables and new Temperature to TwinCat >> feedback section



%% Input : Define variables, should be got from TwinCAT
% data recieve from TwinCAt , process in Matlab and send to TwinCAt again !
asm = NET.addAssembly(strcat(ADS_PathName,ADS_fileName));
import TwinCAT.Ads.*;

% Create instance of class TcAdsClient
adsClt=TwinCAT.Ads.TcAdsClient;

% tcClient.Connect(851)
adsClt.Connect(851);

% Defining variables: for sensors

% 1 variable to indicate the process is running/Stop

err=[];

[m,n]=size(ADS_Vars);

readValue=zeros(m,1);
hVar=zeros(m,1);
Data_type=cell(m,1);


%%
% for saving handles and initialization

Active_Var=0;
for ii=1:m
    
    if ADS_Vars{ii}
        
        Active_Var=Active_Var+1;
        try
            
            hVar(ii) = adsClt.CreateVariableHandle( strcat(TwinCat_file_name,'.',ADS_Vars{ii,2})); %name of the variable
        catch err
            adsClt.Dispose();
            %             disp(err.message);
            msgbox ('Error');
        end
        
        
        Var=adsClt.ReadSymbolInfo(strcat(TwinCat_file_name,'.',ADS_Vars{ii,2}));
        Data_type{ii}=char(Var.Datatype);
        
        
    else
        break
    end
    
end








%%

p=1; % just to enter the loop

while p >0
    
    for ii=1:Active_Var
        
        
        % for getting data Value from TwinCat
        % Create AdsStream instance
        
        switch  Data_type{ii}(6)
            case 'I'
                %  >>> These sentences to read value every time
                dataStreamRead = AdsStream(2); %number of y
                binRead = AdsBinaryReader(dataStreamRead);
                adsClt.Read(hVar(ii),dataStreamRead);
                readValue(ii) = binRead.ReadInt16();
                
                msgbox(sprintf('Value= %f',readValue(ii)));
                
                
                
            case 'R'
                
                % for Real Datatype
                
                %  >>> These sentences to read value every time
                dataStreamRead = AdsStream(8); %number of y
                binRead = AdsBinaryReader(dataStreamRead);
                adsClt.Read(hVar(ii),dataStreamRead);
                readValue(ii) = binRead.ReadDouble();
                
                
                msgbox(sprintf('Value= %f',readValue(ii)));
                
                
        end
        
        
        
        
        
    end
    
    %%
    % program should be run here ??!
    
    
    
    
    
    
    %% Writting in TwinCat variables using output of computational model
    
    for ii=1:Active_Var
        
        
        % for writting data Value from TwinCat
        % Create AdsStream instance
        
        switch  Data_type{ii}(6)
            case 'I'
                % for integer value
                %  >>> These sentences to write value every time
                
                valueToWrite = 3+ii;
                dataStreamWrite= AdsStream(2);
                binWrite = AdsBinaryWriter(dataStreamWrite);
                binWrite.Write(int16(valueToWrite));
                adsClt.Write(hVar(ii),dataStreamWrite);
                
                
                %             adsClt.DeleteVariableHandle(hVar);
                
            case 'R'
                % for Real Datatype
                %  >>> These sentences to write value every time
                
                valueToWrite = 1.025+ii;
                dataStreamWrite= AdsStream(8);
                binWrite = AdsBinaryWriter(dataStreamWrite);
                binWrite.Write(double(valueToWrite));
                adsClt.Write(hVar(ii),dataStreamWrite);
                
                %             adsClt.DeleteVariableHandle(hVar);
        end
        
        
        
        
        
    end
    
    
    
    
    
    
    p=readValue(1);
end

%% Dispose ADS client
adsClt.Dispose();









