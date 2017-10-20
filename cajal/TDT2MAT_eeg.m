cd ~/Desktop/Ca_data/TDT_Behavior

contents = dir;
files = {contents(~[contents.isdir]).name};
idx   = ~cellfun('isempty',regexp(files,'\d+\.mat'));
files = files(idx);

pre_dir = contents(1).folder;

for fid = 1:length(files)
    full_name = files(fid);
    
    cd eeg_analysis
    
    current_file = strtok(full_name,'.');
    
    mkdir(current_file{1});
    
    source_path = strcat(pre_dir, '/', full_name);
    
    des_dir = strcat(pre_dir, '/eeg_analysis/', current_file);
    
    %movefile(full_name, des_dir)
    load(source_path{1})
    
    cd(des_dir{1})
    
    erp = struct();
    try
        erp.srate = Cajal_data.eeg8.streams.EEG8.fs;
        erp.eventCodes = Cajal_data.evnt.scalars.EVNT.data;
        erp.eventTimes = Cajal_data.evnt.scalars.EVNT.ts*1000;
        
        erp.data=cat(1,Cajal_data.eeg8.streams.EEG8.data,Cajal_data.eeg16.streams.EEG6.data);
        erp.pnts = length(erp.data);
        erp.rateAcq = 1/erp.srate;
    catch
    end
    
    save('erp.mat','erp')
    
    clear Cajal_data
    
    cd ../../
    
    
end