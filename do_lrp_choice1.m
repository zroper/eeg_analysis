%attempts to bin trials

load erp;
erp.trial = NaN;

conditions_base = {'Target1';'Target2';'Target3'};

conditions = {'Lat_Target1';'Lat_Target2';'Lat_Target3';...
              'Mid_Target1';'Mid_Target2';'Mid_Target3';...
              };

conditions_LR = {'UTarget1';'LTarget1';'BTarget1';'RTarget1';...
                 'UTarget2';'LTarget2';'BTarget2';'RTarget2';...
                 'UTarget3';'LTarget3';'BTarget3';'RTarget3';...
              };

%conditions_LR = {'LShort1';'LShort2';'LShort4';'LShort8';'LLong1';'LLong2';'LLong4';'LLong8';'RShort1';'RShort2';'RShort4';'RShort8';'RLong1';'RLong2';'RLong4';'RLong8'};
channels = {'VEOG';'HEOG';'F3';'F4';'C3';'C4';'P3';'P4';'PO3';'PO4';'O1';'O2';'OL';'OR';'T3';'T4';'T5';'T6';'Fz';'Cz';'Pz';'L_mas'};
%channels_r = {'F4';'F3';'C4';'C3';'P4';'P3';'PO4';'PO3';'O2';'O1';'OR';'OL';'T4';'T3';'T6';'T5';'Fz';'Cz';'Pz'};

%cue_code = {7 9};
%mem_codes = {11 12 13 14 16 18};


baseline_window = [101:200];%instep
target_window = [201:300];
%mem_window = [101:350];
erp.baseline = [101:200];
erp.mem_window = [200,99];
%erp.mem_window = [100,249];



erp.trialCodes = zeros(1,length(erp.eventCodes));
erp.trialNumCodes = zeros(1,length(erp.eventCodes));
trial_counter = 0;
%GOTTA SORT OUT THE CONDITION ORDER!!!!!!!!
%IMPLEMENT TRIALCODES based on eventcoes!!!!!!!!!
for ec = 1:1:length(erp.eventCodes)-1
    if (erp.eventCodes(ec) >200 && erp.eventCodes(ec) < 209)
        if erp.eventCodes(ec+1) > 100 && erp.eventCodes(ec+1)<197
            trial_counter = trial_counter+1;
            if (erp.eventCodes(ec+2) == 10)
                erp.trialCodes(ec+3) = erp.eventCodes(ec+3)-10;
                erp.trialNumCodes(ec+3) = trial_counter;
            elseif (erp.eventCodes(ec+2) == 20)
                erp.trialCodes(ec+3) = erp.eventCodes(ec+3)-20+4;
                erp.trialNumCodes(ec+3) = trial_counter;
            elseif (erp.eventCodes(ec+2) == 30)
                erp.trialCodes(ec+3) = erp.eventCodes(ec+3)-30+8;
                erp.trialNumCodes(ec+3) = trial_counter;
            end;

        end;
    end;
  
end;

    
%IMPLEMENT TRIALCODES based on eventcoes!!!!!!!!! 

counters = ones(1,length(conditions_LR));
art_free_counter = zeros(1,length(conditions_LR));
art_counter = zeros(1,length(conditions_LR));
blink_counter = zeros(1,length(conditions_LR));
eMove_counter = zeros(1,length(conditions_LR));
blocking_counter = zeros(1,length(conditions_LR));

for ec  = 1:1:length(erp.eventCodes)

    if erp.trialCodes(ec) > 0 %found a trial segment!
            
            fieldname = (char(conditions_LR(erp.trialCodes(ec),:)));
                
            
            trialNumname = [fieldname,'_Num'];
            
            
            artifact = 0;
           
            
            %adjust the timewindow based on condition
            if (erp.trialCodes(ec) >=1 && erp.trialCodes(ec) <13)%Cue
                pre_timepoint = erp.mem_window(1);
                post_timepoint = erp.mem_window(2);
                pre_timepoint_hilb = erp.mem_window(1)+300;%to give ample time for filters!
                post_timepoint_hilb = erp.mem_window(2)+300;%to give ample time for filters!
            end;
                
            
            %check the time range if it's artifact free
            
            if sum(erp.arf.blink(erp.eventTimes(ec)-pre_timepoint:erp.eventTimes(ec)+post_timepoint))>0
                artifact = 1;
                blink_counter(erp.trialCodes(ec))=blink_counter(erp.trialCodes(ec))+1; 
            elseif sum(erp.arf.eMove(erp.eventTimes(ec)-pre_timepoint:erp.eventTimes(ec)+post_timepoint))>0
                artifact = 1;
                eMove_counter(erp.trialCodes(ec))=eMove_counter(erp.trialCodes(ec))+1; 
            elseif sum(sum(erp.arf.blocking(:,erp.eventTimes(ec)-pre_timepoint:erp.eventTimes(ec)+post_timepoint)))>0
                artifact = 1;
                blocking_counter(erp.trialCodes(ec))=blocking_counter(erp.trialCodes(ec))+1; 
            end
            
            if artifact == 0
                art_free_counter(erp.trialCodes(ec))=art_free_counter(erp.trialCodes(ec))+1; 
                erp.trial.(fieldname)(counters(erp.trialCodes(ec)),:,:)= erp.filtered_data(:,erp.eventTimes(ec)-pre_timepoint:erp.eventTimes(ec)+post_timepoint);
                erp.trial_hilb.(fieldname)(counters(erp.trialCodes(ec)),:,:)= erp.filtered_data(:,erp.eventTimes(ec)-pre_timepoint_hilb:erp.eventTimes(ec)+post_timepoint_hilb);
                erp.reject_trial.(fieldname)(counters(erp.trialCodes(ec)),:,:)= nan(22,pre_timepoint+post_timepoint+1);
                erp.reject_trial_hilb.(fieldname)(counters(erp.trialCodes(ec)),:,:)= nan(22,pre_timepoint_hilb+post_timepoint_hilb+1);
            elseif artifact >0
                art_counter(erp.trialCodes(ec))=art_counter(erp.trialCodes(ec))+1; 
                erp.trial.(fieldname)(counters(erp.trialCodes(ec)),:,:)= nan(22,pre_timepoint+post_timepoint+1);
                erp.trial_hilb.(fieldname)(counters(erp.trialCodes(ec)),:,:)= nan(22,pre_timepoint_hilb+post_timepoint_hilb+1);
                
                erp.reject_trial.(fieldname)(counters(erp.trialCodes(ec)),:,:)= erp.filtered_data(:,erp.eventTimes(ec)-pre_timepoint:erp.eventTimes(ec)+post_timepoint);
                erp.reject_trial_hilb.(fieldname)(counters(erp.trialCodes(ec)),:,:)= erp.filtered_data(:,erp.eventTimes(ec)-pre_timepoint_hilb:erp.eventTimes(ec)+post_timepoint_hilb);
            end;
            erp.trialNum.(fieldname)(counters(erp.trialCodes(ec)))=erp.trialNumCodes(ec);
            counters(erp.trialCodes(ec))= counters(erp.trialCodes(ec))+1;
              
        
    end;


end;

%now baseline it!

for condition = 1:1:length(conditions_LR)

    fieldname = (char(conditions_LR(condition,:)));
    
    baseline = repmat(squeeze(mean(erp.trial.(fieldname)(:,:,erp.baseline),3)),[1,1,size(erp.trial.(fieldname),3)]);     
    erp.trial.(fieldname)= erp.trial.(fieldname)-baseline;
    
    reject_baseline = repmat(squeeze(mean(erp.reject_trial.(fieldname)(:,:,erp.baseline),3)),[1,1,size(erp.trial.(fieldname),3)]);     
    erp.reject_trial.(fieldname)= erp.reject_trial.(fieldname)-reject_baseline;

end;

%Now create the erp&ave files
choice1_erp = struct();
reject_choice1_erp = struct();
choice1_erp_ave = struct();

channels = {'VEOG';'HEOG';'F3';'F4';'C3';'C4';'P3';'P4';'PO3';'PO4';'O1';'O2';'OL';'OR';'T3';'T4';'T5';'T6';'Fz';'Cz';'Pz';};
channels_r = {'VEOG';'HEOG';'F4';'F3';'C4';'C3';'P4';'P3';'PO4';'PO3';'O2';'O1';'OR';'OL';'T4';'T3';'T6';'T5';'Fz';'Cz';'Pz';};
channel_r_index = [1 2 4 3 6 5 8 7 10 11 12 11 14 13 16 15 18 17 19 20 21];
for condition = 1:1:length(conditions_base)
    base_condition_name = (char(conditions_base(condition,:)));
    u_condition_name = ['U',base_condition_name];
    l_condition_name = ['L',base_condition_name];
    r_condition_name = ['R',base_condition_name];
    b_condition_name = ['B',base_condition_name];
    
    lat_condition_name = ['Lat_',base_condition_name];
    mid_condition_name = ['Mid_',base_condition_name];
    

    %trial_name = [condition_name,'_trial'];
    %choice1_erp.(trial_name)=[squeeze(erp.trialNum.(condition_name))];

    
    %reject_choice1_erp.(trial_name)=[squeeze(erp.trialNum.(condition_name))];

    
    for channel = 1:1:length(channels)
        channel_name = char(channels(channel,:));
        channel_name_r = char(channels_r(channel,:));
        lat_fieldname = [lat_condition_name,'_',channel_name];
        mid_fieldname = [mid_condition_name,'_',channel_name];
        lat_fieldname_hilb = [lat_condition_name,'_',channel_name,'_hilb'];
        mid_fieldname_hilb = [mid_condition_name,'_',channel_name,'_hilb'];
        %fieldname_l = [condition_name_l,'_',channel_name];
        %fieldname_r = [condition_name_r,'_',channel_name_r];
        

        choice1_erp.(lat_fieldname) = [squeeze(erp.trial.(l_condition_name)(:,channel,:));squeeze(erp.trial.(r_condition_name)(:,channel_r_index(channel),:))];
        choice1_erp.(mid_fieldname) = [squeeze(erp.trial.(u_condition_name)(:,channel,:));squeeze(erp.trial.(b_condition_name)(:,channel,:))];
        choice1_erp.(lat_fieldname_hilb) = [squeeze(erp.trial_hilb.(l_condition_name)(:,channel,:));squeeze(erp.trial_hilb.(r_condition_name)(:,channel_r_index(channel),:))];
        choice1_erp.(mid_fieldname_hilb) = [squeeze(erp.trial_hilb.(u_condition_name)(:,channel,:));squeeze(erp.trial_hilb.(b_condition_name)(:,channel,:))];

        reject_choice1_erp.(lat_fieldname) = [squeeze(erp.reject_trial.(l_condition_name)(:,channel,:));squeeze(erp.reject_trial.(r_condition_name)(:,channel_r_index(channel),:))];
        reject_choice1_erp.(mid_fieldname) = [squeeze(erp.reject_trial.(u_condition_name)(:,channel,:));squeeze(erp.reject_trial.(b_condition_name)(:,channel,:))];
        
        choice1_erp_ave.(lat_fieldname) = nanmean(choice1_erp.(lat_fieldname)(:,:),1);
        choice1_erp_ave.(mid_fieldname) = nanmean(choice1_erp.(mid_fieldname)(:,:),1);


    end;
    
    
    
end;

%summary
for channel = 1:1:length(channels)
    
end;

%summarize artifact
for condition = 1:1:length(conditions_base)
    base_condition_name = (char(conditions_base(condition,:)));
    lat_condition_name = ['Lat_',base_condition_name];
    mid_condition_name = ['Mid_',base_condition_name];
    
    lat_art_free_name = [lat_condition_name,'_art_free'];
    lat_blink_name = [lat_condition_name,'_blink'];
    lat_eMove_name = [lat_condition_name,'_eMove'];
    lat_blocking_name = [lat_condition_name,'_blocking'];
    
    mid_art_free_name = [mid_condition_name,'_art_free'];
    mid_blink_name = [mid_condition_name,'_blink'];
    mid_eMove_name = [mid_condition_name,'_eMove'];
    mid_blocking_name = [mid_condition_name,'_blocking'];
    
    choice1_art.(lat_art_free_name)= (art_free_counter((condition-1)*4+1)+art_free_counter((condition-1)*4+4))/(art_free_counter((condition-1)*4+1)+art_free_counter((condition-1)*4+4)+art_counter((condition-1)*4+1)+art_counter((condition-1)*4+4));
    choice1_art.(lat_blink_name)= (blink_counter((condition-1)*4+1)+blink_counter((condition-1)*4+4));
    choice1_art.(lat_eMove_name)= (eMove_counter((condition-1)*4+1)+eMove_counter((condition-1)*4+4));
    choice1_art.(lat_blocking_name)= (blocking_counter((condition-1)*4+1)+blocking_counter((condition-1)*4+4));
    
    choice1_art.(mid_art_free_name)= (art_free_counter((condition-1)*4+2)+art_free_counter((condition-1)*4+3))/(art_free_counter((condition-1)*4+2)+art_free_counter((condition-1)*4+3)+art_counter((condition-1)*4+2)+art_counter((condition-1)*4+3));
    choice1_art.(mid_blink_name)= (blink_counter((condition-1)*4+2)+blink_counter((condition-1)*4+3));
    choice1_art.(mid_eMove_name)= (eMove_counter((condition-1)*4+2)+eMove_counter((condition-1)*4+3));
    choice1_art.(mid_blocking_name)= (blocking_counter((condition-1)*4+2)+blocking_counter((condition-1)*4+3));
    
end;


%save('erp.mat','erp');
save('choice1_erp.mat','choice1_erp')
%save('reject_enc_erp.mat','reject_enc_erp')
save('choice1_erp_ave.mat','choice1_erp_ave');
save('choice1_art.mat','choice1_art');


 