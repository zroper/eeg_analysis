%attempts to bin trials

load erp;
erp.trial = NaN;
conditions = {'Mem0';'Mem1';'Mem2';'Mem3';...
              };
%conditions_LR = {'LShort1';'LShort2';'LShort4';'LShort8';'LLong1';'LLong2';'LLong4';'LLong8';'RShort1';'RShort2';'RShort4';'RShort8';'RLong1';'RLong2';'RLong4';'RLong8'};
channels = {'VEOG';'HEOG';'F3';'F4';'C3';'C4';'P3';'P4';'PO3';'PO4';'O1';'O2';'OL';'OR';'T3';'T4';'T5';'T6';'Fz';'Cz';'Pz';'L_mas'};
%channels_r = {'F4';'F3';'C4';'C3';'P4';'P3';'PO4';'PO3';'O2';'O1';'OR';'OL';'T4';'T3';'T6';'T5';'Fz';'Cz';'Pz'};

%cue_code = {7 9};
%mem_codes = {11 12 13 14 16 18};


baseline_window = [1:100];%instep
mem_window = [101:600];
%mem_window = [101:350];
erp.baseline = [1:100];
erp.mem_window = [100,499];
%erp.mem_window = [100,249];



erp.trialCodes = zeros(1,length(erp.eventCodes));
erp.trialNumCodes = zeros(1,length(erp.eventCodes));
erp.picCodes = zeros(1,length(erp.eventCodes));
erp.cueCodes = zeros(1,length(erp.eventCodes));
trial_counter = 0;
%GOTTA SORT OUT THE CONDITION ORDER!!!!!!!!
%IMPLEMENT TRIALCODES based on eventcoes!!!!!!!!!
for ec = 1:1:length(erp.eventCodes)-5
    if (erp.eventCodes(ec) >200 && erp.eventCodes(ec) < 211)
        if erp.eventCodes(ec+1) > 100 && erp.eventCodes(ec+1)<191
            trial_counter = trial_counter+1;
            if (erp.eventCodes(ec+2) > 99 && erp.eventCodes(ec+2) <104)
                erp.trialCodes(1,ec+2) = erp.eventCodes(ec+2)-99;
                erp.cueCodes(1,ec+2) = erp.cueCodes(ec+2)-99;
                erp.trialNumCodes(1,ec+2) = trial_counter;
            elseif (erp.eventCodes(ec+2) > 29 && erp.eventCodes(ec+2) <34)%double code!
                %erp.trialCodes(1,ec+2) = erp.eventCodes(ec+)-99;
                %erp.cueCodes(1,ec+2) = erp.cueCodes(ec+2)-99;
                erp.trialNumCodes(1,ec+2) = trial_counter;
            end;

                
            if (erp.eventCodes(ec+3) > 29 && erp.eventCodes(ec+3) < 34)
                if (erp.eventCodes(ec+4) > 19 && erp.eventCodes(ec+4) < 30)
                    if (erp.eventCodes(ec+5) > 9 && erp.eventCodes(ec+5) < 20)
                        if erp.eventCodes(ec)<206
                            erp.picCodes(1,ec+2) = (erp.eventCodes(ec+3)-30)*100+(erp.eventCodes(ec+4)-20)*10+(erp.eventCodes(ec+5)-10);
                            %erp.picCodes(1,ec+3) = (erp.eventCodes(ec+4)-30)*100+(erp.eventCodes(ec+5)-20)*10+(erp.eventCodes(ec+6)-10);
                        elseif erp.eventCodes(ec) >205
                            erp.picCodes(1,ec+2) = 1000+(erp.eventCodes(ec+3)-30)*100+(erp.eventCodes(ec+4)-20)*10+(erp.eventCodes(ec+5)-10);
                            %erp.picCodes(1,ec+3) = 1000+(erp.eventCodes(ec+4)-30)*100+(erp.eventCodes(ec+5)-20)*10+(erp.eventCodes(ec+6)-10);
                        end;
                    end;
                end;
            end;

        end;
    end;
  
end;

    
%IMPLEMENT TRIALCODES based on eventcoes!!!!!!!!! 

counters = ones(1,length(conditions));
art_free_counter = zeros(1,length(conditions));
art_counter = zeros(1,length(conditions));
blink_counter = zeros(1,length(conditions));
eMove_counter = zeros(1,length(conditions));
blocking_counter = zeros(1,length(conditions));

for ec  = 1:1:length(erp.eventCodes)

    if erp.trialCodes(ec) > 0 %found a trial segment!
            
            fieldname = (char(conditions(erp.trialCodes(ec),:)));
                
            
            trialNumname = [fieldname,'_Num'];
            
            
            artifact = 0;
           
            
            %adjust the timewindow based on condition
            if (erp.trialCodes(ec) >=1 && erp.trialCodes(ec) <5)%Cue
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
            erp.cueCode.(fieldname)(counters(erp.trialCodes(ec)))=erp.cueCodes(ec);
            erp.picCode.(fieldname)(counters(erp.trialCodes(ec)))=erp.picCodes(ec);
            counters(erp.trialCodes(ec))= counters(erp.trialCodes(ec))+1;
              
        
    end;


end;

%now baseline it!

for condition = 1:1:length(conditions)

    fieldname = (char(conditions(condition,:)));
    
    baseline = repmat(squeeze(mean(erp.trial.(fieldname)(:,:,erp.baseline),3)),[1,1,size(erp.trial.(fieldname),3)]);     
    erp.trial.(fieldname)= erp.trial.(fieldname)-baseline;
    
    reject_baseline = repmat(squeeze(mean(erp.reject_trial.(fieldname)(:,:,erp.baseline),3)),[1,1,size(erp.trial.(fieldname),3)]);     
    erp.reject_trial.(fieldname)= erp.reject_trial.(fieldname)-reject_baseline;

end;

%Now create the erp&ave files
enc_erp = struct();
reject_enc_erp = struct();
enc_erp_ave = struct();

channels = {'VEOG';'HEOG';'F3';'F4';'C3';'C4';'P3';'P4';'PO3';'PO4';'O1';'O2';'OL';'OR';'T3';'T4';'T5';'T6';'Fz';'Cz';'Pz';'L_mas'};
%channels_r = {'F4';'F3';'C4';'C3';'P4';'P3';'PO4';'PO3';'O2';'O1';'OR';'OL';'T4';'T3';'T6';'T5';'Fz';'Cz';'Pz'};
for condition = 1:1:length(conditions)
    condition_name = (char(conditions(condition,:)));
    %condition_name_l = (char(conditions_LR(condition,:)));
    %condition_name_r = (char(conditions_LR(condition+8,:)));
    trial_name = [condition_name,'_trial'];
    cue_name = [condition_name,'_cue'];
    pic_name = [condition_name,'_pic'];
    enc_erp.(trial_name)=[squeeze(erp.trialNum.(condition_name))];
    enc_erp.(cue_name)=[squeeze(erp.cueCode.(condition_name))];
    enc_erp.(pic_name)=[squeeze(erp.picCode.(condition_name))];
    
    reject_enc_erp.(trial_name)=[squeeze(erp.trialNum.(condition_name))];
    reject_enc_erp.(cue_name)=[squeeze(erp.cueCode.(condition_name))];
    reject_enc_erp.(pic_name)=[squeeze(erp.picCode.(condition_name))];
    
    for channel = 1:1:length(channels)
        channel_name = char(channels(channel,:));
        %channel_name_r = char(channels_r(channel,:));
        fieldname = [condition_name,'_',channel_name];
        fieldname_hilb = [condition_name,'_',channel_name,'_hilb'];
        %fieldname_l = [condition_name_l,'_',channel_name];
        %fieldname_r = [condition_name_r,'_',channel_name_r];
        

        enc_erp.(fieldname) = [squeeze(erp.trial.(condition_name)(:,channel,:))];
        enc_erp.(fieldname_hilb) = [squeeze(erp.trial_hilb.(condition_name)(:,channel,:))];
        reject_enc_erp.(fieldname) = [squeeze(erp.reject_trial.(condition_name)(:,channel,:))];
        
        enc_erp_ave.(fieldname) = nanmean(enc_erp.(fieldname)(:,:),1);


    end;
    
    
    
end;

%summarize artifact
for condition = 1:1:length(conditions)
    condition_name = char(conditions(condition,:));
    art_free_name = [condition_name,'_art_free'];
    blink_name = [condition_name,'_blink'];
    eMove_name = [condition_name,'_eMove'];
    blocking_name = [condition_name,'_blocking'];
    enc_art.(art_free_name)= art_free_counter(condition)/(art_free_counter(condition)+art_counter(condition));
    enc_art.(blink_name)= blink_counter(condition);
    enc_art.(eMove_name)= eMove_counter(condition);
    enc_art.(blocking_name)= blocking_counter(condition);
    
end;


%save('erp.mat','erp');
save('enc_erp.mat','enc_erp')
save('reject_enc_erp.mat','reject_enc_erp')
save('enc_erp_ave.mat','enc_erp_ave');
save('enc_art.mat','enc_art');
% 

 