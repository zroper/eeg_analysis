%attempts to bin trials

load erp;

conditions_Search = {'T_Left_1';'T_Left_2';'T_Left_3';'T_Right_1';'T_Right_2';'T_Right_3';'T_up';'T_down'};
channels = {'OZ1';'blank1';'PO4';'blank2';'PO3';'OR';'OL';'OZ2'};

%baseline_window = [1:50];%instep
%target_window = [51:150];

erp.baseline = [1:200];
erp.mem_window = [200,600];


erp.trialCodes = zeros(1,length(erp.eventCodes));
erp.trialNumCodes = zeros(1,length(erp.eventCodes));
trial_counter = 0;
%GOTTA SORT OUT THE CONDITION ORDER!!!!!!!!
%IMPLEMENT TRIALCODES based on eventcodes!!!!!!!!!
for ec = 1:1:length(erp.eventCodes)-65
    if erp.eventCodes(ec) == 2651 % Target Onset
        search_start = erp.eventTimes(ec+1);
        search_end = search_start+1000;
        searching = 1;
        counter = 1;
        while searching ==1
            current_time = erp.eventTimes(ec+counter);
            if current_time < search_end && erp.eventCodes(ec+counter) == 2600 %Correct Trial
                counter2 = counter;
                while searching ==1
                    if current_time < search_end && erp.eventCodes(ec+counter2) == 5045 %Right Target
                        erp.trialCodes(ec) = 4;
                        searching =0;
                    elseif current_time < search_end && erp.eventCodes(ec+counter2) == 5090 %Right Target
                        erp.trialCodes(ec) = 5;
                        searching =0;
                    elseif current_time < search_end && erp.eventCodes(ec+counter2) == 5135 %Right Target
                        erp.trialCodes(ec) = 6;
                        searching =0;
                    elseif current_time < search_end && erp.eventCodes(ec+counter2) == 5225 %Left Target
                        erp.trialCodes(ec) = 1;
                        searching =0;
                    elseif current_time < search_end && erp.eventCodes(ec+counter2) == 5270 %Left Target
                        erp.trialCodes(ec) = 2;
                        searching =0;
                    elseif current_time < search_end && erp.eventCodes(ec+counter2) == 5315 %Left Target
                        erp.trialCodes(ec) = 3;
                        searching =0;
                    elseif current_time < search_end && erp.eventCodes(ec+counter2) == 5000 %Up Target
                        erp.trialCodes(ec) = 7;
                        searching =0;
                    elseif current_time < search_end && erp.eventCodes(ec+counter2) == 5180 %Down Target
                        erp.trialCodes(ec) = 8;
                        searching =0;                        
                        
                    elseif current_time > search_end
                        searching = 0;
                    end;
                    counter2 = counter2 + 1;
                end
            elseif current_time > search_end
                searching = 0;
            end
            counter = counter + 1;
        end
    end
end

    
%IMPLEMENT TRIALCODES based on eventcodes!!!!!!!!! 

counters = ones(1,length(conditions_Search));
art_free_counter = zeros(1,length(conditions_Search));
art_counter = zeros(1,length(conditions_Search));
blink_counter = zeros(1,length(conditions_Search));
eMove_counter = zeros(1,length(conditions_Search));
blocking_counter = zeros(1,length(conditions_Search));

for ec  = 1:1:length(erp.eventCodes)

    if erp.trialCodes(ec) > 0 %found a trial segment!
            
            fieldname = (char(conditions_Search(erp.trialCodes(ec),:)));
                
            
 %           trialNumname = [fieldname,'_Num'];
            
            
            artifact = 0;
           
            
            %adjust the timewindow based on condition
            if (erp.trialCodes(ec) >=1 && erp.trialCodes(ec) <9)%Cue
                pre_timepoint = erp.mem_window(1);
                post_timepoint = erp.mem_window(2);
                %pre_timepoint_hilb = erp.mem_window(1)+300;%to give ample time for filters!
                %post_timepoint_hilb = erp.mem_window(2)+300;%to give ample time for filters!
            end;
                
            
            %check the time range if it's artifact free
%             
%             if sum(erp.arf.blink(erp.eventTimes(ec)-pre_timepoint:erp.eventTimes(ec)+post_timepoint))>0
%                 artifact = 1;
%                 blink_counter(erp.trialCodes(ec))=blink_counter(erp.trialCodes(ec))+1; 
%             elseif sum(erp.arf.eMove(erp.eventTimes(ec)-pre_timepoint:erp.eventTimes(ec)+post_timepoint))>0
%                 artifact = 1;
%                 eMove_counter(erp.trialCodes(ec))=eMove_counter(erp.trialCodes(ec))+1; 
%             elseif sum(sum(erp.arf.blocking(:,erp.eventTimes(ec)-pre_timepoint:erp.eventTimes(ec)+post_timepoint)))>0
%                 artifact = 1;
%                 blocking_counter(erp.trialCodes(ec))=blocking_counter(erp.trialCodes(ec))+1; 
%             end

            if sum(sum(erp.arf.blocking(:,erp.eventTimes(ec)-pre_timepoint:erp.eventTimes(ec)+post_timepoint)))>0
                artifact = 1;
                blocking_counter(erp.trialCodes(ec))=blocking_counter(erp.trialCodes(ec))+1; 
            end
            
             if artifact == 0
                 art_free_counter(erp.trialCodes(ec))=art_free_counter(erp.trialCodes(ec))+1; 
                 erp.trial.(fieldname)(counters(erp.trialCodes(ec)),:,:)= erp.filtered_data(:,erp.eventTimes(ec)-pre_timepoint:erp.eventTimes(ec)+post_timepoint);
                %erp.trial_hilb.(fieldname)(counters(erp.trialCodes(ec)),:,:)= erp.filtered_data(:,erp.eventTimes(ec)-pre_timepoint_hilb:erp.eventTimes(ec)+post_timepoint_hilb);
                 erp.reject_trial.(fieldname)(counters(erp.trialCodes(ec)),:,:)= nan(8,pre_timepoint+post_timepoint+1);
                %erp.reject_trial_hilb.(fieldname)(counters(erp.trialCodes(ec)),:,:)= nan(22,pre_timepoint_hilb+post_timepoint_hilb+1);
             elseif artifact >0
                 art_counter(erp.trialCodes(ec))=art_counter(erp.trialCodes(ec))+1; 
                erp.trial.(fieldname)(counters(erp.trialCodes(ec)),:,:)= nan(8,pre_timepoint+post_timepoint+1);
                %erp.trial_hilb.(fieldname)(counters(erp.trialCodes(ec)),:,:)= nan(22,pre_timepoint_hilb+post_timepoint_hilb+1);
                
                erp.reject_trial.(fieldname)(counters(erp.trialCodes(ec)),:,:)= erp.filtered_data(8,erp.eventTimes(ec)-pre_timepoint:erp.eventTimes(ec)+post_timepoint);
                %erp.reject_trial_hilb.(fieldname)(counters(erp.trialCodes(ec)),:,:)= erp.filtered_data(:,erp.eventTimes(ec)-pre_timepoint_hilb:erp.eventTimes(ec)+post_timepoint_hilb);
             end;
            %erp.trialNum.(fieldname)(counters(erp.trialCodes(ec)))=erp.trialNumCodes(ec);
            counters(erp.trialCodes(ec))= counters(erp.trialCodes(ec))+1;
              
        
    end;


end;

%now baseline it!

for condition = 1:1:length(conditions_Search)

    fieldname = (char(conditions_Search(condition,:)));
    
    baseline = repmat(squeeze(mean(erp.trial.(fieldname)(:,:,erp.baseline),3)),[1,1,size(erp.trial.(fieldname),3)]);     
    erp.trial.(fieldname)= erp.trial.(fieldname)-baseline;
    
    reject_baseline = repmat(squeeze(mean(erp.reject_trial.(fieldname)(:,:,erp.baseline),3)),[1,1,size(erp.trial.(fieldname),3)]);     
    erp.reject_trial.(fieldname)= erp.reject_trial.(fieldname)-reject_baseline;

end;

%Now create the erp&ave files
Cajal_search_erp = struct();
reject_Cajal_search_erp = struct();
Cajal_search_erp_ave = struct();

channels = {'OZ1';'blank1';'PO4';'blank2';'PO3';'OR';'OL';'OZ2'};
channels_r = {'OZ1';'blank1';'PO4';'blank2';'PO3';'OR';'OL';'OZ2'};
channel_r_index = [1 2 4 3 6 5 8];
for condition = 1:1:length(conditions_Search)
    condition_name = (char(conditions_Search(condition,:)));
    condition_name = [condition_name];


    %trial_name = [condition_name,'_trial'];
    %choice1_erp.(trial_name)=[squeeze(erp.trialNum.(condition_name))];

    
    %reject_choice1_erp.(trial_name)=[squeeze(erp.trialNum.(condition_name))];

    
    for channel = 1:1:length(channels)
        channel_name = char(channels(channel,:));
        %channel_name_r = char(channels_r(channel,:));
        fieldname = [condition_name,'_',channel_name];
        
        Cajal_search_erp.(fieldname) = [squeeze(erp.trial.(condition_name)(:,channel,:));];
        %hm_pd_erp.(lat_fieldname) = [squeeze(erp.trial.(l_condition_name)(:,channel,:));squeeze(erp.trial.(r_condition_name)(:,channel_r_index(channel),:))];
        %hm_pd_erp.(mid_fieldname) = [squeeze(erp.trial.(condition_name)(:,channel,:));squeeze(erp.trial.(b_condition_name)(:,channel,:))];
        %hm_pd_erp.(lat_fieldname_hilb) = [squeeze(erp.trial_hilb.(l_condition_name)(:,channel,:));squeeze(erp.trial_hilb.(r_condition_name)(:,channel_r_index(channel),:))];
        %hm_pd_erp.(mid_fieldname_hilb) = [squeeze(erp.trial_hilb.(condition_name)(:,channel,:));squeeze(erp.trial_hilb.(b_condition_name)(:,channel,:))];
        
        reject_Cajal_search_erp.(fieldname) = [squeeze(erp.reject_trial.(condition_name)(:,channel,:));];
        %reject_hm_pd_erp.(lat_fieldname) = [squeeze(erp.reject_trial.(l_condition_name)(:,channel,:));squeeze(erp.reject_trial.(r_condition_name)(:,channel_r_index(channel),:))];
        %reject_hm_pd_erp.(mid_fieldname) = [squeeze(erp.reject_trial.(condition_name)(:,channel,:));squeeze(erp.reject_trial.(b_condition_name)(:,channel,:))];
        
        
        Cajal_search_erp_ave.(fieldname) = nanmean(Cajal_search_erp.(fieldname)(:,:),1);
        %hm_pd_erp_ave.(lat_fieldname) = nanmean(hm_pd_erp.(lat_fieldname)(:,:),1);
        %hm_pd_erp_ave.(mid_fieldname) = nanmean(hm_pd_erp.(mid_fieldname)(:,:),1);


    end;
    
    
    
end;



%summarize artifact
for condition = 1:1:length(conditions_Search)
    condition_name = (char(conditions_Search(condition,:)));
%     lat_condition_name = ['Lat_',condition_name];
%     mid_condition_name = ['Mid_',condition_name];
%     
    art_free_name = [condition_name,'_art_free'];
    blink_name = [condition_name,'_blink'];
    eMove_name = [condition_name,'_eMove'];
    blocking_name = [condition_name,'_blocking'];
    
%     lat_art_free_name = [lat_condition_name,'_art_free'];
%     lat_blink_name = [lat_condition_name,'_blink'];
%     lat_eMove_name = [lat_condition_name,'_eMove'];
%     lat_blocking_name = [lat_condition_name,'_blocking'];
%     
%     mid_art_free_name = [mid_condition_name,'_art_free'];
%     mid_blink_name = [mid_condition_name,'_blink'];
%     mid_eMove_name = [mid_condition_name,'_eMove'];
%     mid_blocking_name = [mid_condition_name,'_blocking'];

    Cajal_search_art.(art_free_name) = art_free_counter(condition)/(art_free_counter(condition)+art_counter(condition));
    Cajal_search_art.(blink_name) = blink_counter(condition);
    Cajal_search_art.(eMove_name) = eMove_counter(condition);
    Cajal_search_art.(blocking_name) = blocking_counter(condition);
    
%     choice1_art.(lat_art_free_name)= (art_free_counter((condition-1)*4+1)+art_free_counter((condition-1)*4+4))/(art_free_counter((condition-1)*4+1)+art_free_counter((condition-1)*4+4)+art_counter((condition-1)*4+1)+art_counter((condition-1)*4+4));
%     choice1_art.(lat_blink_name)= (blink_counter((condition-1)*4+1)+blink_counter((condition-1)*4+4));
%     choice1_art.(lat_eMove_name)= (eMove_counter((condition-1)*4+1)+eMove_counter((condition-1)*4+4));
%     choice1_art.(lat_blocking_name)= (blocking_counter((condition-1)*4+1)+blocking_counter((condition-1)*4+4));
%     
%     choice1_art.(mid_art_free_name)= (art_free_counter((condition-1)*4+2)+art_free_counter((condition-1)*4+3))/(art_free_counter((condition-1)*4+2)+art_free_counter((condition-1)*4+3)+art_counter((condition-1)*4+2)+art_counter((condition-1)*4+3));
%     choice1_art.(mid_blink_name)= (blink_counter((condition-1)*4+2)+blink_counter((condition-1)*4+3));
%     choice1_art.(mid_eMove_name)= (eMove_counter((condition-1)*4+2)+eMove_counter((condition-1)*4+3));
%     choice1_art.(mid_blocking_name)= (blocking_counter((condition-1)*4+2)+blocking_counter((condition-1)*4+3));
    
end;


%save('erp.mat','erp');
save('Cajal_Search_ERP.mat','Cajal_search_erp')
%save('reject_enc_erp.mat','reject_enc_erp')
save('Cajal_Search_ERP_ave.mat','Cajal_search_erp_ave');
save('Cajal_Search_ERP_art.mat','Cajal_search_art');


 