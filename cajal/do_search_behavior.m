cur_dir = cd;
sub = char(cur_dir(length(cur_dir)-1:length(cur_dir)));

filename1 = ['reward_ext_data1_',sub];
filename2 = ['reward_ext_data2_',sub];
filename3 = ['reward_ext_search_data1_',sub];
filename4 = ['reward_ext_search_data2_',sub];

load(filename1)
load(filename2)
load(filename3)
load(filename4)

% %LET"S DO SEARCHES!!!
TargetColorIndex = reward_ext_search_data1.TargetColorIndex;
search_RT1_counter = zeros(3,4,5);
search_RT1 = nan(3,4,5,300);%shape x targloc(as of color) x response
for t = 1:1:length(reward_ext_search_data1.trial)
    block = reward_ext_search_data1.block(t);
    trial = reward_ext_search_data1.trial(t);
    shape = reward_ext_search_data1.shape(t);
    rewardloc = reward_ext_search_data1.rewardloc(t);
    distlocs = reward_ext_search_data1.distlocs(t,:);
    distcols = reward_ext_search_data1.distcols(t,:);
    targetloc = reward_ext_search_data1.targetloc(t);
    response = reward_ext_search_data1.response(t);
    RT = reward_ext_search_data1.RT(t);
    
    %first reconstruct the color array
    clear color_array;
    for loc = 1:1:4
        if loc == rewardloc
            color_array(loc) = TargetColorIndex(shape);
        else
            dist_index = find(loc==distlocs);
            color_array(loc)=TargetColorIndex(distcols(dist_index));
        end;
    end;
    %figure out which color was target presented on
    target_color = color_array(targetloc);
    
    %let's figure out what color was selected!
    
    if response >0 && response < 100 %correct
        selected_color = color_array(response);
    elseif response == 0
        selected_color = 5;
    elseif response > 100
        selected_color = color_array(response-100);
    end;
    search_RT1_counter(shape,target_color,selected_color)=search_RT1_counter(shape,target_color,selected_color)+1;
    search_RT1(shape,target_color,selected_color,search_RT1_counter(shape,target_color,selected_color)) = RT;
    
end;


%Now congruent incongruent analysis
%Get congruent!
congruent_acc_count = 0;
congruent_acc_rts = [];
congruent_inc_count = 0;
congruent_inc_rts = [];

incongruent_acc_count = 0;
incongruent_acc_rts = [];
incongruent_inc_count = 0;
incongruent_inc_rts = [];

reward_capture_count = 0;
reward_capture_rts = [];
nonreward_capture_count = 0;
nonreward_capture_rts = [];



for shape = 1:1:3
    reward_color = TargetColorIndex(shape);
    nonreward_color = TargetColorIndex(TargetColorIndex ~=reward_color);
    
    for target_color = 1:1:4
        for selected_color = 1:1:4
            
            if target_color == reward_color
                if target_color == selected_color
                    congruent_acc_count = congruent_acc_count+search_RT1_counter(shape,target_color,selected_color);
                    congruent_acc_rts = [congruent_acc_rts, squeeze(search_RT1(shape,target_color,selected_color,:))'];
                elseif selected_color ~= target_color
                    congruent_inc_count = congruent_inc_count+search_RT1_counter(shape,target_color,selected_color);
                    congruent_inc_rts = [congruent_inc_rts, squeeze(search_RT1(shape,target_color,selected_color,:))'];
                end;
            elseif target_color ~= reward_color
                if target_color == selected_color
                    incongruent_acc_count = incongruent_acc_count+search_RT1_counter(shape,target_color,selected_color);
                    incongruent_acc_rts = [incongruent_acc_rts, squeeze(search_RT1(shape,target_color,selected_color,:))'];
                elseif target_color ~= selected_color
                    incongruent_inc_count = incongruent_inc_count+search_RT1_counter(shape,target_color,selected_color);
                    incongruent_inc_rts = [incongruent_inc_rts, squeeze(search_RT1(shape,target_color,selected_color,:))'];
                    
                    if selected_color == reward_color
                        reward_capture_count = reward_capture_count+search_RT1_counter(shape,target_color,selected_color);
                        reward_capture_rts = [reward_capture_rts, squeeze(search_RT1(shape,target_color,selected_color,:))'];
                    elseif selected_color ~= reward_color
                        nonreward_capture_count = nonreward_capture_count+search_RT1_counter(shape,target_color,selected_color);
                        nonreward_capture_rts = [nonreward_capture_rts, squeeze(search_RT1(shape,target_color,selected_color,:))'];
                    end;
                end;
            end;
            
        end;
        
    end;
    
end;

%Now summary!
search_RT_data1.congruent_acc_count = congruent_acc_count;
search_RT_data1.congruent_inc_count = congruent_inc_count;
search_RT_data1.congruent_acc = congruent_acc_count/(congruent_acc_count+congruent_inc_count);
search_RT_data1.congruent_acc_meanRT = nanmean(congruent_acc_rts);
search_RT_data1.congruent_acc_medianRT = nanmedian(congruent_acc_rts);
search_RT_data1.congruent_inc_meanRT = nanmean(congruent_inc_rts);
search_RT_data1.congruent_inc_medianRT = nanmedian(congruent_inc_rts);

search_RT_data1.incongruent_acc_count = incongruent_acc_count;
search_RT_data1.incongruent_inc_count = incongruent_inc_count;
search_RT_data1.incongruent_acc = incongruent_acc_count/(incongruent_acc_count+incongruent_inc_count);
search_RT_data1.incongruent_acc_meanRT = nanmean(incongruent_acc_rts);
search_RT_data1.incongruent_acc_medianRT = nanmedian(incongruent_acc_rts);
search_RT_data1.incongruent_inc_meanRT = nanmean(incongruent_inc_rts);
search_RT_data1.incongruent_inc_medianRT = nanmedian(incongruent_inc_rts);

search_RT_data1.reward_capture_count = reward_capture_count;
search_RT_data1.reward_capture_prop = reward_capture_count/(reward_capture_count+nonreward_capture_count);
search_RT_data1.reward_capture_meanRT = nanmean(reward_capture_rts);
search_RT_data1.reward_capture_medianRT = nanmedian(reward_capture_rts);
search_RT_data1.nonreward_capture_meanRT = nanmean(nonreward_capture_rts);
search_RT_data1.nonreward_capture_medianRT = nanmedian(nonreward_capture_rts);


save('search_RT_data1.mat','search_RT_data1');



%%
% %LET"S DO SEARCHES!!!
TargetColorIndex = reward_ext_search_data2.TargetColorIndex;
search_RT2_counter = zeros(3,4,5);
search_RT2 = nan(3,4,5,300);%shape x targloc(as of color) x response
for t = 1:1:length(reward_ext_search_data2.trial)
    block = reward_ext_search_data2.block(t);
    trial = reward_ext_search_data2.trial(t);
    shape = reward_ext_search_data2.shape(t);
    rewardloc = reward_ext_search_data2.rewardloc(t);
    distlocs = reward_ext_search_data2.distlocs(t,:);
    distcols = reward_ext_search_data2.distcols(t,:);
    targetloc = reward_ext_search_data2.targetloc(t);
    response = reward_ext_search_data2.response(t);
    RT = reward_ext_search_data2.RT(t);
    
    %first reconstruct the color array
    clear color_array;
    for loc = 1:1:4
        if loc == rewardloc
            color_array(loc) = TargetColorIndex(shape);
        else
            dist_index = find(loc==distlocs);
            color_array(loc)=TargetColorIndex(distcols(dist_index));
        end;
    end;
    %figure out which color was target presented on
    target_color = color_array(targetloc);
    
    %let's figure out what color was selected!
    
    if response >0 && response < 100 %correct
        selected_color = color_array(response);
    elseif response == 0
        selected_color = 5;
    elseif response > 100
        selected_color = color_array(response-100);
    end;
    search_RT2_counter(shape,target_color,selected_color)=search_RT2_counter(shape,target_color,selected_color)+1;
    search_RT2(shape,target_color,selected_color,search_RT2_counter(shape,target_color,selected_color)) = RT;
    
end;


%Now congruent incongruent analysis
%Get congruent!




for shape = 1:1:3
    
    congruent_acc_count = 0;
    congruent_acc_rts = [];
    congruent_inc_count = 0;
    congruent_inc_rts = [];

    incongruent_acc_count = 0;
    incongruent_acc_rts = [];
    incongruent_inc_count = 0;
    incongruent_inc_rts = [];

    reward_capture_count = 0;
    reward_capture_rts = [];
    nonreward_capture_count = 0;
    nonreward_capture_rts = [];

    reward_color = TargetColorIndex(shape);
    nonreward_color = TargetColorIndex(TargetColorIndex ~=reward_color);
    
    for target_color = 1:1:4
        for selected_color = 1:1:4
            
            if target_color == reward_color
                if target_color == selected_color
                    congruent_acc_count = congruent_acc_count+search_RT2_counter(shape,target_color,selected_color);
                    congruent_acc_rts = [congruent_acc_rts, squeeze(search_RT2(shape,target_color,selected_color,:))'];
                elseif selected_color ~= target_color
                    congruent_inc_count = congruent_inc_count+search_RT2_counter(shape,target_color,selected_color);
                    congruent_inc_rts = [congruent_inc_rts, squeeze(search_RT2(shape,target_color,selected_color,:))'];
                end;
            elseif target_color ~= reward_color
                if target_color == selected_color
                    incongruent_acc_count = incongruent_acc_count+search_RT2_counter(shape,target_color,selected_color);
                    incongruent_acc_rts = [incongruent_acc_rts, squeeze(search_RT2(shape,target_color,selected_color,:))'];
                elseif target_color ~= selected_color
                    incongruent_inc_count = incongruent_inc_count+search_RT2_counter(shape,target_color,selected_color);
                    incongruent_inc_rts = [incongruent_inc_rts, squeeze(search_RT2(shape,target_color,selected_color,:))'];
                    
                    if selected_color == reward_color
                        reward_capture_count = reward_capture_count+search_RT2_counter(shape,target_color,selected_color);
                        reward_capture_rts = [reward_capture_rts, squeeze(search_RT2(shape,target_color,selected_color,:))'];
                    elseif selected_color ~= reward_color
                        nonreward_capture_count = nonreward_capture_count+search_RT2_counter(shape,target_color,selected_color);
                        nonreward_capture_rts = [nonreward_capture_rts, squeeze(search_RT2(shape,target_color,selected_color,:))'];
                    end;
                end;
            end;
            
        end;
        
    end;
    
    search_RT_data2.congruent_acc_count(shape) = congruent_acc_count;
    search_RT_data2.congruent_inc_count(shape) = congruent_inc_count;
    search_RT_data2.congruent_acc(shape) = congruent_acc_count/(congruent_acc_count+congruent_inc_count);
    search_RT_data2.congruent_acc_meanRT(shape) = nanmean(congruent_acc_rts);
    search_RT_data2.congruent_acc_medianRT(shape) = nanmedian(congruent_acc_rts);
    search_RT_data2.congruent_inc_meanRT(shape) = nanmean(congruent_inc_rts);
    search_RT_data2.congruent_inc_medianRT(shape) = nanmedian(congruent_inc_rts);

    search_RT_data2.incongruent_acc_count(shape) = incongruent_acc_count;
    search_RT_data2.incongruent_inc_count(shape) = incongruent_inc_count;
    search_RT_data2.incongruent_acc(shape) = incongruent_acc_count/(incongruent_acc_count+incongruent_inc_count);
    search_RT_data2.incongruent_acc_meanRT(shape) = nanmean(incongruent_acc_rts);
    search_RT_data2.incongruent_acc_medianRT(shape) = nanmedian(incongruent_acc_rts);
    search_RT_data2.incongruent_inc_meanRT(shape) = nanmean(incongruent_inc_rts);
    search_RT_data2.incongruent_inc_medianRT(shape) = nanmedian(incongruent_inc_rts);

    search_RT_data2.reward_capture_count(shape) = reward_capture_count;
    search_RT_data2.reward_capture_prop(shape) = reward_capture_count/(reward_capture_count+nonreward_capture_count);
    search_RT_data2.reward_capture_meanRT(shape) = nanmean(reward_capture_rts);
    search_RT_data2.reward_capture_medianRT(shape) = nanmedian(reward_capture_rts);
    search_RT_data2.nonreward_capture_meanRT(shape) = nanmean(nonreward_capture_rts);
    search_RT_data2.nonreward_capture_medianRT(shape) = nanmedian(nonreward_capture_rts);

end;

%Now summary!



save('search_RT_data2.mat','search_RT_data2');


