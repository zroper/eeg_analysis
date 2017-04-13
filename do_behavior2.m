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

%first figure out the modification for each shape
reward_shape = reward_ext_data2.rewardmodification(1);
baseline_shape = reward_ext_data2.rewardmodification(2);
punishment_shape = reward_ext_data2.rewardmodification(3);

reward_color = reward_ext_data1.TargetColorIndex(reward_shape);
baseline_color = reward_ext_data1.TargetColorIndex(baseline_shape);
punishment_color = reward_ext_data1.TargetColorIndex(punishment_shape);
no_reward_color = reward_ext_data1.TargetColorIndex(4);


%do for reward_ext1!
choice_RT1 = nan(8,96,3,5);%block x reward_manip x selected_color x trial
for t = 1:1:length(reward_ext_data1.RT)
    block = reward_ext_data1.block(t);
    trial = reward_ext_data1.trial(t);
    shape = reward_ext_data1.shape(t);
    rewardloc = reward_ext_data1.rewardloc(t);
    response = reward_ext_data1.response(t);
    distlocs = reward_ext_data1.distlocs(t,:);
    distcols = reward_ext_data1.distcols(t,:);
    RT = reward_ext_data1.RT(t);
    if response == 0
        response_color = 5;
    elseif response > 0 && response < 5
        response_color = reward_ext_data1.TargetColorIndex(shape);
    elseif response > 100
        response = response-100;
        response_color = reward_ext_data1.TargetColorIndex(distcols(find(distlocs==response)));
    end;
    
    
    %figure out responded color!
    
    choice_RT1(block,trial,shape,response_color)=RT;

end;

%Now squeeeeze!
for block = 1:1:8
    for reward = 1:1:3
        for response = 1:1:5
            relevant_index = [];
            relevant_index = find(isnan(squeeze(choice_RT1(block,:,reward,response)))==0);
            if length(relevant_index)>0
                choice_RT_data1.count(block,reward,response) = length(relevant_index);
                choice_RT_data1.prop(block,reward,response) = length(relevant_index)/32;
                choice_RT_data1.meanRTs(block,reward,response) = mean(choice_RT1(block,relevant_index,reward,response));
                choice_RT_data1.medianRTs(block,reward,response) = median(choice_RT1(block,relevant_index,reward,response));
            else
                choice_RT_data1.count(block,reward,response) = 0;
                choice_RT_data1.prop(block,reward,response) = 0;
                choice_RT_data1.meanRTs(block,reward,response) = nan;
                choice_RT_data1.medianRTs(block,reward,response) = nan;
            end;
        end;
            
            
        
    end;
end;

%Do summary too...



save('choice_RT_data1.mat','choice_RT_data1');

%do for reward_ext2!
choice_RT2 = nan(4,96,3,5);%block x reward_manip x selected_color x trial
for t = 1:1:length(reward_ext_data2.RT)
    block = reward_ext_data2.block(t);
    trial = reward_ext_data2.trial(t);
    shape = reward_ext_data2.shape(t);
    rewardloc = reward_ext_data2.rewardloc(t);
    response = reward_ext_data2.response(t);
    distlocs = reward_ext_data2.distlocs(t,:);
    distcols = reward_ext_data2.distcols(t,:);
    RT = reward_ext_data2.RT(t);
    if response == 0
        response_color = 5;
    elseif response > 0 && response < 5
        response_color = reward_ext_data2.TargetColorIndex(shape);
    elseif response > 100
        response = response-100;
        response_color = reward_ext_data1.TargetColorIndex(distcols(find(distlocs==response)));
    end;
    
    
    %figure out responded color!
    
    choice_RT2(block,trial,shape,response_color)=RT;

end;

%Now squeeeeze!
for block = 1:1:4
    for reward = 1:1:3
        for response = 1:1:5
            relevant_index = [];
            relevant_index = find(isnan(squeeze(choice_RT2(block,:,reward,response)))==0);
            if length(relevant_index)>0
                choice_RT_data2.count(block,reward,response) = length(relevant_index);
                choice_RT_data2.prop(block,reward,response) = length(relevant_index)/32;
                choice_RT_data2.meanRTs(block,reward,response) = mean(choice_RT2(block,relevant_index,reward,response));
                choice_RT_data2.medianRTs(block,reward,response) = median(choice_RT2(block,relevant_index,reward,response));
            else
                choice_RT_data2.count(block,reward,response) = 0;
                choice_RT_data2.prop(block,reward,response) = 0;
                choice_RT_data2.meanRTs(block,reward,response) = nan;
                choice_RT_data2.medianRTs(block,reward,response) = nan;
            end;
        end;
            
            
        
    end;
end;

%Do summary too...

save('choice_RT_data2.mat','choice_RT_data2');




% %LET"S DO SEARCHES!!!
search_RT1_counter = zeros(3,4,5);
search_RT1 = nan(3,4,5,300);%shape x targloc(as of color) x response
for t = 1:1:length(reward_ext_search_data1.trial)
    block = reward_ext_search_data1.block(t);
    trial = reward_ext_search_data1.trial(t);
    shape = reward_ext_search_data1.shape(t);
    rewardloc = reward_ext_search_data1.rewardloc(t);
    distlocs = reward_ext_search_data1.distlocs(t,:);
    distcols = reward_ext_search_data1.distcols(t,:);
    targetloc = reward_ext_search_data1.distcols(t);
    response = reward_ext_search_data1.response(t);
    RT = reward_ext_search_data1.RT(t);
    
    %let's figure out what color was selected!
    
    if response >0 && response < 100 %correct
        selected_color = response;
    elseif
    
    search_RT1(block,trial,shape,targetloc,rewardloc,1) = 
    
end;
