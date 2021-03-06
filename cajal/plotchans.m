conditions_Search = {'T_Left_1';'T_Left_2';'T_Left_3';'T_Right_1';'T_Right_2';'T_Right_3';'T_up';'T_down'};
channels = {'chan1';'chan2';'chan3';'chan4';'chan5';'chan6';...
                'chan7';'chan8';'chan9';'chan10';'chan11';'chan12';...
                'chan13';'chan14';'chan15';'chan16';'chan17';'chan18'};

figure
hold
for chan = 1:length(channels)
        clear datacond
        datasum = zeros(1,801);
        for pos = 4:6
            datacond=eval(['Cajal_search_erp_ave.',conditions_Search{pos},'_',channels{chan}]);
            datasum = datasum+datacond;
        end
        data = datasum;
        plot(data)
end

% clear datacond
% datasum = zeros(1,801);
% for chan = 1:length(channels)
%     for pos = 4:6
%         datacond=eval(['Cajal_search_erp_ave.',conditions_Search{pos},'_',channels{chan}]);
%         datasum = datasum+datacond;
%     end
% end
% 
% data = datasum;
% plot(data)