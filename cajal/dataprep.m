%data prep!

clear erp;
dir = cd;
%sub_dir = char(dir(length(dir)-1:length(dir)));
%datafile = [dir,'/',sub_dir,'_monkey_human_pd.rdf'];

ninputchan = 8

%convert raw file!
erp.data = Cajal_data.eeg(1:ninputchan,:);
%erp.data = erp.data*0.00819;



% Build artifact rejection file
erp.arf = build_arf_wood;

% Artifact Rejection
erp = arf_wood(erp);
L_mas_data = erp.data(ninputchan,:);
L_mas_data = repmat(L_mas_data,ninputchan,1);
erp.ave_ref_data = erp.data-(L_mas_data)/2;
%erp.ave_ref_data = erp.data;
erp.filtered_data = eegfilt(erp.ave_ref_data,250,0,30);
disp('filt complete');
save('erp.mat','erp')
% Normalize waveforms
%erp = normERPs(erp);

%
do_erp;

