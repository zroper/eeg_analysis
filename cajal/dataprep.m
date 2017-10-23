%data prep!
load erp.mat
erp.data = erp.data(any(erp.data,2),:); %removes all zeroed-out channels

%clear erp;
dir = cd;
%sub_dir = char(dir(length(dir)-1:length(dir)));
%datafile = [dir,'/',sub_dir,'_monkey_human_pd.rdf'];

ninputchan = 18;

%convert raw file!
%erp.data = Cajal_data.eeg(1:ninputchan,:);
%erp.data = erp.data*0.00819;
% erp.srate = Cajal_data.eeg8.streams.EEG8.fs;
% erp.eventCodes = Cajal_data.evnt.scalars.EVNT.data;
% erp.eventTimes = Cajal_data.evnt.scalars.EVNT.ts*1000;
% 
% erp.data=cat(1,Cajal_data.eeg8.streams.EEG8.data,Cajal_data.eeg16.streams.EEG6.data);
% erp.pnts = length(erp.data);
% erp.rateAcq = 1/erp.srate;

% Build artifact rejection file
erp.arf = build_arf_wood;

% Artifact Rejection
erp = arf_wood(erp);
%L_mas_data = erp.data(ninputchan,:);
%L_mas_data = repmat(L_mas_data,ninputchan,1);
erp.ave_ref_data = double(erp.data);%-(L_mas_data)/2;
%erp.ave_ref_data = erp.data;
erp.filtered_data = eegfilt(erp.ave_ref_data,erp.srate,0,100);
disp('filt complete');
save('erp.mat','erp')
% Normalize waveforms
%erp = normERPs(erp);

%
do_erp;

