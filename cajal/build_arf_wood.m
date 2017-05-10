function arf = build_arf_wood

% Define channels
arf.chans = 1:24;
%arf.chanLabels = {'OZ1';'blank1';'PO4';'blank2';'PO3';'OR';'OL';'OZ2'};
arf.chanLabels = {'chan1';'chan2';'chan3';'chan4';'chan5';'chan6';...
                'chan7';'chan8';'chan9';'chan10';'chan11';'chan12';...
                'chan13';'chan14';'chan15';'chan16';'chan17';'chan18';...
                'chan19';'chan20';'chan21';'chan22';'chan23';'chan24'};
% Define artifact rejection criteria
block = 1;
rHEM = 1000;%30; % horizontal eye movement
rVEM = 1000;%100; % vertical eye movement

for c = 1:1:length(arf.chans)
    if c>0
        thresh(c) = block;
        criter{c} = 'block';
%     elseif c == 1
%         thresh(c) = rVEM;
%         criter{c} = 'blink';
%     elseif c == 2
%         thresh(c) = rHEM;
%         criter{c} = 'eyemove';
    end
    arf.thresh(c) = thresh(c);
    arf.criter{c} = criter{c};
end

        