function arf = build_arf_wood

% Define channels
arf.chans = 1:22;
arf.chanLabels = {'VEM','HEM','F3','F4','C3','C4','P3','P4','PO3','PO4','O1','O2','OL','OR','T3','T4','T5','T6','Fz','Cz','Pz','L_mas'};

% Define artifact rejection criteria
block = 1;
rHEM = 30; % horizontal eye movement
rVEM = 100; % vertical eye movement

for c = 1:1:length(arf.chans)
    if c>2
        thresh(c) = block;
        criter{c} = 'block';
    elseif c == 1
        thresh(c) = rVEM;
        criter{c} = 'blink';
    elseif c == 2
        thresh(c) = rHEM;
        criter{c} = 'eyemove';
    end
    arf.thresh(c) = thresh(c);
    arf.criter{c} = criter{c};
end

        