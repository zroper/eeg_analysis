function erp = arf_wood(erp)
% function for rejecting trials

for i = 1:size(erp.data,1) % go through each channel and apply correct arf criteria
%     if i == 1 % VEM channel
%         % Check for blinks
%         erp.arf.blink = ppa_wood(erp,i);
%     end;
%     
% %     if i == 1  % HEM channel
% %         % Check for eye movements
% %         erp.arf.veMove = step_wood(erp,i);
% %     end;
%         
%     if i == 2  % HEM channel
%         % Check for eye movements
%         erp.arf.eMove = step_wood(erp,i);
%     end;
    
    
    if i >0 && i < 9
        if i ~=2
            if i~=4
                erp.arf.blocking(i,:) = blocking(erp,i);
            end
        end;
    end
end
