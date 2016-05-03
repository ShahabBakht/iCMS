classdef WidebandManager < handle
    % Abstract Base Class for handling wideband data. 
    
    
    
    
    properties(GetAccess=public, SetAccess=protected)
        samplingRate; %Sampling rate, in Hz
        nSamples % Total number of samples per channel
        channels; % Per-Channel information
        
        
        refGroups; % Index of reference groups
    end
    
    properties(GetAccess=protected, SetAccess=protected)
        ;
    end
    
    methods(Abstract=true)                      
    	data = getDataByTime(self, chan_index, t_start, t_stop); 
        refGroup = getRefGroup(self, chan_index);   
        n = channelCount(self); %Total number of channels
    end            
end



