% To do:
% 1- Add the filtering
% 2- Separating different stimulation frequencies

classdef iCMSanalysis < handle
    properties
        SamplingRate = 30000;
    end
    
    methods
        function self = iCMSanalysis()
            
        end
        function [OnStimData, PreData, PostData] = getDataChunks(self,wb,stim_times,delta_pre,delta_post,t0)
            SR = self.SamplingRate;
            NumEvents = length(stim_times);
            channel = 1:2;
            for ecount = 1:NumEvents
                StimStartTime = stim_times(ecount).rpStartTime;
                StimEndTime = stim_times(ecount).rpEndTime;
                data_temp = wb.getDataByTime(channel, StimStartTime,StimEndTime);
                if size(data_temp,2) > stim_times(ecount).train_dur * SR/1000
                    data_temp = data_temp(:,1:end-1);
                end
                
                OnStimData(ecount,:,:) = data_temp;
                
                PreEndTime = StimStartTime - self.t0;
                PreStartTime = PreEndTime - delta_pre;
                data_temp = wb.getDataByTime(channel, PreStartTime,PreEndTime-1/SR);
                PreData(ecount,:,:) = data_temp;
                
                PostStartTime = StimEndTime + t0;
                PostEndTime = PostStartTime + delta_post;
                data_temp = wb.getDataByTime(channel, PostStartTime,PostEndTime-1/SR);
                PostData(ecount,:,:) = data_temp;
                
            end
            
        end
        function signal_bp = bpfilter(self,signal,L,H)
            SR = self.SamplingRate;
            
            
        end
    end
end



