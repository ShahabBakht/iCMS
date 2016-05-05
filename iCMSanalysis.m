% To add:
% 1- modifying the cut to conditions function
%   1-1 Separating different stimulation frequencies
%   1-2 Making the function work with wideband,lfp, and spikes data
% 2- Are the sorted spikes really spikes?

classdef iCMSanalysis < handle
    properties
        SamplingRate = 30000;
        channels;
    end
    
    methods
        function self = iCMSanalysis(Channels)
            self.channels = Channels;
        end
        function [OnStimData, PreData, PostData] = getDataChunks(self,wb,stim_times,delta_pre,delta_post,t0)
            SR = self.SamplingRate;
            NumEvents = length(stim_times);
            
            for ecount = 1:NumEvents
                StimStartTime = stim_times(ecount).rpStartTime;
                StimEndTime = stim_times(ecount).rpEndTime;
                data_temp = wb.getDataByTime(self.channels, StimStartTime,StimEndTime);
                if size(data_temp,2) > stim_times(ecount).train_dur * SR/1000
                    data_temp = data_temp(:,1:end-1);
                end
                
                OnStimData(ecount,:,:) = data_temp;
                
                PreEndTime = StimStartTime - t0;
                PreStartTime = PreEndTime - delta_pre;
                data_temp = wb.getDataByTime(self.channels, PreStartTime,PreEndTime-1/SR);
                PreData(ecount,:,:) = data_temp;
                
                PostStartTime = StimEndTime + t0;
                PostEndTime = PostStartTime + delta_post;
                data_temp = wb.getDataByTime(self.channels, PostStartTime,PostEndTime-1/SR);
                PostData(ecount,:,:) = data_temp;
                
            end
            
        end
        
        function LFP = getLFP(self,wb)
            % This extracts LFP signals from N channels of wideband recordings
            % The downsampling rate, low-pass and high-pass edges of the
            % band-pass filter are hardcoded.
            
            scale        = 30;      % this scale changes the sampling rate to 1k
            low          = 0.1;     % Hz
            high         = 200;     % Hz
            
            
            sr = self.SamplingRate / scale;
            Ch = self.channels;
            
            data_temp = wb.getDataByTime(Ch);
            neur_ds = Resample(self,data_temp,scale);
            clear data_temp
            LFP = bpfilter(self,neur_ds',low,high,sr);
            
            fprintf('getting LFP --> ')
            fprintf('\n')
       
        end
        
        
        
        function signal_ds = Resample(self,neur,scale)
            % neur: NxT matrix --> N number of channels, and T the time
            % length of the signals.
            fprintf('Resampling --> ')
            signal_ds = resample(neur',1,scale);
        end
        
        function signal_bp = bpfilter(self,neur,low,high,sr)
            % neur: NxT matrix --> N number of channels, and T the time
            % length of the signals.
            %signal_bp=bpfilter(self,detrend(ne(:,:))',20,30);
            fprintf('Band-pass filtering --> ')
            nyq = sr/2;
          
            %NEu=zeros(size(neur));
            % [z,p,k]=ellip(6,0.5,20,[low/nyq high/nyq]);
            [z,p,k]=butter(9,[low/nyq high/nyq]);
            [sos,g]=zp2sos(z,p,k);
            H=dfilt.df2sos(sos,g);
            
            
            sig=ones(size(neur,1),size(neur,2)*2).*repmat(mean(neur,2),1,size(neur,2)*2);
            sig(:,(ceil(size(neur,2)-size(neur,2)/2)):(ceil(size(neur,2)-size(neur,2)/2))+size(neur,2)-1)=(neur);
            
            for ii = 1:size(neur,1)
                sig1=filter(H,sig(ii,:));
                sig1=fliplr(sig1);
                sig1=filter(H,sig1);
                sig1=fliplr(sig1);
                NEu(ii,:)=sig1(1,(ceil(size(neur,2)-size(neur,2)/2)):(ceil(size(neur,2)-size(neur,2)/2))+size(neur,2)-1);
%                 NEu(ii,:) = filtfilt(sos,g,sig(ii,:));
            end
            
            signal_bp = NEu;
        end
    end
end



