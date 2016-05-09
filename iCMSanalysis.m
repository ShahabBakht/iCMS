% To add/solve:
% 2- add analysis:
%   2-2 correlation
%   2-3 coherence
%   2-4 lfp spectrum

classdef iCMSanalysis < handle
    properties
        SamplingRate = 30000;
        channels;
        LFP
        lfpSamplingRate = 1000; 
        spikeTimes
        stimTimes
    end
    
    methods
        
        function self = iCMSanalysis(Channels,spikes_file,stim_file)
            self.channels = Channels;
            load(spikes_file);
            load(stim_file);
            self.spikeTimes = spiketimes;
            self.stimTimes = stim_times;
            
        end
        
        function [OnStimData, PreStimData, PostStimData] = getDataChunks(self,wb,type,delta_pre,delta_post,t0,goodUnitIdx_file)
            % This function gets any type of neural data (wideband, LFP,
            % spikes), and cut them to pre-, post-, and on-stimulation, as
            % well as different conditions (different stimulation
            % frequencies).
            % The input data will be loaded automatically based on the
            % 'type' input variable.
            % The output data structure: A{Fx1} cell array --> F is the number
            % of conditions (e.g. stimulation frequencies). Ai (the ith 
            % element of A) is a three-dimensional matrix. Ai(MxNxT) matrix -->
            % M: number of trials (repetitions of stimulation with same 
            % frequency), N: number of channels (or neurons in the case of spikes),
            % T: time length
            
            switch type
                case 'lfp'
                    SR = self.lfpSamplingRate;
                case 'wb' 
                    SR = self.SamplingRate;
                case 'spk'
                    SR = self.SamplingRate;
            end

            
            NumEvents = length(self.stimTimes);
            allStimFreqs = nan(NumEvents,1);
            for ecount = 1:NumEvents
                allStimFreqs(ecount) = self.stimTimes(ecount).train_freq;
            end
            uniqueStimFreqs = unique(allStimFreqs);
            NumConditions = length(uniqueStimFreqs); % number of conditions - now only the stimulation frequencies
            
            OnStimData = cell(NumConditions,1);
            PreStimData = cell(NumConditions,1);
            PostStimData = cell(NumConditions,1);
            counter = zeros(NumConditions,1);
            for ecount = 1:NumEvents
                fprintf(['ripple number ',num2str(ecount), ' \n'])
                
                % populating the on-stimulation cell array
                StimStartTime = self.stimTimes(ecount).rpStartTime;
                StimEndTime = self.stimTimes(ecount).rpEndTime;
                if strcmp(type,'lfp')
                    data_temp = getLFP(self,wb,StimStartTime,StimEndTime);
                elseif strcmp(type,'wb')
                    data_temp = wb.getDataByTime(self.channels, StimStartTime,StimEndTime);
                elseif strcmp(type,'spk')
                    data_temp = getSpikes(self,goodUnitIdx_file,StimStartTime,StimEndTime);
                end
                
                if ~strcmp(type,'spk')
                    if size(data_temp,2) > self.stimTimes(ecount).train_dur * SR/1000
                        data_temp = data_temp(:,1:end-1);
                    end
                end
                
                thisFreq = self.stimTimes(ecount).train_freq;
                whichCondition = thisFreq == uniqueStimFreqs;
                counter(whichCondition) = counter(whichCondition) + 1;
                if strcmp(type,'lfp') || strcmp(type,'wb') 
                    OnStimData{whichCondition}(counter(whichCondition),:,:) = data_temp;
                else
                    OnStimData{whichCondition}{counter(whichCondition),:,:} = data_temp;
                end

                
                
                %  populating the pre-stimulation cell array
                PreEndTime = StimStartTime - t0;
                PreStartTime = PreEndTime - delta_pre;
                if strcmp(type,'lfp')
                    data_temp = getLFP(self,wb,PreStartTime,PreEndTime-1/SR);
                elseif strcmp(type,'wb')
                    data_temp = wb.getDataByTime(self.channels, PreStartTime,PreEndTime-1/SR);
                elseif strcmp(type,'spk')
                    data_temp = getSpikes(self,goodUnitIdx_file,PreStartTime,PreEndTime-1/SR);
                end

                if strcmp(type,'lfp') || strcmp(type,'wb') 
                    PreStimData{whichCondition}(counter(whichCondition),:,:) = data_temp;
                else
                    PreStimData{whichCondition}{counter(whichCondition),:,:} = data_temp;
                end

                
                % populating the post-stimulation cell array
                PostStartTime = StimEndTime + t0;
                PostEndTime = PostStartTime + delta_post;
                if strcmp(type,'lfp')
                    data_temp = getLFP(self,wb,PostStartTime,PostEndTime-1/SR);
                elseif strcmp(type,'wb')
                    data_temp = wb.getDataByTime(self.channels, PostStartTime,PostEndTime-1/SR);
                elseif strcmp(type,'spk')
                    data_temp = getSpikes(self,goodUnitIdx_file,PostStartTime,PostEndTime-1/SR);
                end

                if strcmp(type,'lfp') || strcmp(type,'wb') 
                    PostStimData{whichCondition}(counter(whichCondition),:,:) = data_temp;
                else
                    PostStimData{whichCondition}{counter(whichCondition),:,:} = data_temp;
                end

                
            end
            
        end
        function goodSpikes = getGoodSpikes(self,goodUnitIdx_file)
            load(goodUnitIdx_file);
            goodSpikes = self.spikeTimes(goodunitsindex);
            
        end
        
        function SpikesinWindow = getSpikes(self,goodUnitIdx_file,startTime,endTime)
%             load(spikes_file);
            goodSpikes = getGoodSpikes(self,goodUnitIdx_file);
            spiketimes = goodSpikes;
            NumNeurons = length(spiketimes); % spiketimes is loaded from sorted spikes files
            SpikesinWindow = cell(size(spiketimes))';
            for neurcount = 1:NumNeurons
                
                thisNeuron = spiketimes{neurcount};
                thisNeuroninWindow = thisNeuron((thisNeuron >= startTime * 1000) & (thisNeuron < endTime * 1000));
                SpikesinWindow{neurcount} = thisNeuroninWindow;
                clear thisNeuroninWindow
            end
            
        end
        function LFPinWindow = getLFP(self,wb,startTime,endTime)
            % This function gives the LFP signal between the prespicified
            % temporal windows. If the start and end time are not
            % specified, it gives the whole data.
            
            if isempty(self.LFP)
            self.extractLFP(wb);
            end
            if startTime == 0 
                LFPinWindow = self.LFP(:,1:floor(endTime * self.lfpSamplingRate));
            elseif startTime ~= 0
                LFPinWindow = self.LFP(:,floor(startTime * self.lfpSamplingRate):floor(endTime * self.lfpSamplingRate));
            elseif ~exist('var','startTime') && ~exist('var','endTime')
                LFPinWindow = self.LFP;
            end
            
        end
        
        function extractLFP(self,wb)
            % This extracts LFP signals from N channels of wideband recordings
            % The downsampling rate, low-pass and high-pass edges of the
            % band-pass filter are hardcoded.
            
            scale        = 30;      % this scale changes the sampling rate to 1k
            low          = 0.01;     % Hz
            high         = 200;     % Hz
            
            
            sr = self.SamplingRate / scale;
            Ch = self.channels;
            for chcount = Ch
                fprintf(['channel #',num2str(chcount), ': '])
                data_temp = wb.getDataByTime(chcount);
                neur_ds = Resample(self,data_temp,scale);
                clear data_temp
                self.LFP(chcount,:) = bpfilter(self,neur_ds',low,high,sr);
            end
            self.lfpSamplingRate = sr;
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
            fprintf('Band-pass filtering --> \n')
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
        function spikesCount = getSpikesCount(self,SpikeTimes)
            % This function gets the spike times, and determines the spike
            % counts within prespecified bins. 
            % The structure of SpikeTimes input is assumed to be similar to
            % the output of getDataChunks function.
            % timelength and binwidth in seconds
            
            NumConditions = length(SpikeTimes);
            spikesCount = cell(NumConditions,1);
            for condcount = 1:NumConditions
                thisCondData = SpikeTimes{condcount};
                NumEvents = length(thisCondData);
                for ecount = 1:NumEvents
                    thisEventData = thisCondData{ecount};
                    NumNeurons = length(thisEventData);
                    for neurcount = 1:NumNeurons
                        thisSpikeTrain = thisEventData{neurcount};
                        spikesCount{condcount}(ecount,neurcount) = length(thisSpikeTrain);
                
                    end
                end
            end
            
        end
        function Corr = getCorrelation(self,Data)
            % this function, for any input data (time-series), calculates
            % the correlation between the channels/neurons. 
            % It is assumed that the input Data has the structure of the
            % getDataChunks function output.
            NumConditions = length(Data);
            Corr = cell(NumConditions,1);
            for condcount = 1:NumConditions
                thisCond = Data{condcount};
                numEvents = size(thisCond,1);
                numChannels = size(thisCond,2);
                if numChannels == 1
                    error('More than one channel is needed to calculate the correlation matrix');
                end
                
                for ecount = 1:numEvents
                    thisEvent = squeeze(thisCond(ecount,:,:));
                    Corr{condcount}(ecount,:,:) = corrcoef(thisEvent');
                    
                end
            end
            
            
            
        end
    end
end



