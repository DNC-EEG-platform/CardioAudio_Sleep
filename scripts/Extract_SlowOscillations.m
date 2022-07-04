%% Description
% The following script extracts Slow Oscillations for the cardio audio 
% experiment described in https://doi.org/10.1101/2022.03.03.482861 based 
% on the method described by Ngo et al., 2013 (DOI: 10.1016/j.neuron.2013.03.006) 
% & Besedovsky et al., 2017 (DOI: 10.1038/s41467-017-02170-3 ) adapted from:
% https://data.mendeley.com/datasets/c872f83pdz/4
%
% Slow oscillations detected by identified positive-to-negative zero
% crossing in timeseries if they are within a certain latency to each other 
% (e.g. here 0.833 to 2s here) of each other and if the negative peak within a zero 
% crossing satisfies specific amplitude thresholds set by the user (e.g. here
% 1.25 lower than mean negative amplitude and 1.25 higher than mean
% positive to negative amplitude difference)
%
% %
% Requires modification of data paths according to user
%
% %
% Dependencies
% Fieldtrip toolbox (https://www.fieldtriptoolbox.org/download)
%
% %
% Input
%
% data = fieldtrip data structure containing segmented condition-specific
% subject-wise EEG data.
% data is the output of fieldtrip's ft_preprocessing that has been:
% 1. high pass filtered at 3.5Hz using FIR fieldtrip
% 2. downsampled to 100Hz
% 3. cleared from artefacted trials (labeled by the preprocessing of the
% entire dataset for each subject)
% data structure with fields:
%         label: {65×1 cell}         --- electrode labels (i.e. 62 EEG, 2 EOG, 1 ECG)
%    sampleinfo: [Yx3 double]        --- matrix containing timestamps for
%                                        beginning timeframe (column 1), end timeframe (column 2) 
%                                        and block number (column 3) for Y number of trials (rows)                                          
%         trial: {1xY cell}          --- Y cells containing trial-wise actual timeseries data 
%                                        (channel x time)
%          time: {1xY cell}          --- Y cells containing trial corresponding time info in seconds
%                                        (channel x time)
%           cfg: [1×1 struct]        --- cfg from fieldTrip (log-keeping), can be empty at beginning
%       fsample: 100                 --- sampling rate of downsampled recording
%
% loaded in this script as e.g. '../S_ft_filt_slowosc_N2.mat'
%
% Output
%
% evt_info = trial x channel structure containing the condition-specific subject-wise 
%            detected SO information
% each structure entry (for each trial and channel) with fields:
%        numEvt: [1x1 vector]        --- contains the number of events (N) for each trial
%     startTime: [1XN double]        --- matrix containing start timestamps of N detected SOs for each trial     
%       endTime: [1XN double]        --- matrix containing end timestamps of N detected SOs for each trial 
%      duration: [1XN double]        --- matrix containing duration in seconds of N detected SOs for each trial 
%       minTime: [1XN double]        --- matrix containing negative peak timestamps of N detected SOs for each trial
%       maxTime: [1XN double]        --- matrix containing positive peak timestamps of N detected SOs for each trial
%   minroigTime: [1XN double]        --- matrix containing negative peak timestamps of N detected SOs for each trial
%                                        for the entire length of each block for later sleep stage labelling                                         
%   maxroigTime: [1XN double]        --- matrix containing positive peak timestamps of N detected SOs for each trial
%                                        for the entire length of each block for later sleep stage labelling                                         
%        minAmp: [1XN double]        --- matrix containing negative peak amplitude of N detected SOs for each trial
%        maxAmp: [1XN double]        --- matrix containing positive peak amplitude of N detected SOs for each trial
%
% %
%     Copyright (C) 2022, Andria Pelentritou
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% %
%
%% Code
clear all;
OSflag=1; %1 = Linux; 2 = Windows

%clear out fieldtrip dependencies
restoredefaultpath
if OSflag==1
    addpath /mnt/HDD1/CardioAudio_sleep/scripts/fieldtrip-20201205/; 
    addpath /mnt/HDD1/CardioAudio_sleep/scripts/ScriptsLinux/;
    addpath /mnt/HDD1/CardioAudio_sleep/data/
elseif OSflag==2
	addpath D:\CardioAudio_sleep\scripts\fieldtrip-20201205\; 
    addpath D:\CardioAudio_sleep\scripts\Scripts\;
    addpath D:\CardioAudio_sleep\data\;
end

%restore fieldtrip defaults
ft_defaults

%valid subjects
subjs=[3:6 8:26]; %included subjects
%valid conditions
condition{1}='S'; %synch
condition{2}='A'; %asynch
condition{3}='O'; %isoch
condition{4}='B'; %baseline
%required sampling rate
fsample=100;
%length of positive-to-negative zero crossing (one full SO)
thresh_lngth=[0.833 2];
%negative peak amplitude threshold for accepting an SO
slct_thresh=1.25;
%channel list for SO identification
slct_chans={'Cz'};

%set path
for sub=1:length(subjs)
    for cond=1:length(condition)
        if OS_flag==1
            sub_dir = ['/mnt/HDD1/CardioAudio_sleep/data/sleep/s',num2str(subjs(sub)),'/process/']; 
        elseif OS_flag==2
            sub_dir = ['D:\CardioAudio_sleep\data\sleep\s',num2str(subjs(sub)),'\process\']; 
        end
cd(sub_dir)    

filename = sprintf('%s_ft_filt_slowosc_N2.mat',condition{cond}); % 3.5Hz low pass FIR filtered data
load(filename,'data')

%only include selected channel data
cfg=[];
cfg.channel=slct_chans;
data=ft_selectdata(cfg,data);

evt_info=[];
evt_summary=[];

%identify positive-to-negative zero crossings
pos2negZC = cell(size(data.trial,2), size(slct_chans,1));
neg2posZC = cell(size(data.trial,2), size(slct_chans,1));
for i_trl=1:size(data.trial,2)
    for i_ch=1:size(slct_chans,1)
    signdata=sign(data.trial{1,i_trl}(i_ch,:));
    signdiff=diff(signdata);
    
    if size(signdiff,2)<2
    else
    pos2negZC{i_trl,i_ch} = uint32(find(signdiff < 0));
    neg2posZC{i_trl,i_ch} = uint32(find(signdiff > 0));
    if or(isempty(pos2negZC{i_trl,i_ch})==1,isempty(neg2posZC{i_trl,i_ch})==1)
    else
    if pos2negZC{i_trl,i_ch}(1) > neg2posZC{i_trl,i_ch}(1) 
       neg2posZC{i_trl,i_ch}(1) = [];
    end
    end
    end
    end
end

%set threshold for valid SOs
SOnegVal = cell(size(data.trial,2), size(slct_chans,1));
SOnegPos = cell(size(data.trial,2), size(slct_chans,1));
SOposVal = cell(size(data.trial,2), size(slct_chans,1));
SOposPos = cell(size(data.trial,2), size(slct_chans,1));
SOnegorigPos= cell(size(data.trial,2), size(slct_chans,1));
SOposorigPos= cell(size(data.trial,2), size(slct_chans,1));
SOpk2pkVal = cell(size(data.trial,2), size(slct_chans,1));
SOlngth = cell(size(data.trial,2), size(slct_chans,1));
for i_trl=1:size(data.trial,2)
    for i_ch=1:size(slct_chans,1)
    for i_ZC=1:numel(pos2negZC{i_trl,i_ch})-1
        [SOnegVal{i_trl,i_ch}(i_ZC), SOnegPos{i_trl,i_ch}(i_ZC)]=min(data.trial{1,i_trl}(i_ch,pos2negZC{i_trl,i_ch}(i_ZC):neg2posZC{i_trl,i_ch}(i_ZC)));
        [SOposVal{i_trl,i_ch}(i_ZC), SOposPos{i_trl,i_ch}(i_ZC)]=max(data.trial{1,i_trl}(i_ch,neg2posZC{i_trl,i_ch}(i_ZC):pos2negZC{i_trl,i_ch}(i_ZC+1)));
        SOpk2pkVal{i_trl,i_ch}(i_ZC)=abs(SOnegVal{i_trl,i_ch}(i_ZC))+SOposVal{i_trl,i_ch}(i_ZC);
        SOnegPos{i_trl,i_ch}(i_ZC)=SOnegPos{i_trl,i_ch}(i_ZC)+pos2negZC{i_trl,i_ch}(i_ZC)-1;
        SOposPos{i_trl,i_ch}(i_ZC)=SOposPos{i_trl,i_ch}(i_ZC)+neg2posZC{i_trl,i_ch}(i_ZC)-1;
        SOnegorigPos{i_trl,i_ch}(i_ZC)=SOnegPos{i_trl,i_ch}(i_ZC)+data.sampleinfo(i_trl,1);
        SOposorigPos{i_trl,i_ch}(i_ZC)=SOposPos{i_trl,i_ch}(i_ZC)+data.sampleinfo(i_trl,1);
        SOlngth{i_trl,i_ch}(i_ZC)=(double(pos2negZC{i_trl,i_ch}(i_ZC+1)-pos2negZC{i_trl,i_ch}(i_ZC)))/fsample;
    end
    end
end

SOcritVal   = SOnegVal;
SOcritVal2  = SOpk2pkVal;

vardata_min	= cell(size(slct_chans,1),1);
vardata_max	= cell(size(slct_chans,1),1);

for i_trl=1:size(data.trial,2)
    for i_ch=1:size(slct_chans,1)
        counter=size(vardata_min{i_ch,1},2);
    for i_evt=1:numel(SOcritVal{i_trl,i_ch})
        if SOlngth{i_trl,i_ch}(i_evt)>=thresh_lngth(1,1) && SOlngth{i_trl,i_ch}(i_evt)<=thresh_lngth(1,2)
            counter = counter+1;
            vardata_min{i_ch,1}(counter)=abs(SOcritVal{i_trl,i_ch}(i_evt));
            vardata_max{i_ch,1}(counter)=abs(SOcritVal2{i_trl,i_ch}(i_evt));
        end
    end
    end
end

%select valid SOs based on threshold
for i_ch=1:size(slct_chans,1)
    thresh_min(i_ch)=mean(vardata_min{i_ch,1},2);
    thresh_max(i_ch)=mean(vardata_max{i_ch,1},2);

    SOthresh_min(i_ch)=thresh_min(i_ch).*slct_thresh;
    SOthresh_max(i_ch)=thresh_max(i_ch).*slct_thresh;
end

for i_trl=1:size(data.trial,2)
    for i_ch=1:size(slct_chans,1)
    counter_evt=0;
for f_evt = 1:numel(SOcritVal{i_trl,i_ch})                                                   
	f_evt_idx = pos2negZC{i_trl,i_ch}(f_evt) - round(1*fsample) : pos2negZC{i_trl,i_ch}(f_evt+1) + round(1*fsample);
	f_evt_idx = f_evt_idx(f_evt_idx > 0 & f_evt_idx <= data.sampleinfo(i_trl,2));

    if abs(SOcritVal{i_trl,i_ch}(f_evt)) >= SOthresh_min(i_ch)
        test_crit=abs(SOcritVal2{i_trl,i_ch}(f_evt)) >= SOthresh_max(i_ch);
        if test_crit
            if SOlngth{i_trl,i_ch}(f_evt)>=thresh_lngth(1,1) && SOlngth{i_trl,i_ch}(f_evt) <= thresh_lngth(1,2)
                counter_evt=counter_evt+1;
                evt_info(i_trl,i_ch).startTime(counter_evt)=double(pos2negZC{i_trl,i_ch}(f_evt));
                evt_info(i_trl,i_ch).endTime(counter_evt)=double(pos2negZC{i_trl,i_ch}(f_evt+1));
                
                evt_info(i_trl,i_ch).duration(counter_evt)=SOlngth{i_trl,i_ch}(f_evt);
                
                evt_info(i_trl,i_ch).minTime(counter_evt)=SOnegPos{i_trl,i_ch}(f_evt);
                evt_info(i_trl,i_ch).maxTime(counter_evt)=SOposPos{i_trl,i_ch}(f_evt);
                
                evt_info(i_trl,i_ch).minorigTime(counter_evt)=SOnegorigPos{i_trl,i_ch}(f_evt);
                evt_info(i_trl,i_ch).maxorigTime(counter_evt)=SOposorigPos{i_trl,i_ch}(f_evt);
                
                evt_info(i_trl,i_ch).minAmp(counter_evt)=SOnegVal{i_trl,i_ch}(f_evt);
                evt_info(i_trl,i_ch).maxAmp(counter_evt)=SOposVal{i_trl,i_ch}(f_evt);
            end
        end
    end
end
evt_info(i_trl,i_ch).numEvt = counter_evt;   % Save number of detected SO
    end
end
        filename = sprintf('%s_SO_N2.mat',condition{cond}); %save data
        save(filename,'evt_info','data','-v7.3')
    end
end