%% Description
%
% The following script performs the control analyses for the cardio audio 
% experiment described in https://doi.org/10.1101/2022.03.03.482861
% SR, RS, SS, RR intervals and variabilities are calculated to assess 
% whether the experimental design was successfuly administered
%
% %
% Requires modification of data paths according to user
%
% %
% Inputs
%
% FOR EACH AUDITORY CONDITION AND PARTICIPANT
% T variable is a table array of dimension Yx4, where Y is the number of 
% trials per experimental block and 4 columns as follows:
% T.trial = trial number per block 
% T.cond = binary variable (1 for sound trial; 2 for omission trial)
% T.lsound = actual sound trigger onset in samples
% T.lr = closest r peak onset in samples (to stimulus onset)
% loaded in this script as e.g. '../trg_S1.mat'
% example
% trial     cond    lsound      lr
% 1         1       9173        9105
% 2         1       10064       9999
% 3         1       10984       10919
% 4         1       11945       11875
% 5         2       NaN         13895
%
% FOR BASELINE CONDITION AND PARTICIPANT
% rr variable (Nx1 vector) containing all RR intervals
% calculated as the difference between all offline detected time stamps of 
% R peaks (from raw ECG data timeseries)
% loaded in this script as e.g. ..'/ecg_B1.mat'
%
% %
% Output
% 
% awake.SR = mean SR interval
% awake.RS = mean RS interval
% awake.SS = mean SS interval
% awake.RR = mean RR interval
% awake.SR, awake.RS, awake.SS variable (3xN matrix) with rows for each 
% condition (synch, asynch, isoch) and columns for N sujects
% awake.RR variable (4xN matrix) with rows for each 
% condition (synch, asynch, isoch, baseline) and columns for N sujects
%
% awake.varSR = mean SR variability
% awake.varRS = mean RS variability
% awake.varSS = mean SS variability
% awake.varRR = mean RR variability
% variability refers to SEM of corresponding intervals 
% awake.varSR, awake.varRS, awake.varSS variable (3xN matrix) with rows for each 
% condition (synch, asynch, isoch) and columns for N sujects
% awake.varRR variable (4xN matrix) with rows for each 
% condition (synch, asynch, isoch, baseline) and columns for N sujects
%
% awake.SRall = SRall; %all SR interval values 
% awake.RSall = RSall; %all RS interval values
% awake.SSall = SSall; %all SS interval values
% awake.RRall = final_RRall; %all RR interval values
% awake.SRall, awake.RSall, awake.SSall cell structure (3xN) with rows for each 
% condition (synch, asynch, isoch) and columns for N sujects
% awake.RRall cell structure (4xN) with rows for each 
% condition (synch, asynch, isoch, baseline) and columns for N sujects
% each cell is a Yx1 vector containing all available interval values
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
%% Code
clear all;
rng(1,'twister'); %fixed randomization

OSflag=1; %1 = Linux; 2 = Windows
% set path
if OSflag==1
    mpath=['/mnt/HDD1/CardioAudio_sleep/data/awake/']; %
    mpath_analysis=['/mnt/HDD1/CardioAudio_sleep/data/'];
elseif OSflag==2
    mpath=['D:\sleep\data\awake\'];
    mpath_analysis=['D:\sleep\'];
end

%valid subjects
subjs=[1:26]; %included subjects
%valid conditions
condition{1}='S'; %synch
condition{2}='A'; %asynch
condition{3}='O'; %isoch
%sampling rate
Fs=1200;
%block number
block_num=[1:6]; %included blocks

compute_means=1; %1 for yes; 0 for no
exclude_outliers=1; %1 for yes; 0 for no
 
            
%% Check for changes across subjects
if compute_means==1
for k=1:length(condition)
    for j=1:length(subjs)
        
        Ts=[]; %sound onset
        SoundR=[]; %sound to R interval
        dTs=[]; % sound to sound interval
        RSound=[]; %R to sound interval
        RRint=[]; %%R to R interval
        condstage=[]; %block number

        for block=block_num
            
            if OSflag==1
            	load([mpath,'s',num2str(subjs(j)),'/process/trg_',condition{k},num2str(block),'.mat']) % trigger info before artefact rejection
            elseif OSflag==2
                load([mpath,'s',num2str(subjs(j)),'\process\trg_',condition{k},num2str(block),'.mat']) % trigger info before artefact rejection
            end
            
            %prepare T variable
            s = T.lsound(T.cond==1);
            ss = diff(T.lsound);
            ss_avg = nanmean(ss);

            T.lomiss = nan(size(T.lsound));
            ind = find(T.cond==2);
            ind=ind(ind>1);
            T.lomiss(ind) = T.lsound(ind-1)+ss_avg;
            clear ind
            
            %compute different control analyses variables
            Ts=[Ts; T.lsound]; %S latencies
            
            io=find(T.cond == 2);
            io(1)=[];
 
            SoundR=[SoundR; T.lr(io)-T.lsound(io-1)]; %SR

            ir=find(T.cond == 1);
            ir(1)=[];
            
            RSound=[RSound; T.lsound(ir)-T.lr(ir)]; %RS
            
            RRint=[RRint; diff(T.lr)]; %RR
            
            dTs=[dTs; diff(T.lsound)]; %SS     
            
            %keep SS block information for blockwise SS variability calculation
            keep_length=length(diff(T.lsound));
            cond_blocks = ones(keep_length,1);
            cond_blocks(:,1)=block;
            condstage= [condstage; cond_blocks];

            clear io ir cond_blocks keep_length
            
        end
        
            if exclude_outliers==1
            %identify and remove outliers
            dTsplot=dTs; %SS outliers - currently not used further
            dTsplot(find(isnan(dTs) == 1))=[]; 
            [dTsClean,II1]=rmoutliers(dTsplot);
            SoundRplot=SoundR; %SR outliers
            SoundRplot(find(isnan(SoundR) == 1))=[];
            [SoundRClean,II2]=rmoutliers(SoundRplot);
            RSoundplot=RSound; %RS outliers
            RSoundplot(find(isnan(RSound) == 1))=[];
            [RSoundClean,II3]=rmoutliers(RSoundplot);
            RRintplot=RRint; %RR outliers
            RRintplot(find(isnan(RRint) == 1))=[];
            [RRintClean,II4]=rmoutliers(RRintplot);

            end
        

        if exclude_outliers==1
        %keep all values
        SRall{k,j}=SoundRClean;
        RSall{k,j}=RSoundClean;
        %SSall{k,j}=dTsClean;
        SSall{k,j}=dTs; %using the non outlier corrected data
        RRall{k,j}=RRintClean;
        
        SR(k,j)=nanmean(SoundRClean); %SR interval
        varSR(k,j)=(nanstd(SoundRClean))/sqrt(length(SoundRClean)); %SR variability

            %SS variability calculated for each block
            unique_blocks=unique(condstage);
            for i_un=1:length(unique_blocks)
            iunique=find(condstage == unique_blocks(i_un));
            sl_var=dTs(iunique);
            var_temp=nanstd(sl_var)/sqrt(sum(~isnan(sl_var)));
            end
            varSS(k,j)=nanmean(var_temp);
            var_temp=[];
        
        SS(k,j)=nanmean(dTs); %SS interval
        varSS(k,j)=nanmean(dTs); 
        varSS(k,j)=sqrt(varSS(k,j)); %SS variability using the non outlier corrected data 
        
        RS(k,j)=nanmean(RSoundClean); %RS interval
        varRS(k,j)=(nanstd(RSoundClean))/sqrt(length(RSoundClean)); %RS variability

        RR(k,j)=nanmean(RRintClean); %RR interval
        varRR(k,j)=(nanstd(RRintClean))/sqrt(length(RRintClean)); %RR variability

        
        else %if exclude_outliers==0
        SRall{k,j}=SoundR;
        SSall{k,j}=dTs;
        RSall{k,j}=RSound;
        RRall{k,j}=RRint;
        
        SR(k,j)=nanmean(SoundR); %SR interval
        varSR(k,j)=(nanstd(SoundR))/sqrt(length(SoundR)); %SR variability

            %SS variability calculated for each block
            unique_blocks=unique(condstage);
            for i_un=1:length(unique_blocks)
            iunique=find(condstage == unique_blocks(i_un));
            sl_var=dTs(iunique);
            var_temp=nanstd(sl_var)/sqrt(sum(~isnan(sl_var)));
            end
            varSS(k,j)=nanmean(var_temp);
            var_temp=[];
        
        SS(k,j)=nanmean(dTs); %SS interval
        varSS(k,j)=nanmean(dTs); 
        varSS(k,j)=sqrt(varSS(k,j)); %SS variability using the non outlier corrected data 
        
        RS(k,j)=nanmean(RSound); %RS interval
        varRS(k,j)=(nanstd(RSound))/sqrt(length(RSound)); %RS variability
        
        RR(k,j)=nanmean(RRint); %RR interval
        varRR(k,j)=(nanstd(RRint))/sqrt(length(RRint)); %RR variability
        
        end

        clear Ts Tr dTs SoundR RSound RRint dTsClean SoundRClean RSoundClean RRintClean
    
    end
end

%include baseline data for RR info
for j=1:length(subjs)
    
    clear irr
    if OSflag==1
        load([mpath,'s',num2str(subjs(j)),'/process/ecg_B1.mat']) % ECG info before artefact rejection
    elseif OSflag==2
        load([mpath,'s',num2str(subjs(j)),'\process\ecg_B1.mat']) % ECG info before artefact rejection
    end
    
     if exclude_outliers==1
     %identify and remove outliers
     	rrbplot=rr;            
     	rrbplot(find(isnan(rr) == 1))=[];
     	[rrbClean,II1]=rmoutliers(rrbplot);
     
        %ensure baseline trials match condition trials
        max_btrl=max([length(RRall{1,j}); length(RRall{2,j}); length(RRall{3,j})]); %synch asynch isoch trial numbers
        if length(rrbClean)>max_btrl
            final_btrl=(randperm(length(rrbClean),max_btrl))';
            final_trlb=rrbClean(final_btrl,1);
        else
            final_trlb=rrbClean;
        end
        
     else %if exclude outliers==0
        %ensure baseline trials match condition trials
        max_btrl=max([length(RRall{1,j}); length(RRall{2,j}); length(RRall{3,j})]); %synch asynch isoch trial numbers
        if length(rr)>max_btrl
            final_btrl=(randperm(length(rr),max_btrl))';
            final_trlb=rr(final_btrl,1);
        else
            final_trlb = rr;
        end
     end  
     
        %store baseline info
        RRall{4,j}=final_trlb; 
        
        clear rr
        
        %use equal trials for RR where trial length is of relevance
        min_trl=min([length(RRall{1,j}); length(RRall{2,j}); length(RRall{3,j}); length(RRall{4,j})]); %synch asynch isoch baseline trial numbers
        for k=1:4 %number of conditions
        final_RRtrl=(randperm(length(RRall{k,j}),min_trl))';
        final_RRall{k,j}=RRall{k,j}(final_RRtrl,1);
        final_RR(k,j)=nanmean(final_RRall{k,j});
        final_varRR(k,j)=(nanstd(final_RRall{k,j}))/sqrt(length(final_RRall{k,j}));

        end
end

%store relevant data for statistical analysis
awake.SR = SR; %mean SR interval
awake.RS = RS; %mean RS interval
awake.SS = SS; %mean SS interval
awake.RR = final_RR; %mean RR interval
awake.SRall = SRall; %all SR interval values 
awake.RSall = RSall; %all RS interval values
awake.SSall = SSall; %all SS interval values
awake.RRall = final_RRall; %all RR interval values
awake.varSR = varSR; %mean SR variability
awake.varRS = varRS; %mean RS variability
awake.varSS = varSS; %mean SS variability
awake.varRR = final_varRR; %mean RR variability
if exclude_outliers==1
    save awake_CA_all awake
else
    save awake_CA_withoutliers_all awake
end
end