%% Description
%
% The following script extracts the RR intervals prior, during and after
% sound omissions to identify cardiac omission responses from ECG signals
% for the cardio audio experiment as described in 
% https://doi.org/10.1101/2022.03.03.482861
%
% %
% Requires modification of data paths according to user
%
% %
% Inputs
%
% FOR EACH EXPERIMENTAL CONDITION AND PARTICIPANT
% T variable is a table array of dimension Yx3, where Y is the number of 
% trials per experimental block and 3 columns as follows:
% T.trial = trial number per block 
% T.cond = binary variable (1 for sound trial; 2 for omission trial)
% T.lsound = actual sound trigger onset in samples
% loaded in this script as e.g. '../trg_S1.mat'
% example
% trial     cond    lsound
% 1         1       9173
% 2         1       10064
% 3         1       10984
% 4         1       11945
% 5         2       NaN
%
% r variable (Nx1 vector) containing all offline detected time stamps of 
% R peaks (from raw ECG data timeseries)
% loaded in this script as e.g. '../ecg_S1.mat'
%
% %
% Output
%
% awake structure containing RR intervals for omission trials as follows:
% awake.RR_omis_m1 = average RR interval for trial prior to omission
% awake.RR_omis = average RR interval for omission trial
% awake.RR_omis_p1 = average RR interval for trial after omission
% awake.RR_omis_p2 = average RR interval for two trials after omission
% awake.RR_omis_XX variable (3xN matrix) with rows for each condition (synch,
% asynch, isoch) and columns for N subjects

% awake.RRall_omis_m1 = all RR intervals for trial prior to omission
% awake.RRall_omis = all RR intervals for omission trial
% awake.RRall_omis_p1 = all RR intervals for trial after omission
% awake.RRall_omis_p2 = all RR intervals for two trials after omission
% awake.RRall_omis_XX cell structure (3xN) with rows for each condition (synch,
% asynch, isoch) and columns for N sujects
% each cell is a Yx1 vector containing all available RR interval values
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
%set path
if OSflag==1
    mpath=['/mnt/HDD1/CardioAudio_sleep/data/awake/'];
    mpath_analysis=['/mnt/HDD1/CardioAudio_sleep/data/'];
elseif OSflag==2
    mpath=['D:\CardioAudio_sleep\data\awake\'];
    mpath_analysis=['D:\CardioAudio_sleep\'];
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
        Tomis=[]; %R to R interval during omission
        Tomis_m1=[]; %R to R interval before omission
        Tomis_p1=[]; %R to R interval after omission
        Tomis_p2=[]; %2 R to R interval after omission
        
        for block=block_num
            
            if OSflag==1
            	load([mpath,'s',num2str(subjs(j)),'/process/trg_',condition{k},num2str(block),'.mat']) % trigger info before artefact rejection
                load([mpath,'s',num2str(subjs(j)),'/process/ecg_',condition{k},num2str(block),'.mat']) % ecg info before artefact rejection
            elseif OSflag==2
                load([mpath,'s',num2str(subjs(j)),'\process\trg_',condition{k},num2str(block),'.mat']) % trigger info before artefact rejection
                load([mpath,'s',num2str(subjs(j)),'\process\ecg_',condition{k},num2str(block),'.mat']) % ecg info before artefact rejection
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
            
            %remove SOSO trials
            T_temp = T;
            is = [];
            io=find(T.cond(2:end-3) == 2);
            io=io+1;
            for i=1:length(io)
            if T.cond(io(i,1))==2 && T.cond((io(i,1))+2)==2        
                is = [is;io(i,1)-1;io(i,1);io(i,1)+1;io(i,1)+2];
            end
            end
            
            T_temp.lomiss(is,:) = NaN;
            T = T_temp;
            clear T_temp io is;
            
            %start at 3 trial
            io=find(isnan(T.lomiss(3:end-3)) == 0);    
            io=io+2;
            
            %identify omission R peaks
            T_temp = T;
            for i=1:length(io)
                valom = round(T_temp.lomiss(io(i,1)));
                [~,idx]=min(abs(r-valom));
                minval = r(idx);
                if minval > valom
                    Tomis = [Tomis; r(idx) - r(idx-1)];
                    Tomis_m1=[Tomis_m1; r(idx-1) - r(idx-2)]; 
                    Tomis_p1=[Tomis_p1; r(idx+1) - r(idx)]; 
                    Tomis_p2=[Tomis_p2; r(idx+2) - r(idx+1)]; 
                else
                    %elseif minval < valom
                    Tomis = [Tomis; r(idx+1) - r(idx)];
                    Tomis_m1=[Tomis_m1; r(idx) - r(idx-1)]; 
                    Tomis_p1=[Tomis_p1; r(idx+2) - r(idx+1)]; 
                    Tomis_p2=[Tomis_p2; r(idx+3) - r(idx+2)]; 
                end
            end          
            clear io r T_temp           
            
        end
            
            %identify and remove outliers
            if exclude_outliers==1
            Tomisplot=Tomis;
            Tomisplot(find(isnan(Tomis) == 1))=[];
            [TomisClean,II1]=rmoutliers(Tomisplot);
            Tomism1plot=Tomis_m1;
            Tomism1plot(find(isnan(Tomis_m1) == 1))=[];
            [Tomism1Clean,II2]=rmoutliers(Tomism1plot);
            Tomisp1plot=Tomis_p1;
            Tomisp1plot(find(isnan(Tomis_p1) == 1))=[];
            [Tomisp1Clean,II3]=rmoutliers(Tomisp1plot);
            Tomisp2plot=Tomis_p2;
            Tomisp2plot(find(isnan(Tomis_p2) == 1))=[];
            [Tomisp2Clean,II4]=rmoutliers(Tomisp2plot);

            end

        
        if exclude_outliers==1
        %ensure equal numbers of distance to omission trials following outlier removal    
        min_omistrl = min([length(TomisClean); length(Tomism1Clean); length(Tomisp1Clean); length(Tomisp2Clean)]);
        %Store subject-wise outlier free data
        RRall_omis{k,j}=TomisClean(1:min_omistrl,1);
        RRall_omis_m1{k,j}=Tomism1Clean(1:min_omistrl,1);
        RRall_omis_p1{k,j}=Tomisp1Clean(1:min_omistrl,1);
        RRall_omis_p2{k,j}=Tomisp2Clean(1:min_omistrl,1);
        
        %store group outlier free data
        RR_omis(k,j)=mean(RRall_omis{k,j},1);
        RR_omis_m1(k,j)=mean(RRall_omis_m1{k,j},1);
        RR_omis_p1(k,j)=mean(RRall_omis_p1{k,j},1);
        RR_omis_p2(k,j)=mean(RRall_omis_p2{k,j},1);
        
        else %don't exclude outliers
        %ensure equal numbers of distance to omission trials following outlier removal    
        min_omistrl = min([length(Tomis); length(Tomis_m1); length(Tomis_p1); length(Tomis_p2)]);
        %Store subject-wise data
        RRall_omis{k,j}=Tomis(1:min_omistrl,1);
        RRall_omis_m1{k,j}=Tomis_m1(1:min_omistrl,1);
        RRall_omis_p1{k,j}=Tomis_p1(1:min_omistrl,1);
        RRall_omis_p2{k,j}=Tomis_p2(1:min_omistrl,1);
        
        %store group data
        RR_omis(k,j)=mean(RRall_omis{k,j},1);
        RR_omis_m1(k,j)=mean(RRall_omis_m1{k,j},1);
        RR_omis_p1(k,j)=mean(RRall_omis_p1{k,j},1);
        RR_omis_p2(k,j)=mean(RRall_omis_p2{k,j},1);
       
        end
        
        clear Ts Tr dTs soundR Rsound RRint
    
    end
end

%store relevant data for statistical analysis
awake = [];
awake.RR_omis_m1 = RR_omis_m1;
awake.RR_omis = RR_omis;
awake.RR_omis_p1 = RR_omis_p1;
awake.RR_omis_p2 = RR_omis_p2;
awake.RRall_omis_m1 = RRall_omis_m1;
awake.RRall_omis = RRall_omis;
awake.RRall_omis_p1 = RRall_omis_p1;
awake.RRall_omis_p2 = RRall_omis_p2;
if exclude_outliers==1
    save awake_CR_all awake
else
    save awake_CR_withoutliers_all awake
end
end