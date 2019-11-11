
%% This script performs GED AND CFC analyses
% Finds theta phase of maximum power RT correlation for each subject
% Uses data around the phase of maximum correlation for the S matrix,
% broadband data from thetamaxfrequ+ 2 to 80Hz as the R matrix.
% Theta phase-frequency power is then computed (phases are shifted so that
% phase = 0 is actually the phase of maximum correlation for a given
% subject.

% This version of the script applies a Cross Validation approach that uses
% 90 % of the trials to define the spatial filter that is only applied to
% the 10% remaining data. The scripts loops so that in the end all of the
% data is used as a the test sample.


%% Subjects
subjects = {
    'pp01'      [ ];
    'pp02'      [ ];
    'pp03'      [ ];
    'pp04'      [ ];
    'pp05'      [ ];
    'pp06'      [ ];
    'pp07'      [ ];
    'pp08'      [ ];
    'pp09'      [ ];
    'pp11'      [ ];
    'pp12'      [ ];
    'pp13'      [ ];
    'pp14'      [ ];
    'pp15'      [ ];
    'pp16'      [ ];
    'pp17'      [ ];
    'pp18'      [ ];
    'pp19'      [ ];
    'pp20'      [ ];
    'pp21'      [ ];
    'pp23'      [ ];
    'pp24'      [ ];
    'pp25'      [ ];
    'pp26'      [ ];
    'pp27'      [ ];
    'pp28'      [ ];
    'pp29'      [ ];
    'pp30'      [ ];
    }; % pp 10 and 22 excluded

%% Import fitted max and min correlation phase values

% load phases of min and max correlations

maxminph = zeros(size(subjects, 1),2,5); % Initialize matrix that will contain phase bins

for condi = 1:5
    maxph = load(['emp_maxminph13_300-1200_' num2str(condi) '.txt']); % first column is phase with min correlation, second is max
    maxminph(:,:,condi) = maxph;
end

phases = linspace(-pi,pi,20); % vector containing the different phase values

maxminphases = zeros(size(subjects, 1),2,5); % Initialize matrix that will contain real phase values

for celi = 1:(28*2*5)
    maxminphases(celi) = phases(maxminph(celi));
end
maxminphases = double(maxminphases);

%% Convolution parameters - used for TF decomposition after GED

% Soft coded convolution parameters for phase-frequency power --> to see
% high frequency power (>10Hz) according to theta component phase.
min_frex  =  10;
max_frex  = 50;
num_frex  = 40;
nCycRange = [6 12];

frex = logspace(log10(min_frex),log10(max_frex),num_frex);

% Soft coded convolution parameters for time-frequency power
min_frex_tf  =  1;
max_frex_tf  = 50;
num_frex_tf  = 60;
nCycRange_tf = [4 12];

frex_tf = logspace(log10(min_frex_tf),log10(max_frex_tf),num_frex_tf);

times2save = -300:15:1200; % for temporal downsampling at 67 Hz
basetime   = [ -500 -200 ]; % baseline period


% Initialize TF matrix of the maxph component for all subjects (to plot average TF maps)

tfmaxph_avg = zeros(10, size(subjects, 1), 5, num_frex,length(times2save)); % Added a first dimension of 10 for the 10 iterations


% Initialize the maxph component matrix for all subjects (to plot average topomap)

maps_maxminph_avg = zeros(10, size(subjects, 1), 64); % Added a first dimension of 10 for the 10 iterations

% Initialize output matrix for phase-frequency power map (sub, nfreq, nbin)

nbin = 51; % instead of 50, needs to be odd to have 0 at the center
tf_phaseresolved = zeros (10, size(subjects, 1), 5, size(frex,2), nbin); % Added a first dimension of 10 for the 10 iterations
tf_phaseresolved_shifted = zeros(10, size(subjects, 1), 5, size(frex,2), nbin); % for when the phases are shifted so that 0 is maxcorr phase
% Added a first dimension of 10 for the 10 iterations

%%
for subno=1:size(subjects,1)
    
    %% Load and initialize
    
    
    %% Load data
    
    load([  subjects{subno,1} '_GED.mat' ])
    
    EEG.data = double(EEG.data); % double converts to double precision
    
    for condi = 1:5
        EEG.data(:,:,confconds==condi) = bsxfun(@minus,EEG.data(:,:,confconds==condi),squeeze(mean(EEG.data(:,:,confconds==condi),3)));
    end % subtracts the phase-locked signal
    
    
    % reject a really weird trial
    if subno == 5
        EEG.thetadata(:,:,133) = []; % suppress trial 133 for pp5
        EEG.data(:,:,133) = []; % suppress trial 133 for pp5
        confconds(133) = [];
    else
    end
    
    thetadata = double(EEG.thetadata);
    
    
    % Time to save and baseline indices
    times2saveidx = dsearchn(EEG.times',times2save');
    baseidx = dsearchn(EEG.times',basetime');
    
    
    % number of trials to test
    ntest = round((10*size(EEG.data, 3)/100)); % 10 to select 10% of the data for testing
    
    % Create a variable that increases of ntest at each iteration
    n2test = 0;
    
    for iteri = 1:10
        
        if iteri <= 9
            trial4test = n2test+1:n2test+ntest;
        else
            trial4test = n2test+1:size(EEG.data, 3);
        end
        
        ntrial = 1:size(EEG.data, 3); % create a trial variable that will be used to isolate the trials to test
        ntrial(trial4test) = 0;
        trial4filt = logical(ntrial);
        
       
        %% Convolution parameters
        
        % wavelet parameters for the theta phase-frequency power analysis
        s = logspace(log10(nCycRange(1)),log10(nCycRange(2)),num_frex)./(2*pi.*frex);
        t = -2:1/EEG.srate:2;
        halfwave = floor((length(t)-1)/2);
        
        % convolution params
        nData = size(trial4test,2)*EEG.pnts;
        nWave = length(t);
        nConv = nData+nWave-1;
        
        % wavelets 
        cmwX = zeros(num_frex,nConv);
        for fi=1:num_frex
            cmw = fft(  exp(1i*2*pi*frex(fi).*t) .* exp( (-t.^2)/(2*s(fi)^2) )  ,nConv);
            cmwX(fi,:) = cmw ./ max(cmw); % amplitude-normalize in the frequency domain
        end
        
        
        % wavelet parameters for the time-frequency power analysis
        s = logspace(log10(nCycRange_tf(1)),log10(nCycRange_tf(2)),num_frex_tf)./(2*pi.*frex_tf);
        
        % wavelets
        cmwX_tf = zeros(num_frex_tf,nConv);
        for fi=1:num_frex_tf
            cmw_tf = fft(  exp(1i*2*pi*frex_tf(fi).*t) .* exp( (-t.^2)/(2*s(fi)^2) )  ,nConv);
            cmwX_tf(fi,:) = cmw_tf ./ max(cmw_tf); % amplitude-normalize in the frequency domain
        end
        
        %% GED
        
        % Redefine confcond so it corresponds to the trials used for
        % filt/test
        confconds2 = confconds(trial4filt);
        
        % Filter-Hilbert to get the angle time series according to theta
        % max frequency --> Needed to define the time windows around max
        % phase that will be used for GED
        
        if iteri <= 9 % different according to iterations, because ndata is not the same size for each iteration
            thfiltdat = reshape( hilbert( filterFGx2(reshape(thetadata(:,:,trial4filt),1,[]),EEG.srate,EEG.th_maxfreq,5)' )',   EEG.pnts, size(thetadata, 3)-ntest); % Hilbert decomposition of the data according to max theta freq
        else
            thfiltdat = reshape( hilbert( filterFGx2(reshape(thetadata(:,:,trial4filt),1,[]),EEG.srate,EEG.th_maxfreq,5)' )',   EEG.pnts, size(thetadata(:,:,trial4filt), 3)); % Hilbert decomposition of the data according to max theta freq
        end
        
        tidx = dsearchn(EEG.times',[300 1200]');% same time window as the one used when investigating power-rt corr in phase bins;
        thfiltdat = thfiltdat(tidx(1):tidx(2),:); % select theta data from that time window
        
        % get the phase angles time series
        thfiltdat_ang = angle(thfiltdat);
        
        % defines the size of the time window selection around the target phase
        % depends on theta max freq and corresponds to 1/4 of a cycle
        timewin = round(((1000/EEG.th_maxfreq)/4)*EEG.srate/1000);
        
        [Smatrix] = deal(zeros(64)); % Initialize one empty 64*64 matrix, do the same for the Rmatrix if based on the phase of min correlation
        
        tempdat0 = EEG.data(1:64,tidx(1):tidx(2),trial4filt); % filter from theta max freq+2 to 80Hz CAREFUL TAKE EEG.data and NOT thetadata !! Otherwise you'll only have the theta component
        tempdat = eegfilt(tempdat0, EEG.srate, EEG.th_maxfreq+2, 80); % eegfilt returns a nchannel*(timepoints*trials) matrix --> needs reshaping
        tempdat = reshape(tempdat, 64, size(tempdat0, 2), size(thfiltdat, 2));
        
        for triali = 1:size(thfiltdat, 2)
            
            condi = confconds2(triali);
                       
            maxph = maxminphases(subno,2,condi); % gets the phase of max correlation
            
            % find target phases
            
            maxphidx = find(diff(sign(diff(abs(thfiltdat_ang(:,triali)-maxminphases(subno,2,condi)))))>0)+1; % finds the indices of the local minima
            % (subtract the target phase from the angle time series, take the abs, values closest to 0 are the targets)
            maxphidx(maxphidx<timewin+1) = [];
            maxphidx(maxphidx>size(thfiltdat,1)-timewin) = [];
            
            
            [Sdata] = deal([]); % creates one 0*0 matrix

            for ti=1:length(maxphidx)
                Sdata = [Sdata zscore(tempdat(:,maxphidx(ti)-timewin:maxphidx(ti)+timewin,triali),[],2) ]; % zscore so we don't have to subtract the mean for GED
            end
            
            if ~isempty(Sdata) % ~ indicates 'NOT'; 1 if Sdata is full, 0 if empty
                Smatrix = Smatrix + (Sdata*Sdata')/size(Sdata,2); % Sum Smatrix to itself at each iteration. It will be divided by the number of trials so we get the mean
            end
                      
        end
        
        Smatrix = Smatrix./triali;
        
        % Define the R matrix as the broadband data between thetamaxfreq+2 and
        % 80
        td    = reshape(tempdat(:,:,:),64,[] );
        td    = bsxfun(@minus,td,mean(td,2));
        Rmatrix = (td*td')/size(td,2);  % reference covariance matrix
        
        shr = .001; % regularization parameter
        Rmatrix = (1-shr)*Rmatrix + shr*mean(eig(Rmatrix))*eye(size(Rmatrix));
        
        [evecs_maxminph,evals] = eig( Smatrix, Rmatrix );
        
        [~,sidx] = sort(diag(evals));
        
        maps_maxminph = inv(evecs_maxminph');
        
        % force sign of components --> so that topographical maps always show
        % positive values
        for ci=1:64
            [~,idx] = max(abs(maps_maxminph([9 10 11 47 46 45 44 14 13 12 48 49 50 51], ci))); % find strongest weight by finding the max in all the FCs and Cs electrodes
            maps_maxminph(:,ci) = maps_maxminph(:,ci) * sign(maps_maxminph(idx,ci)); % force to positive sign
        end
        
        % create data and remove components with eyeblinks
        
        comps2keep = 50:64;
        
        maxphdata  = reshape( (reshape(EEG.data(1:64,:,trial4test),64,[])'*evecs_maxminph(:,comps2keep))' ,[length(comps2keep) EEG.pnts size(trial4test, 2)]); % CRUCIAL TO USE trial4test HERE !
        
        maps_maxminph = maps_maxminph(:, comps2keep);
        
        % select components without eyeblinks
        comps2save = (mean(maps_maxminph([11 12 19 47 48 32 46 49 56],:)) - mean(maps_maxminph([1 33 34],:))) > 0; % creates a vector of logical values
        % that compares if activity in central electrodes is higher
        maps_maxminph = maps_maxminph(:, comps2save);                                                                                                                                           % than activity around the eyes (1 if true)
        
        maxphdata = maxphdata(comps2save,:,:);
        
        maxphcompidx = size(maxphdata,1); % Component to keep is the one with highest eigenvalue
        
        %     % Plot the topographical maps of the forward problem
        %
        %     figure(1), clf
        %
        %     for topi = 1:size(maps_maxminph, 2)
        %
        %         subplot(5,3,topi)
        %         topoplot(maps_maxminph(:,topi),EEG.chanlocs,'numcontour', 6, 'gridscale', 200, 'maplimits', 'maxmin')
        %         title([  'Comp ' num2str(topi)  ])
        %
        %     end
        
        maps_maxminph_avg(iteri, subno, :) = maps_maxminph(:, maxphcompidx);
        
                   
        % maxphase component TF power according to theta phase
        
        
        %% Convolution to get theta phase-frequency and time-frequency power
        
        % Redefine confcond so it corresponds to the trial used for
        % filt/test
        
        confconds2 = confconds(trial4test);
        
        % Redefine thfiltdat so that it only uses trial4test : CAREFUL, here we extract theta phase from the midfrontal thetadata !
        thfiltdat = reshape( hilbert( filterFGx2(reshape(thetadata(:,:,trial4test),1,[]),EEG.srate,EEG.th_maxfreq,5)' )',   EEG.pnts,size(trial4test, 2)); % Hilbert decomposition of the data according to max theta freq
        thfiltdat_ang = angle(thfiltdat);
        
        % Phases
        phases = linspace(-pi,pi,51);
        
        % For theta phase-frequency power
        
        % fft of data from maxphase GED : maxphdata
        eegX = fft( reshape(maxphdata(maxphcompidx,:,:),1,nData) ,nConv);
        
        % loop over frequencies
        for fi=1:num_frex
            
            as = ifft( eegX.*cmwX(fi,:) );
            as = as(halfwave+1:end-halfwave);
            as = reshape(as,EEG.pnts,size(trial4test, 2));
            
            
            % condition-average baseline
            basepow = mean(mean( abs(as(baseidx(1):baseidx(2),:)).^2,2),1);
            
            % power matrix
            trialpow = zeros(size(as, 2), nbin);
            
            
            % loop over trials
            for triali = 1:size(as, 2)
                
                binedges = discretize(thfiltdat_ang(:,triali),linspace(-pi,pi,52));
                
                % loop over phase bins
                for bini = 1:nbin
                    
                    trialpow(triali,bini) = nanmean(abs(as(binedges==bini,triali)).^2);
                    
                end % end of phase bin loop
                
            end % end of trial loop
            
            
            for condi=1:5
                
                % phase - frequency power
                tf_phaseresolved(iteri, subno,condi,fi,:) = 10*log10( nanmean(trialpow(confconds2==condi,:),1) ./ basepow ); % confconds==condi returns a logical, tests whether confconds equals condi
                
                maxph =  maxminphases(subno,2,condi);
                
                if maxph == -pi
                    phaseidx = 1;
                elseif maxph == pi
                    phaseidx = size(phases,2);
                else
                    phaseidx = find(diff(sign(diff(abs(phases-maxph))))>0)+1; % Find the closest phase value to max corr phase
                end
                
                
                zeroph = find(phases(1,:)==0); % Find zero phase on the phase vector used for analyses
                
                bindiff = phaseidx - zeroph; % The index shift to use
                
                
                tf_phaseresolved_shifted(iteri, subno, condi, fi, :) = circshift(tf_phaseresolved(iteri, subno, condi, fi,:), - bindiff); % shifts the tfmap so that phase 0 is the phase of the subject/condition maximum power-RT correlation
                
                
            end % end condition loop
            
        end % end frequencies loop
        
        
        
        n2test = n2test+ntest;
        disp(['iteration n°' num2str(iteri)])
        
    end % end cross validation loop
    
    disp(['subject n°' num2str(subno)])
    
end % end subject loop



%% Plot and permutation

maps_maxminph_avg2 = squeeze(mean(maps_maxminph_avg, 1));
tf_phaseresolved_shifted2 = squeeze(mean(tf_phaseresolved_shifted, 1));

% Plot averaged topographical map
                 
figure(2), clf
topoplot(mean(maps_maxminph_avg2, 1), EEG.chanlocs, 'numcontour', 6, 'gridscale', 200, 'maplimits', 'maxmin')
title('Averaged across participants')

% show laplacian-transformed topoplot to compare with study 2
maps_maxminph_avg_lap = laplacian_perrinX(maps_maxminph_avg2',[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]); % transpose maps_maxminph_avg because first dimension must be channels
figure(2), clf
topoplot(mean(maps_maxminph_avg_lap', 1), EEG.chanlocs, 'numcontour', 6, 'gridscale', 200, 'maplimits', 'maxmin')
title('Averaged across participants')


% % Plot averaged TF power
% figure
% climdb  = [-1 1];
% contourf(times2save,frex_tf,squeeze(mean(mean(tfmaxph_avg, 1),2)),40,'linecolor','none')
% set(gca, 'clim', climdb, 'FontSize',8,'yscale','log','ytick',logspace(log10(min(frex_tf)),log10(max(frex_tf)),6),'yticklabel',round(logspace(log10(min(frex_tf)),log10(max(frex_tf)),6)*10)/10)
% title('TF power of the component')
% colormap jet
% colorbar
% ylabel(colorbar, 'dB')
% ylabel('Frequency')
% xlabel('Time (ms)')
% hold on 
% 


%zscore each frequency line

tf_z = zeros(28, 5, 40, 51);
for subno = 1:28
    for condi = 1:5
        for freqi = 1:size(frex,2)
            tf_z(subno,condi,freqi,:) = zscore(tf_phaseresolved_shifted2(subno, condi, freqi,:));
        end
    end
end




figure(4), clf
climdb  = [-0.3 0.3];
contourf(linspace(-pi,pi,51),frex,squeeze(nanmean(squeeze(nanmean(tf_z(:,:,:,:), 1)), 1)),40,'linecolor','none')
set(gca, 'clim', climdb, 'FontSize',8,'yscale','log','ytick',logspace(log10(min(frex)),log10(max(frex)),6),'yticklabel',round(logspace(log10(min(frex)),log10(max(frex)),6)*10)/10)
title('\theta phase-resolved power')
colormap jet
colorbar
ylabel(colorbar, 'dB  from baseline')
ylabel('Frequency')
xlabel('\theta phase (rad)')



%% Permutation testing

% NULL HYPOTHESIS WOULD BE THAT POWER IS THE SAME AT EACH PHASES, IN WHICH
% CASE PHASES CAN BE RANDOMLY SWAPPED (cut at one phase and everything
% that's after that phase is put at the beginning of the distribution)

% Get averaged across subject map

tf_z_avg = squeeze(nanmean(squeeze(nanmean(tf_z(:,:,:,:), 1)), 1));

% p-value
pval = 0.001;

% convert p-value to Z value
zval = abs(norminv(pval));

% number of permutations
n_permutes = 10000;

% initialize null hypothesis maps
permmaps = zeros(n_permutes,num_frex,size(tf_z_avg, 2));

% number of phase bins
nphase = size(tf_z_avg, 2);

% Copy tf map in second matrix for permutation

phases = linspace(-pi, pi, 51);


% generate maps under the null hypothesis
for permi = 1:n_permutes

      cutpoint = randsample (10:size(tf_z_avg, 2)-10, 1);

%     randphaseidx = dsearchn(phases', cutpoint');
    tf_z_avg_2 = [tf_z_avg(:,cutpoint:end) tf_z_avg(:,1:cutpoint-1)];
               
    % compute the "difference" map
    % what is the difference under the null hypothesis?
    permmaps(permi,:,:) = tf_z_avg_2 - tf_z_avg ;
end


% show non-corrected thresholded maps

% compute mean and standard deviation maps
mean_h0 = squeeze(mean(permmaps));
std_h0  = squeeze(std(permmaps));

% now threshold real data...
% first Z-score
zmap = (tf_z_avg-mean_h0) ./ std_h0;

% threshold image at p-value, by setting subthreshold values to 0
zmap(abs(zmap)<zval) = 0;

figure(5), clf
contourf(linspace(-pi,pi,51),frex,tf_z_avg,40,'linecolor','none')
hold on
contour(linspace(-pi,pi,51),frex, logical(zmap),1,'linecolor','k');
climdb  = [-0.3 0.3];

set(gca, 'clim', climdb, 'FontSize',8,'yscale','log','ytick',logspace(log10(min(frex)),log10(max(frex)),6),'yticklabel',round(logspace(log10(min(frex)),log10(max(frex)),6)*10)/10)
title('\theta phase-resolved power')
colormap jet
colorbar
ylabel(colorbar, 'dB  from baseline')
ylabel('Frequency')
xlabel('\theta phase (rad)')


% CONTINUE WITH CLUSTER CORRECTION


%% corrections for multiple comparisons

% initialize matrices for cluster-based correction
max_cluster_sizes = zeros(1,n_permutes);
% ... and for maximum-pixel based correction
max_val = zeros(n_permutes,2); % "2" for min/max

% loop through permutations
for permi = 1:n_permutes
    
    % take each permutation map, and transform to Z
    threshimg = squeeze(permmaps(permi,:,:));
    threshimg = (threshimg-mean_h0)./std_h0;
    
    % threshold image at p-value
    threshimg(abs(threshimg)<zval) = 0;
    
    
    % find clusters (need image processing toolbox for this!)
    islands = bwconncomp(threshimg);
    if numel(islands.PixelIdxList)>0
        
        % count sizes of clusters
        tempclustsizes = cellfun(@length,islands.PixelIdxList);
        
        % store size of biggest cluster
        max_cluster_sizes(permi) = max(tempclustsizes);
    end
    
    
    % get extreme values (smallest and largest)
    temp = sort( reshape(permmaps(permi,:,:),1,[] ));
    max_val(permi,:) = [ min(temp) max(temp) ];
    
end

%% show histograph of maximum cluster sizes

figure(3), clf
hist(max_cluster_sizes,100);
xlabel('Maximum cluster sizes'), ylabel('Number of observations')
title('Expected cluster sizes under the null hypothesis')


% find cluster threshold (need image processing toolbox for this!)
% based on p-value and null hypothesis distribution
cluster_thresh = prctile(max_cluster_sizes,100-(100*pval));

%% plots with multiple comparisons corrections

% now find clusters in the real thresholded zmap
% if they are "too small" set them to zero
islands = bwconncomp(zmap);
for i=1:islands.NumObjects
    % if real clusters are too small, remove them by setting to zero!
    if numel(islands.PixelIdxList{i}==i)<cluster_thresh
        zmap(islands.PixelIdxList{i})=0;
    end
end

figure(6), clf
contourf(linspace(-pi,pi,51),frex,tf_z_avg,40,'linecolor','none')
hold on
contour(linspace(-pi,pi,51),frex, logical(zmap),1,'linecolor','k');

set(gca, 'clim', climdb, 'FontSize',8,'yscale','log','ytick',logspace(log10(min(frex)),log10(max(frex)),6),'yticklabel',round(logspace(log10(min(frex)),log10(max(frex)),6)*10)/10)
title('\theta phase-resolved power')
colormap jet
colorbar
ylabel(colorbar, 'dB  from baseline')
ylabel('Frequency')
xlabel('\theta phase (rad)')


%% Extract power from significant cluster


% isolate biggest cluster

numPixels = cellfun(@numel,islands.PixelIdxList);
[biggest,idx] = max(numPixels);
zmap2=zmap*0;
zmap2(islands.PixelIdxList{idx}) = 1;

% check on plot
figure(8), clf
contourf(linspace(-pi,pi,51),frex,tf_z_avg,40,'linecolor','none')
hold on
contour(linspace(-pi,pi,51),frex, logical(zmap2),1, 'Linewidth', 2, 'linecolor','k');

set(gca, 'clim', climdb, 'FontSize',8,'yscale','log','ytick',logspace(log10(min(frex)),log10(max(frex)),6),'yticklabel',round(logspace(log10(min(frex)),log10(max(frex)),6)*10)/10)
title('\theta phase-resolved power')
colormap jet
colorbar
ylabel(colorbar, 'dB  from baseline')
ylabel('Frequency')
xlabel('\theta phase (rad)')



clust = logical(zmap2);
clusterpw = zeros(size(tf_z, 2), size(tf_z, 1));

for condi = 1:size(tf_z, 2)
    for subno = 1:size(tf_z, 1)
        tempdat =  squeeze(tf_z(subno, condi, :, : ));
        clusterpw(condi, subno) = mean(tempdat(clust));
    end
end

colors=[0    0.4470    0.7410; % Define colors
    0.8500    0.3250    0.0980; 
    0.9290    0.6940    0.1250; 
    0.4940    0.1840    0.5560; 
    0.4660    0.6740    0.1880];

figure(15), clf

std_dev = std(clusterpw, 1,  2);

 for condi = 1:5
     
bar(condi, mean(clusterpw (condi,:)), 'EdgeColor', 'none')
hold on
errorbar(condi, mean(clusterpw (condi,:)),std_dev(condi), 'Color', 'k' )
hold on
 end
set(gca, 'xtick', [1 2 3 4 5], 'xticklabel', {'No conflict', 'Little conflict' 'Lot of conflict', 'Mixed correct', 'Error'})
xtickangle(45)

% test condition differences with anova
[p,tbl,stats] = anova1(clusterpw')




