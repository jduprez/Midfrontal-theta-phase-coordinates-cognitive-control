%% Midfrontal theta phase

% This script does the 	analyses from the paper  'Midfrontal theta phase coordinates behaviorally relevant brain computations during cognitive control' (https://doi.org/10.1016/j.neuroimage.2019.116340)
% It is largely based on the equations and scripts accompanying the book 'Analyzing Neural Time Series Data' by Mike X. Cohen (mikexcohen.com).
% Requires the EEGlab toolbox for topographical plots
% Requires the FilterFGX2 function for narrow-band filtering
% Works with EEG data structures in the EEGLAB format
% This version works on the first dataset described in the paper and should be adjusted to your personal needs (especially the participant list and the GED component selection section).
% The data from the paper can be found at : data.donders.ru.nl

% Version 2 12/09/2019
% duprez.joan@gmail.com
% ## https://duprezjoan.wixsite.com/joanduprez ##



clear %close all




%% PARTICIPANT LIST


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


%% SOFT CODED CONVOLUTION PARAMETERS

min_frex  =  1;
max_frex  = 50;
num_frex  = 60;
nCycRange = [4 12];

times2save = -300:15:1200; % for temporal downsampling (at 67 Hz here) : only for plotting
basetime   = [ -500 -200 ]; % baseline period
ncond = 5; % number of conditions
nchan = 64; % number of channels
nbin = 20; % number of bins for phase-resolved analyses





%% INITIALIZE OUTPUT MATRICES

% for subject averaged topographical maps
mapsTh_avg = zeros(nchan,size(subjects,1));

% for subject averaged TF maps
tfTheta_avg=zeros(ncond,num_frex,length(times2save),size(subjects,1));

% for subject averaged peak theta power frequency
maxtheta=zeros(1,size(subjects,1));

% for subject averaged RT-power correlation
rtcorr_sub = zeros(size(subjects,1),ncond,length(times2save));

% for subject averaged phase-resolved RT-power correlation
rtcorrsPhase_sub = zeros(size(subjects,1), ncond, nbin);



% loop over subjects
for subno=2:size(subjects,1)
    
    


    %% LOAD FILE AND CONVOLUTION
        
    % load data
    
    load([  subjects{subno,1} '_ready.mat' ])

    
    EEG.data = double(EEG.data); % converts to double precision
    % apply Laplacian transform if needed
    %EEG.data(1:64,:,:) = laplacian_perrinX(EEG.data(1:64,:,:),[EEG.chanlocs(1:64).X],[EEG.chanlocs(1:64).Y],[EEG.chanlocs(1:64).Z]);
    
    % remove condition-specific ERP (bsxfun also does it in a channel specific manner)
    for condi = 1:5
        EEG.data(:,:,confconds==condi) = bsxfun(@minus,EEG.data(:,:,confconds==condi),squeeze(mean(EEG.data(:,:,confconds==condi),3))); % subtracts the ERP
    end 
    
    % setup convolution parameters
    
    frex = logspace(log10(min_frex),log10(max_frex),num_frex); % logarithmically spaced frequency
    
    % wavelet parameters
    s = logspace(log10(nCycRange(1)),log10(nCycRange(2)),num_frex)./(2*pi.*frex); % standard deviation of the Gaussian (with a logarithmic increase in the number of cycles)
    t = -2:1/EEG.srate:2;
    halfwave = floor((length(t)-1)/2);
    
    % convolution parameters
    nData = EEG.trials*EEG.pnts;
    nWave = length(t);
    nConv = nData+nWave-1;
    
    % wavelets
    cmwX = zeros(num_frex,nConv);
    for fi=1:num_frex
        cmw = fft(  exp(1i*2*pi*frex(fi).*t) .* exp( (-t.^2)/(2*s(fi)^2) )  ,nConv);
        cmwX(fi,:) = cmw ./ max(cmw); % amplitude-normalize in the frequency domain
    end



    %% GET THETA FREQUENCY
    
    % find peak midfrontal theta frequency
    
    % variables for TF

    % indices of the time-points to save
    times2saveidx = dsearchn(EEG.times',times2save'); 

    % indices for the baseline
    baseidx = dsearchn(EEG.times',basetime'); 
    
    % initialize trial-average output matrix
    tfFCz = zeros(num_frex,length(times2save));
    

    % Convolution and power extraction

    % fft of data, reshape in 1 time series of concatenated trials    
    eegX = fft( reshape(EEG.data(47,:,:),1,nData) ,nConv); % 47 is the FCz electrode

    % loop over frequencies
    for fi=1:num_frex

		% inverse fft of the multiplication of both signal and wavelets fourier power spectra
        as = ifft( eegX.*cmwX(fi,:) ); 

        % trim the result of convolution (since nConv = nData+nWave-1)
        as = as(halfwave+1:end-halfwave); 

        % reshape in time*trials
        as = reshape(as,EEG.pnts,EEG.trials); 
        
        % condition-average baseline
        basepow = mean(mean( abs(as(baseidx(1):baseidx(2),:)).^2,2),1);
        
        % Mean power (dB) for conflict conditions
        tfFCz(fi,:) = 10*log10( mean(abs(as(times2saveidx,confconds==2 | confconds==3 | confconds==4)).^2,2) ./ basepow ); 

    end % end of frequency loop

    
    % define time-frequency window to look for theta peak power
    thetatime = dsearchn(times2save',[0 800]');
    thetafreq = dsearchn(frex',[3 10]');
    
    % get the index of the peak frequency and peak time in tfFCz
    [maxfreq,maxtime] = ind2sub(size(tfFCz),find(tfFCz==max(reshape(tfFCz(thetafreq(1):thetafreq(2),thetatime(1):thetatime(2)),1,[])))); 

    maxtheta(subno) = frex(maxfreq); % store max theta power frequency
    
    


    %% GED FOR THETA

    % Indices to use in EEG.times
    times2use = times2saveidx(maxtime)-round(300*EEG.srate/1000):times2saveidx(maxtime)+round(250*EEG.srate/1000); %  times2saveidx(maxtime): index of peak theta power time
   
    tempdat   = EEG.data(1:64,:,:); % temporary data to be used for GED


    % Create the S matrix

    % Narrow-band filter (at peak theta power frequency) via frequency-domain Gaussian
    td    = filterFGx2(tempdat,EEG.srate,frex(maxfreq),3); 

	% reshape into 64 concatenated trials 
    td    = reshape( td(:,times2use,:),64,[] ); 

    % subtract the mean (DO NOT FORGET, data must have zero mean before computing the covariance matrix)
    td    = bsxfun(@minus,td,mean(td,2)); 

    % theta filtered channel covariance matrix
    covAt = (td*td')/size(td,2);  
    

    % Create the R matrix

	% reshape the broabdand data into 64 concatenated trials 
    td    = reshape( tempdat(:,times2use,:),64,[] );

    % subtract the mean (DO NOT FORGET, data must have zero mean before computing the covariance matrix)
    td    = bsxfun(@minus,td,mean(td,2));

    % reference covariance matrix (broadband, unfiltered signal)
    covBB = (td*td')/size(td,2);  

    % GED
    [evecsTh,evals] = eig( covAt,covBB ); % The columns of evecsTh are the spatial filters, the column with the largest evals 
                                          % is the one that maximally differentiate covAt and covBB
                                                         
   	
   	% sort evecs according to evals 
    [~,sidx] = sort(diag(evals));
    evecsTh  = evecsTh(:,sidx); 

    % forward model to see activation patterns
    mapsTh   = inv(evecsTh'); 
    
    % force sign of components so that topographical maps always show positive values
    for ci=1:64
        [~,idx] = max(abs(mapsTh(:,ci))); % find strongest weight
        mapsTh(:,ci) = mapsTh(:,ci) * sign(mapsTh(idx,ci)); % force to positive sign
    end
    
        
    % create data and remove components with eyeblinks
    
    % keep the 15 components with the highest eigenvalues
    comps2keep = 50:64;
    
    thetadata  = reshape( (reshape(EEG.data(1:64,:,:),64,[])'*evecsTh(:,comps2keep))' ,[length(comps2keep) EEG.pnts EEG.trials]); 
    % reshape(EEG.data(1:64,:,:),64,[]) gives a 64*(concatenated trials) matrix then we transpose it and multiply it by the 64*15 
    % eigenvector matrix (evecsTh(:,comps2keep))) which gives a (concatenated trials)*15 matrix of weighted concatenated trials. 
    % We take the transpose of that which gives a 15 (component)*(concatenated trials) matrix, which we reshape in a 15 (comp)*n(timepoints)*m(trials) matrix
    
    mapsTh = mapsTh(:,comps2keep); % overwrites mapsTh with only the components we want to keep


    % keep the compponents without eyeblinks : DEPENDS ON YOUR ELECTRODE MONTAGE
    comps2save = (mean(mapsTh([11 12 19 47 48 32 46 49 56],:)) - mean(mapsTh([1 33 34],:))) > 0; % vector of logical values
    % compares if activity in central electrodes is higher that activity at electrodes near the eyes

    % only keep components without eye blinks in the forward model and spatially transformed data matrices
    mapsTh     = mapsTh(:,comps2save);
    thetadata  = thetadata(comps2save,:,:);


    % Plot forward models to inspect activation patterns    ## UNCOMMENT TO PLOT ##

    % figure(1), clf

    % for topi = 1:size(mapsTh, 2)
       
    %      subplot(5,3,topi)
    %      topoplot(mapsTh(:,topi),EEG.chanlocs,'numcontour', 6, 'gridscale', 100)
    %      title([  'Comp ' num2str(topi)  ])
       
    %  end
    
    % By default keep last component (component with the highest eigenvalue)
    thetacompidx = size(mapsTh, 2);
%    
%     hold on
%     subplot(5,3, thetacompidx(1))
%     title([  'Comp ' num2str(thetacompidx(1))  ], 'Color', 'r') % red title for the component with highest conflict-no conflict difference
% 


% force component selection : only do after subject-specific inspection of the activation patterns to choose the theta component if the component with the highest eigenvalue is not the good one
% in that case :
% if strcmpi(subjects{subno,1},'filename'), thetacompidx= ' chosen component' ; end 
 

    mapsTh_avg(:,subno) = mapsTh(:,thetacompidx); % store forward problem in subject averaged matrix

   

    % keep the selected component
    thetadata = thetadata(thetacompidx,:,:);
    EEG.thetadata = thetadata; % store in EEG struct

    
    % keep theta peak frequency for further analysis
    EEG.th_maxfreq = frex(maxfreq);


    %% THETA: TF DECOMPOSITION TO INSPECT TIME-FREQUENCY POWER OF THE COMPONENT
    
        % fft of data
         eegX = fft( reshape(thetadata,1,nData) ,nConv);
        
        % loop over frequencies
        for fi=1:num_frex
            
            as = ifft( eegX.*cmwX(fi,:) );
            as = as(halfwave+1:end-halfwave);
            as = reshape(as,EEG.pnts,EEG.trials);
            
            % condition-average baseline
            basepow = mean(mean( abs(as(baseidx(1):baseidx(2),:)).^2,2),1);
            
            % loop over conditions
            for condi=1:5

                tfTheta_avg(condi,fi,:,subno) = 10*log10( mean(abs(as(times2saveidx,confconds==condi)).^2,2) ./ basepow );

            end % end condition loop
            
        end % end frequencies loop

   % Plot subject TF maps to inspect power    ## UNCOMMENT TO PLOT ##

%     figure(2), clf
% 
%     climdb  = [-2.5 2.5];
%     for condi=1:5
%         subplot(3,2,condi)
%         contourf(times2save,frex,squeeze(tfTheta_avg(condi,:,:,1)),40,'linecolor','none')
%         set(gca, 'clim', climdb, 'FontSize',8,'yscale','log','ytick',logspace(log10(min(frex)),log10(max(frex)),6),'yticklabel',round(logspace(log10(min(frex)),log10(max(frex)),6)*10)/10)
%         title([ '\theta component TF power for cond ' num2str(condi)])
%         colormap jet
%         colorbar
%         ylabel(colorbar, 'dB  from baseline')
%     end



    %% PHASE-RESOLVED BRAIN-BEHAVIOR RELATIONSHIP

    % Filter-Hilbert data according to max theta freq defined above     
    thfiltdat = reshape( hilbert( filterFGx2(reshape(thetadata(:,:),1,[]),EEG.srate,EEG.th_maxfreq,5)' )',   EEG.pnts,EEG.trials); 


    % time series of theta power-RT correlations 
      
          % loop over time points
          for ti=1:length(times2save)

              % loop over conditions
              for condi=1:ncond

               % Spearman correlation 
               rtcorr_sub(subno,condi,ti) = corr(abs(thfiltdat(times2saveidx(ti),confconds==condi))'.^2 ,rt(confconds==condi)','Type','Spearman','rows','complete');

           end % end of condition loop

       end % end of time points loop


    % theta power-RT correlations in phase bins

    % initialize correlation matrix (conditions * phase bins)
    rtcorrsPhase = zeros(ncond,nbin);

    % get time indices of the time window of interest (based on overall power-RT significance against 0)
    tidx = dsearchn(EEG.times',[300 1200]');


        % loop over conditions
        for condi=1:ncond       
        
        % get the phase angles time series
        thphasdat = angle( thfiltdat(tidx(1):tidx(2),confconds==condi) );

        % get the power time series
        thpowrdat = abs  ( thfiltdat(tidx(1):tidx(2),confconds==condi) ).^2;


        % get phase-specific power

        % initialize matrix containing mean power at each phase bin
        trialpow = zeros(size(thpowrdat, 2), nbin);
        
            % loop over trials
            for triali = 1:size(thpowrdat, 2)

            binedges = discretize(thphasdat(:,triali),linspace(-pi,pi,21));

                % loop over phase bins
                for bini = 1:nbin

                trialpow(triali,bini) = nanmean(thpowrdat(binedges==bini,triali));

                end % end of phase bin loop

            end % end of trial loop
        
        
            % now correlate with RT

            % loop over phase bins
            for bini=1:nbin
            
            % Spearman correlation
            rtcorrsPhase_sub(subno,condi,bini) = corr(trialpow(:,bini),rt(confconds==condi)','Type','Spearman','rows','complete'); % correlation
        
            end % end of phase bin loop

        end % end of condition 

save( [  subjects{subno,1} '_GED.mat' ] , 'accuracy', 'condlabels', 'condmarks', 'confconds', 'corrTimes', 'EEG', 'emgonset', 'respmarks', 'rt')

subno
end % end of subject loop

