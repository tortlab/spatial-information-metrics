function [fr, occup, spatial_metrics] = spatial_information(x,x_srate, spk, trials, varargin)
% SPATIAL_INFORMATION computes the spatial information metrics as described
% in Souza et al. 2018, On information metrics for spatial coding.
% Neuroscience (in press).
%
%
% [FR, occup, spatial_metrics] = spatial_information(X,X_SRATE,SPK,TRIALS) 
%    where X is the animal position in the linear track, X_SRATE is the
%    sampling rate of the position, SPK is a N-by-1 cell array where SPK{1}
%    has a vector with the spiketimes (in sec) of neuron 1 and TRIALS is a
%    T-by-2 matrix where TRIALS(I,1) and TRIALS(I,2) are the position
%    indexes in X where the Ith trial begins and ends, respectively. 
% 
%    The function returns
%       FR: a T-by-NBIN-by-N matrix with the firing rates of the N neurons
%          over each trial T and NBIN spatial bins.
%       OCCUP: a T-by-NBIN matrix with the occupancy of the animal.
%       SPATIAL_METRICS: a struct with the six spatial metrics (MI, Isec,
%           Ispike, NormMI, NormIsec and NormIspike) computed for each
%           neuron.
% 
%   Example of use:
%        load linear_track_sample
%        [fr, occup, sm] = spatial_information(x,x_srate,spk,trials);
%        subplot(2,2,1)
%        imagesc(squeeze(fr(:,:,1)))
%        xlabel('Space bins'); ylabel('Trials'); title('Firing rate of neuron 1')
%        subplot(2,2,2)
%        imagesc(occup)
%        xlabel('Space bins'); ylabel('Trials'); title('Animal occupancy')
%        subplot(2,2,3)
%        bar([sm.MI;  sm.Isec;  sm.Ispike;]')
%        xlabel('Neuron')
%        ylabel('Metric')
%        legend({'MI', 'Isec', 'Ispike'})
%        subplot(2,2,4)
%        bar([ sm.NormMI;  sm.NormIsec; sm.NormIspike]')
%        xlabel('Neuron')
%        ylabel('Metric')
%        legend({'Norm. MI', 'Norm. Isec', 'Norm. Ispike'})
% 
% 
% [...] = spatial_information(..., NBIN, NBINFR, NSHUFFLES) also defines
%    the number of spatial bins used NBIN (default is 25), the number
%    NBINFR of firing rate bins used to compute the MI (default is 8) and
%    the number NSHUFFLES of shufflings done to compute the normalized
%    metrics (default is 100)
% 
% B. C. Souza,
% Brain Institute, Natal, Brazil,
% January, 2018.



nbin=25;
if length(varargin)>4
    nbin = varargin{5}; 
end

nbinFR=8;
if length(varargin)>5
    nbinFR = varargin{6}; 
end

nshuffles=100;
if length(varargin)>6
    nshuffles= varargin{7}; 
end


% idxtrial, referring to x (position) and spikecount
idxtrial={};
for itrial=1:length(trials)
    idxtrial{itrial,1} = trials(itrial,1) : trials(itrial,2) ;
end

bins= linspace(min(x),max(x),nbin+1); %%%%%%%%%%%%%%%%%%%%%%%%%% 101 bins
bins(end)=bins(end)*1.0001;

xbin = NaN(size(x));
for ibin=1:length(bins)-1
    xbin(x>=bins(ibin) & x<bins(ibin+1)) = ibin;
end


x_time = (1:length(x))/x_srate;
spkcount=zeros(length(spk), length(x_time));
for icell=1:length(spk)
    spkcount(icell,:) = histc(spk{icell},x_time);
end
clear spkcountT xbinT
% spkcount, x and x_time from each Trial

spkcountT = cell(length(trials), length(spk));
xbinT = cell(length(trials));
for itrial=1:length(trials)
    for icell=1:size(spkcount,1)
        spkcountT{itrial,icell}= spkcount(icell,idxtrial{itrial,1})';
    end
%     xT{itrial,1}= x(idxtrial{itrial,1});
    xbinT{itrial,1}= xbin(idxtrial{itrial,1});
%     xtimeT{itrial,1}= x_time(idxtrial{itrial,1})';
end


% FR(#trials,#spacebins,#neuron) and occupSK(same dim1 and dim2 as FR)

fr=NaN( length(trials) , length(bins)-1 , size(spkcount,1) );
occup=fr(:,:,1);
for itrial=1:length(trials)
    for ibin=(1:length(bins)-1)
        for icell=1:size(spkcount,1)
            aux = spkcountT{itrial,icell}(xbinT{itrial,1}==ibin);
            if icell==1
                occup(itrial, ibin) = (length(aux)/x_srate);
            end
            fr(itrial, ibin, icell) = sum(aux) / occup(itrial, ibin);
        end
    end
end

MI=NaN(size(spkcount,1),1);
Isec=NaN(size(spkcount,1),1);
Ispike=NaN(size(spkcount,1),1);
NormMI=NaN(size(spkcount,1),1);
NormIsec=NaN(size(spkcount,1),1);
NormIspike=NaN(size(spkcount,1),1);


ShannonMI_p=NaN(size(spkcount,1),1);
SkaggsSec_p=NaN(size(spkcount,1),1);
SkaggsSpike_p=NaN(size(spkcount,1),1);

clear spatial_metrics

for icell=1:size(spkcount,1)
    output = info_metrics( fr(:,:,icell), [], nbinFR, occup, nshuffles);
    spatial_metrics.MI(icell)=output.ShannonMI;
    spatial_metrics.Isec(icell)=output.SkaggsSec;
    spatial_metrics.Ispike(icell)=output.SkaggsSpike;
    spatial_metrics.NormMI(icell)=output.ShannonMI_Zsh;
    spatial_metrics.NormIsec(icell)=output.SkaggsSec_Zsh;
    spatial_metrics.NormIspike(icell)=output.SkaggsSpike_Zsh;
end

end