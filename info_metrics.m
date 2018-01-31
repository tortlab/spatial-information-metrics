function [output FRbin] = info_metrics( FR, stimulus, nbinFRSH, occupSK, nshuffles)
% INFO_METRICS calculates Shannon and Skaggs information metrics
%
% Syntax: output = info_metrics( FR, stimulus, nbinFRSH, occupSK, nshuffles)
%
% 1) FR: firing rate (ntrialsXnspacebins or nbinsX1, defining the auxiliar
%    variable 'stimulus')
% 2) stimulus: use [] if FR(ntrials,nspacebins) or same dimension as FR. 
% 3) nbinFRMI: number of bins for FR, for Shannon's mutual infomation
%     calculation.
% 4) occupSK: occupation % (or time or bins)  (same dimension as FR), for
%     Skaggs's information calculation.
% 5) nshuffles: number of shufflings for generate surrogate distribution
%     (method: shuffle on pooled spatial bins). For speed, use 0 during
%     testing.
%
% Example:
%     FR=zeros(10,100);                 %basalFR=0Hz
%     FR(:,45:55)=10; FR(:,49:51)=20;   %fieldFR=10 to 20 to 10 Hz
%     stimulus=[];                      %spatial bins are defined in FR
%     nbinFRSH=10;                      %10 FRbins
%     occupSK=ones(size(FR));           %equioccupated
%     nshuffles=100;                    %100 trialFR shufflings
%     output = info_metrics( FR, stimulus, nbinFRSH, occupSK, nshuffles)
%
% R. Pavão
% Brain Institute

if isempty(stimulus)
    spatial_bin=1:size(FR,2);
    stimulus=repmat(spatial_bin,size(FR,1),1);
end
FRtemp=FR(:); FRtemp=FRtemp(~isnan(FRtemp));
stimulus=stimulus(:); stimulus=stimulus(~isnan(FRtemp));

%if nshuffles~=0, disp('SHANNON''s Mutual Information...'); end
FRbin=NaN(size(FRtemp));
%tempedges = quantile(FRtemp,linspace(0,1,nbinFRSH+1)); tempedges(end)=tempedges(end)+0.01;
tempedges=[min(FRtemp) quantile(FRtemp(FRtemp>min(FRtemp)),linspace(0,1,nbinFRSH))]; 
tempedges(end)=tempedges(end)+0.01;


for b=1:length(tempedges)-1
    FRbin(FRtemp>=tempedges(b) & FRtemp<tempedges(b+1))=b;
end
output.ShannonMI = MutualInformation(stimulus,FRbin);

%if nshuffles~=0, disp('SKAGGS''s Information...'); end
stim_list=unique(stimulus);
averageFR=NaN(size(stim_list'));
p=NaN(size(stim_list'));
for s=1:length(stim_list)
    averageFR(s)=nanmean(FRtemp(stimulus==stim_list(s)));
    p(s)=nansum(occupSK(stimulus==stim_list(s)))/nansum(occupSK(:));
end
averageFRnorm = averageFR/nansum(p.*averageFR); %normalize across all trials
nonzero=averageFR~=0;

% bits/s
output.SkaggsSec   =    nansum(p(nonzero) .*     averageFR(nonzero)   .*   log2(averageFRnorm(nonzero)));
% bits/spike
output.SkaggsSpike =    nansum(p(nonzero) .* averageFRnorm(nonzero)   .*   log2(averageFRnorm(nonzero)));

if nshuffles>0
    %disp('Shuffling...');
    output.ShannonMI_shuffle   = NaN(1,nshuffles);
    output.SkaggsSec_shuffle   = NaN(1,nshuffles);
    output.SkaggsSpike_shuffle = NaN(1,nshuffles);
    for s=1:nshuffles
        shuffleoutput = info_metrics( shuffle(FR(:)), stimulus, nbinFRSH, occupSK, 0); %recursive
        output.ShannonMI_shuffle(s)   = shuffleoutput.ShannonMI;
        output.SkaggsSec_shuffle(s)   = shuffleoutput.SkaggsSec;
        output.SkaggsSpike_shuffle(s) = shuffleoutput.SkaggsSpike;
    end
    mu=mean(output.ShannonMI_shuffle); si=std(output.ShannonMI_shuffle);
    output.ShannonMI_p = sum(output.ShannonMI>output.ShannonMI_shuffle)/nshuffles;
    output.ShannonMI_Zsh=(output.ShannonMI-mu)/si;

    mu=mean(output.SkaggsSec_shuffle); si=std(output.SkaggsSec_shuffle);
    output.SkaggsSec_p = sum(output.SkaggsSec>output.SkaggsSec_shuffle)/nshuffles;
    output.SkaggsSec_Zsh=(output.SkaggsSec-mu)/si;
    
    mu=mean(output.SkaggsSpike_shuffle); si=std(output.SkaggsSpike_shuffle);
    output.SkaggsSpike_p = sum(output.SkaggsSpike>output.SkaggsSpike_shuffle)/nshuffles;
    output.SkaggsSpike_Zsh=(output.SkaggsSpike-mu)/si;
end

%if nshuffles~=0, disp('Done!'); end

end

