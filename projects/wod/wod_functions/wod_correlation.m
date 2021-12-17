function wod_correlation(config,rongeur,per_rat) 
 
analysis_names = {'timefreq_wod', 'timefreq_wod_timenorm', 'timefreq_baseline','timefreq_wod_blcorrected', 'timefreq_wod_timenorm_blcorrected', 'timefreq_baseline_blcorrected','log_timefreq_wod', 'log_timefreq_wod_timenorm', 'log_timefreq_baseline','log_timefreq_wod_blcorrected', 'log_timefreq_wod_timenorm_blcorrected','log_timefreq_baseline_blcorrected'};
S = {'CC','HPC', 'NC', 'PTA', 'S1', 'TH'};

if ispc
    savedir = '\\lexport\iss01.charpier\analyses\wod\Sofia\images';
elseif isunix
    savedir = '/network/lustre/iss01/charpier/analyses/wod/Sofia/images';
end


for idata = 1%:length(analysis_names)
    for nrat=1:length(config)
        if nrat==4
            continue
        end
        
        chan_selec={};
        
        for icat=1:size(config{nrat}.classe.structure,2)
           
            for itrial=1:size(rongeur(nrat).catego(icat).wod_time_ini,2)  %size(config{nrat}.directorylist{1},2) marche pas
                
                if config{nrat}.classe.structure(icat).trial >1
                    continue
                end
                
                rongeur(nrat).catego(icat).wod_time_ini;
                rongeur(nrat).catego(icat).chan_depth;
              
                cfgtemp.channel   = rongeur(nrat).catego(icat).chan_ini{1,1};
                chan_selec{icat}  =  cfgtemp.channel 
            end 
        end
    end
               


 cfgtemp=[];
 cfgtemp.channel=chan_selec;
    hasjack   = 0; %or 1 specifying whether the Repetitions represent leave-one-out samples
    complex   = 'abs';%, 'angle', 'real', 'imag', 'complex', 'logabs' for post-processing of coherency
    feedback  = 'text';%'none', 'text', 'textbar' type of feedback showing progress of computation
    dimord    =[1 2 3];%specifying how the input matrix should be interpreted
%   powindx   = required if the input data contain linearly indexed
%               channel pairs. should be an Nx2 matrix indexing on each
%               row for the respective channel pair the indices of the
%               corresponding auto-spectra
   pownorm   = true %flag that specifies whether normalisation with the
%               product of the power should be performed (thus should
%               be true when correlation/coherence is requested, and
%               false when covariance or cross-spectral density is
%               requested).
 [c, v, n] = ft_connectivity_corr([40 cfgtemp.channel 25], 'dimord',[1 2 3] )
 %voila
 x=1
 
    
end