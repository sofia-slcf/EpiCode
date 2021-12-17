function pnh_cluster(ipatient)

%% Analysis script for Valerio Frazzini's study on electrophysiology of pathoanatomy
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
% requires bandpassFilter.m from Mario
% requires releaseDec2015 from Neuralynx website

%% Add path

restoredefaultpath

%% Stephen's paths

if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/pnh/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/binaries'));
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/subaxis
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\pnh
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared\utilities
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\subaxis
end

ft_defaults
feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

% load settings
config = pnh_setparams;

% read muse markers
% MuseStruct{ipatient} = readMuseMarkers(config{ipatient}, false);

% align markers
% MuseStruct{ipatient} = alignMuseMarkersXcorr(config{ipatient}, MuseStruct{ipatient}, false);

% add seizures
% MuseStruct{ipatient} = updateMarkers(config{ipatient}, MuseStruct{ipatient}, {'CriseStart', 'CriseEnd'});

% % read LFP data (trial-by-trial)
% LFP{ipatient} = readLFP(config{ipatient}, MuseStruct{ipatient}, true);
% 
% % make TFR
TFR{ipatient} = TFRtrials(config{ipatient}, true);

% for i = 0:0.1:0.5
%     cfg_fig = [];
%     cfg_fig.colorscale = i;
%     Figure2(cfg_fig);
% end


%  
% % write parameters for spyking circus
% %     [MuseStruct{ipatient}] = updateBadMuseMarkers(config{ipatient}, MuseStruct{ipatient});
% %     writeSpykingCircusDeadfiles(config{ipatient}, MuseStruct_aligned{ipatient}, true);
% %     writeSpykingCircusParameters_new(config{ipatient});
% %     [filelist, sampleinfo, timestamps, hdr] = writeSpykingCircusFileList(config{ipatient}, true);
% %     [~, samples_separate, ~, ~] = writeSpykingCircusFileList(config{ipatient}, false);
% 
% % read spike data from Phy as one continuous trial
% SpikeRaw{ipatient} = readSpikeRaw_Phy_new(config{ipatient}, true);
% 
% % segment into trials based on IED markers
% SpikeTrials_timelocked{ipatient}    = readSpikeTrials_MuseMarkers(config{ipatient}, MuseStruct{ipatient}, SpikeRaw{ipatient}, true);
% 
% % get spike density and stats vs. baseline
% SpikeDensity_timelocked{ipatient}   = spikeTrialDensity(config{ipatient}, SpikeTrials_timelocked{ipatient}, true);
% 
% % plot intervals between IEDs
% intervals{ipatient} = inter_trial_intervals(config{ipatient}, MuseStruct{ipatient}, true);

% % plot LFP timecourse examples for article
% config{ipatient}.plot.name      = {'patternA','patternB','patternC'};
% config{ipatient}.plot.ncols     = 4;
% config{ipatient}.plot.postfix   = "_IEDs";
% plotTimeCoursesExamples(config{ipatient});
% 
% config{ipatient}.plot.unit{1}   = ones(size(config{ipatient}.plot.unit{1})) * -1; % all individually
% 
% if ipatient == 2
%     config{ipatient}.plot.name      = {'FA','ES'};
% else
%     config{ipatient}.plot.name      = {'PSW','FA','ES'};
% end
% 
% config{ipatient}.plot.ncols     = 4;
% config{ipatient}.plot.postfix   = "_IEDs";
% plot_patterns_multilevel_examples2(config{ipatient});
% 
% % plot only seizure timecourse examples for article
% config{ipatient}.plot.name      = "seizureA";
% config{ipatient}.plot.ncols     = 1;
% config{ipatient}.plot.postfix   = "_seizureA";
% plotTimeCoursesExamples(config{ipatient});
% 
% config{ipatient}.plot.unit{1}   = ones(size(config{ipatient}.plot.unit{1})) * 0; % all individually
% config{ipatient}.plot.name      = "SEIZURE";
% config{ipatient}.plot.ncols     = 1;
% config{ipatient}.plot.postfix   = "_SEIZURE";
% plot_patterns_multilevel_examples2(config{ipatient});
% 
% % plot firing change for each pattern
% plotstats(config{ipatient});
