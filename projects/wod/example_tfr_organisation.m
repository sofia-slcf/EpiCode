cfg = config{irat};

analysis_names = {'timefreq_wod', 'timefreq_wod_timenorm', 'timefreq_baseline','timefreq_wod_blcorrected', 'timefreq_wod_timenorm_blcorrected', 'timefreq_baseline_blcorrected','log_timefreq_wod', 'log_timefreq_wod_timenorm', 'log_timefreq_baseline','log_timefreq_wod_blcorrected', 'log_timefreq_wod_timenorm_blcorrected','log_timefreq_baseline_blcorrected'};
idata = 1;
data_temp = load(fullfile(cfg.datasavedir,[cfg.prefix,analysis_names{idata},'.mat']));

itrial = 1;

for ifield = string(fieldnames(data_temp.timefreq_wod{itrial}))'
    if isempty(data_temp.timefreq_wod{itrial}.(ifield))
        data_temp.timefreq_wod{itrial} = rmfield(data_temp.timefreq_wod{itrial}, ifield);
    end
end

temp =  struct2cell(data_temp.timefreq_wod{itrial});

cfgtemp = [];
cfgtemp.parameter  = 'powspctrm';
tfr_all = ft_appendfreq(cfgtemp, temp{:});

cfgtemp = [];
cfgtemp.channel = {'E18LFP'};
cfgtemp.channel = {'*bLFP'};
cfgtemp.channel = {'all', '-E30LFP'};
cfgtemp.latency = [0 10];
data = ft_selectdata(cfgtemp, LFP{1}.WoD_short);

cfgtemp = [];
cfgtemp.frequency = [0 10]; %Hz
cfgtemp.latency = [0 30]; %s
cfgtemp.channel = {'*bLFP'};
data = ft_selectdata(cfgtemp, tfr_all);

% average over chan
data_avg = mean(tfr_all.powspctrm, 1, 'omitnan');
data_avg = permute(data_avg, [2 3 1]);
plot(data_avg);

%FFT (sans dimension temporelle)
data_fft = mean(tfr_all.powspctrm, 3, 'omitnan');
plot(data_fft);