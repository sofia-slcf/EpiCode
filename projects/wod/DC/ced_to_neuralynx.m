

try %en local
    scriptpath = matlab.desktop.editor.getActiveFilename;
catch %cluster
    scriptpath = mfilename('fullpath');
end

epicodepath = [fileparts(fileparts(fileparts(scriptpath))), filesep];

addpath (genpath([epicodepath,'development']))
addpath (genpath([epicodepath,'shared']))
addpath (genpath([epicodepath,'external']))
addpath (genpath([epicodepath,'templates']))
addpath (genpath([epicodepath,'projects', filesep, 'wod']))
addpath (genpath([epicodepath,'projects', filesep, 'dtx']))
addpath (genpath([epicodepath,'projects', filesep, 'wod',filesep,'wod_functions']))
addpath

if ispc
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
    addpath \\lexport\iss01.charpier\echanges\scripts-paul\Spike2_vers_MATLAB
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
    addpath /network/lustre/iss01/charpier/echanges/scripts-paul/Spike2_vers_MATLAB

end

ft_defaults

config = eval(configscript);


chanlist= ["DC", "Vm"]; %chanlist double DC
chanlist= ["Vm","EEG-S1-L","Im"]; %chanlist Intra

for iprot= 13:15%size(config,2)
    datapath = char(fullfile(config{iprot}.rawdir,config{iprot}.directorylist{1}));
    
    CEDStruct = readCEDevents(datapath);
    
    for channame = chanlist
        data = readCEDwaveforms(datapath, channame);
        
        
        [folder, file, extension] = fileparts(datapath);
        Neurlynx_datapath= fullfile(folder,'Neuralynx_data');
        
        %Resample channels to have same samplerate for both
        cfgtemp= [];
        cfgtemp.resamplefs      = config{iprot}.DC.resamplefs;
        cfgtemp.detrend         = 'no';
        data                    = ft_resampledata(cfgtemp,data);
        
        
        
        header.filename         = fullfile(Neurlynx_datapath, [file,'_', char(channame), '.ncs']);
        header.date             = CEDStruct.starttime;
        header.dateEnd          = CEDStruct.endtime;
        header.nChans           = 1;
        header.chanunit         = data.chanunit{1};
        header.chantype         = 'waveform';
        header.label            = data.label;
        header.Fs               = data.fsample;
        header.FirstTimeStamp   = 0;
        header.TimeStampPerSample = 1;
        
        
        if ~isfolder(Neurlynx_datapath)
            mkdir(Neurlynx_datapath);
        end
        
        
        fprintf('\nWriting data in Neuralynx format with Fieldtrip : %s \n',header.filename);
        ft_write_data(header.filename,data.trial{1},'header',header, 'dataformat','neuralynx_ncs');
        
        
        
    end %channame
end %iprot






%% Save data as MATLAB structures and test filtering for each rat and channel

for iprot= 1:size(config,2)
    datapath = char(fullfile(config{iprot}.rawdir,config{iprot}.directorylist{1}));
    [folder, file, extension] = fileparts(datapath);
    Neurlynx_datapath= fullfile(folder,'Neuralynx_data');
    CEDStruct = readCEDevents(datapath);
    
    if iprot >= 13
        chanlist= "Vm";
    end
    for channame= chanlist
        
        datapath_raw = fullfile(Neurlynx_datapath,sprintf('%s_%s.%s',config{iprot}.directorylist{1}{1},channame,'ncs'));
        datapath_filt = fullfile(Neurlynx_datapath,sprintf('%s_%s_%s.%s','filt',config{iprot}.directorylist{1}{1},channame,'ncs'));
        
        cfgtemp = [];
        cfgtemp.dataset = datapath_raw;
        cfgtemp.demean = 'yes';
        cfgtemp.baselinewindow = [CEDStruct.markers.VentOff.synctime CEDStruct.markers.VentOff.synctime+10];
        data_raw = ft_preprocessing(cfgtemp);
        
        cfgtemp = [];
        cfgtemp.dataset = datapath_filt;
        data_filt = ft_preprocessing(cfgtemp);
        
        %save as matlab structures
        mat_datapath= fullfile(config{iprot}.datasavedir,'Matlab_struct');
        
        if ~isfolder(mat_datapath)
            mkdir(mat_datapath);
        end
        
        %cut data
        t1=(CEDStruct.markers.VentOff.synctime-5);
        t2=(CEDStruct.markers.VentOn.synctime+100);
        
        cfgtemp=[];
        cfgtemp.trials='all';
        cfgtemp.channel='all';
        cfgtemp.latency= [t1 t2];
        data_raw= ft_selectdata(cfgtemp,data_raw);
        
        cfgtemp=[];
        cfgtemp.trials='all';
        cfgtemp.channel='all';
        cfgtemp.latency= [t1 t2];
        data_filt= ft_selectdata(cfgtemp,data_filt);
        
        clear t1 t2
        
        %Filter data with lpfilter
        cfgtemp = [];
        cfgtemp.lpfilter= 'yes';
        cfgtemp.lpfreq=config{iprot}.DC.lpfilter ;
        cfgtemp.lpfilttype = 'but';
        cfgtemp.lpinstabilityfix = 'reduce';
        data_raw = ft_preprocessing(cfgtemp, data_raw);
        
        cfgtemp= [];
        cfgtemp.lpfilter= 'yes';
        cfgtemp.lpfreq= config{iprot}.DC.lpfilter;
        cfgtemp.lpfilttype = 'but';
        cfgtemp.lpinstabilityfix = 'reduce';
        data_filt = ft_preprocessing(cfgtemp, data_filt);
     
        save(fullfile(mat_datapath,sprintf('%s_%s.%s',config{iprot}.directorylist{1}{1},channame,'mat')),'data_raw');
        save(fullfile(mat_datapath,sprintf('%s_%s_%s.%s',config{iprot}.directorylist{1}{1},channame,'filt','mat')),'data_filt');
        save(fullfile(mat_datapath,sprintf('%s_%s_%s.%s',config{iprot}.directorylist{1}{1},channame,'events','mat')),'CEDStruct');
        
    end %channame
end %iprot

