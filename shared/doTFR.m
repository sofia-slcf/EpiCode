function [TFR] = doTFR(cfg,MuseStruct_micro,MuseStruct_macro,force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [TFR_micro, TFR_macro] = doTFR(cfg,MuseStruct_macro_micro,MuseStruct_macro_macro,dat_micro,dat_macro,force)



% Note:
% Names of markers that contain a space (' ') or minus ('-') will be
% replaced by an underscore ('_').
%
% Dependencies: writeMuseMarkers.m, dir2.m, recent FieldTrip version
%
% (c) Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




fname_out = fullfile(cfg.datasavedir,[cfg.prefix,'TFR.mat']);

if exist(fname_out,'file') && force == false
    fprintf('************************************\n');
    fprintf('*** Loading precomputed TFR data ***\n');
    fprintf('************************************\n\n');
    
    load(fname_out,'TFR');
    
else
    fprintf('********************************\n');
    fprintf('*** (re-) computing TFR data ***\n');
    fprintf('********************************\n\n');
    
    
    for ipart = 1 : size(MuseStruct_macro,2)
        
        % just define MuseStruct_macro here for only one part to simplify code
        % and similarity with no-parts
        fprintf('\n*** Starting on part %d ***\n',ipart)
        
        % process channels separately
        chan_counter = 1;   
        clear chandat
        
        for ichan = 1 : size(cfg.TFR.channel,2)
            
            clear dirdat
            
            % loop over all directories (time), concatinating channel
            for idir = 1 : size(MuseStruct_macro{ipart},2)
                
                % find corresponding file
                d                         = dir2(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.TFR.channel{ichan},'*.ncs']));
                
                % save filenames to cfg (output of function)
                cfg.fnames_ncs_TFR{ipart}{idir,ichan} = fullfile(d.folder,d.name);
                
                % load data
                cfgtemp                   = [];
                cfgtemp.dataset           = fullfile(d.folder,d.name);
                fprintf('LOADING: %s\n',cfgtemp.dataset);
                clear fname
                fname{1}                  = cfgtemp.dataset;
                dirdat{idir}              = ft_read_neuralynx_interp(fname);
                
                % downsample data
                cfgtemp                   = [];
                cfgtemp.resamplefs        = 100;
                dirdat{idir}              = ft_resampledata(cfgtemp,dirdat{idir});
                
                % truncate label to make them equal over files
                dirdat{idir}.label{1}     = dirdat{idir}.label{1}(end-6:end); % can be replaced by circus.channel
                
                % save sampleinfo
                hdr  = ft_read_header(fullfile(d.folder,d.name));
                cfg.sampleinfo_TFR{ipart}{ichan}(idir,:) = [1 hdr.nSamples];
                
            end % idir
            
            
            % concatinate data over files
            chandat{chan_counter} = dirdat{1};
            for idir = 2 : length(MuseStruct_macro{ipart})
                fprintf('Concatinating directory %d, channel %d\n',idir, ichan);
                chandat{chan_counter}.trial{1} = [chandat{chan_counter}.trial{1} dirdat{idir}.trial{1}];
                chandat{chan_counter}.time{1}  = [chandat{chan_counter}.time{1} (dirdat{idir}.time{1} + chandat{chan_counter}.time{1}(end))];
            end
            
            
            
            %
            %     % create filename for concatinated data
            %     temp        = dir(fullfile(MuseStruct_macro{ipart}{idir}.directory,['*',cfg.TFR.channel{ichan},'.ncs']));
            %     hdrtemp     = ft_read_header(fullfile(MuseStruct_macro{ipart}{idir}.directory, temp.name));
            %     clear fname
            %     subjdir     = cfg.prefix(1:end-1);
            %     partdir     = ['p',num2str(ipart)];
            %     filename    = [cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.labels.micro{ichan},'.ncs'];
            %     fname       = fullfile(cfg.datasavedir,subjdir,partdir,filename);
            
            
            chan_counter = chan_counter + 1;
        end % ichan
        
        dat{ipart} = ft_appenddata([],chandat{:});
        
        % TFR
        % time frequency analysis
        cfgtemp                         = [];
        cfgtemp.channel                 = 'all'; %ichannel;
        cfgtemp.method                  = 'mtmconvol';
        cfgtemp.output                  = 'pow';
        cfgtemp.taper                   = 'hanning';
        cfgtemp.pad                     = 'nextpow2';
        cfgtemp.keeptrials              = 'yes';
        cfgtemp.foi                     = 1:0.1:30;
        cfgtemp.t_ftimwin               = ones(size(cfgtemp.foi))*60;
        
        cfgtemp.toi                     = 0:30:dat{ipart}.time{1}(end);
        TFR{ipart}                      = ft_freqanalysis(cfgtemp,dat{ipart});
        
        
    end % ipart
        save(fname_out,'TFR','-v7.3');
end % if file already exists

