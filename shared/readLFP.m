function [LFP] = readLFP(cfg,MuseStruct,force,savedat)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [dat_micro, dat_macro] = readLFP(cfg,MuseStruct_micro,MuseStruct_macro,force,savedat)
%
% Reads data from macro and micro electrodes, epoched according to markers extracted from Muse,
% and downsampled to same samplerate.
%
% Necessary fields:





% Note:
% Names of markers that contain a space (' ') or minus ('-') will be
% replaced by an underscore ('_').
%
% Dependencies: writeMuseMarkers.m, dir2.m, recent FieldTrip version
%
% (c) Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




fname_out = fullfile(cfg.datasavedir,[cfg.prefix,'data_aligned.mat']);

if exist(fname_out,'file') && force == false
    fprintf('************************************\n');
    fprintf('*** Loading precomputed LFP data ***\n');
    fprintf('************************************\n\n');
    
    load(fname_out,'dat_micro','dat_macro');
    
else
    fprintf('********************************\n');
    fprintf('*** (re-) computing LFP data ***\n');
    fprintf('********************************\n\n');
    
    
    % select those markers to load
    markerlist = [];
    for i = 1 : size(cfg.name,2)
        if ismember(cfg.name{i},cfg.name)
            markerlist = [markerlist, i];
        end
    end
    
    for ipart = 1:length(MuseStruct)
        
        for imarker = markerlist
            
            hasmarker = false(length(MuseStruct{ipart}),1);
            
            for idir = 1:length(MuseStruct{ipart})
                if isfield(MuseStruct{ipart}{idir},'markers')
                    if isfield(MuseStruct{ipart}{idir}.markers,(cfg.muse.startend{imarker,1}))
                        if isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}),'synctime')
                            
%                             % select MICRO files
%                             micro_filenrs = [];
%                             for ifile = 1 : size(MuseStruct{ipart}{idir}.filenames,2)
%                                 for ilabel = 1 : size(cfg.labels.micro,2)
%                                     if ~isempty(strfind(MuseStruct{ipart}{idir}.filenames{ifile},cfg.labels.micro{ilabel}))
%                                         filenrs       = [micro_filenrs, ifile];
% %                                         microlabel{ifile}   = cfg.labels.micro{ilabel};
%                                     end
%                                 end
%                             end
                            
                            for ifile = 1 : size(cfg.LFP.channel,2)
                                
                                % to deal with missing data
                                temp                    = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.LFP.channel{ifile},'*.ncs']));
                                fname{1}                = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},temp.name);
                                dat                     = ft_read_neuralynx_interp(fname);
                                
                                % filter with FT
                                cfgtemp                 = [];
                                cfgtemp.hpfilter        = cfg.LFP.hpfilter;
                                cfgtemp.hpfreq          = cfg.LFP.hpfreq;
                                dat_filt                = ft_preprocessing(cfgtemp,dat);
                                clear dat
                                
                                % create Fieldtrip trl
                                hdr                     = ft_read_header(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},temp.name));
                                Startsample             = [];
                                Endsample               = [];
                                Stage                   = [];
                                trialnr                 = [];
                                
                                for ievent = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime,2)
                                    
                                    ss  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(ievent) * hdr.Fs);
                                    idx = find(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,2}).synctime * hdr.Fs >= ss,1,'first');
                                    es  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,2}).synctime(idx) * hdr.Fs);
                                    
                                    if ~isempty(es) % && (es - ss) * hdr_micro.Fs < 4 %% find a way to check for Paul's data
                                        
                                        Startsample(ievent) = ss + cfg.epoch.toi{imarker}(1) * hdr.Fs - cfg.epoch.pad(imarker) * hdr.Fs;
                                        Endsample(ievent)   = es + cfg.epoch.toi{imarker}(2) * hdr.Fs + cfg.epoch.pad(imarker) * hdr.Fs;
                                        trialnr(ievent)     = ievent;
                                        Stage(ievent)       = -1;
                                        
                                        % find overlap with hypnogram markers
                                        for hyplabel = {'PHASE_1','PHASE_2','PHASE_3','REM','AWAKE','NO_SCORE'}
                                            if isfield(MuseStruct{ipart}{idir}.markers,[cell2mat(hyplabel),'__START__'])
                                                for i = 1 : size(MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).synctime,2)
                                                    y1 = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).synctime(i) * hdr.Fs;
                                                    y2 = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__END__']).synctime(i) * hdr.Fs;
                                                    if (y1 < ss) && (ss < y2)
                                                        fprintf('Found "%s" overlapping with "%s" : adding to trialinfo: ',cfg.name{imarker},cell2mat(hyplabel));
                                                        switch cell2mat(hyplabel)
                                                            case 'PHASE_1'
                                                                fprintf('%d\n',1);
                                                                Stage(ievent) = 1;
                                                            case 'PHASE_2'
                                                                fprintf('%d\n',2);
                                                                Stage(ievent) = 2;
                                                            case 'PHASE_3'
                                                                fprintf('%d\n',3);
                                                                Stage(ievent) = 3;
                                                            case 'REM'
                                                                fprintf('%d\n',4);
                                                                Stage(ievent) = 4;
                                                            case 'AWAKE'
                                                                fprintf('%d\n',0);
                                                                Stage(ievent) = 0;
                                                            case 'NO_SCORE'
                                                                fprintf('%d\n',0);
                                                                Stage(ievent) = 0;
                                                            otherwise
                                                                error('Unexpected label name in Hypnogram\n');
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end % ~isempty(es)
                                end
                                Offset                          = ones(size(Endsample)) * (cfg.epoch.toi{imarker}(1) - cfg.epoch.pad(imarker)) * hdr.Fs;
                                cfgtemp                         = [];
                                cfgtemp.trl                     = round([Startsample; Endsample; Offset]');
                                cfgtemp.trl(:,4)                = 1:size(cfgtemp.trl,1); 
                                cfgtemp.trl(:,5)                = trialnr;
                                cfgtemp.trl(:,6)                = idir; 
                                cfgtemp.trl(:,7)                = Stage;
                                cfgtemp.trl                     = cfgtemp.trl(Startsample > 0 & Endsample < hdr.nSamples,:); % so not to read before BOF or after EOFs
                                
                                cfgtemp.dataset                 = fname{1};
                                filedat{ifile}                  = ft_redefinetrial(cfgtemp,dat_filt);
                                clear dat_filt
                                
                                % downsample data and correct baseline
                                cfgtemp                         = [];
                                cfgtemp.resamplefs              = cfg.LFP.resamplefs;
                                if strcmp(cfg.LFP.baseline,'no')
                                    cfgtemp.demean              = 'no';
                                else
                                    cfgtemp.demean              = 'yes';
                                    cfgtemp.baselinewindow      = cfg.LFP.baselinewindow{imarker};
                                end
                                filedat{ifile}                  = ft_resampledata(cfgtemp,filedat{ifile});
                                
                                % same label over files
                                filedat{ifile}.label{1}         = cfg.LFP.channel{ifile};
                                
                                % flag for averaging
                                hasmarker(idir)                 = true;
                                
                                %                                     catch
                                %                                         fprintf('problems with file %s\n',fname{1});
                                %                                         hasdata_micro(ifile) = false;
                                %                                 end
                                
                            end

                            % concatinate channels
                            cfgtemp                             = [];
                            cfgtemp.keepsampleinfo              = 'no';
                            dirdat{idir}                        = ft_appenddata(cfgtemp,filedat{:});
                            clear filedat*
                        end
                    end
                end
            end % idir
            
            % concatinate data of different datasets (over trials)
            LFP{ipart}{imarker} = ft_appenddata([],dirdat{find(hasmarker_macro)});
            clear dirdat*
            
            % add samplerate
            LFP{ipart}{imarker}.fsample = cfg.LFP.resamplefs;
            
        end % imarker
        
    end % ipart
end

if savedat
    save(fname_out,'LFP','-v7.3');
end
