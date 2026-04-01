% O-information analysis of fMEG
% Joel Frohlich
% 25.03.24
%
% See Varley et al. 2024 BioRxiv for O-information background
% https://www.biorxiv.org/content/10.1101/2024.02.23.581375v1.full.pdf

%% Compute Omega (X) for MEG with Fairy processing algorithm
clearvars
close all
dbstop if error

folder_OP = dir('O_information/SampleData_OP');
folder_OP = folder_OP(4:3:end); % every third file is of interest (mMCG2-fMCG2.ds), sorry for hard coding
seed = round(rand(1,length(folder_OP))*1000);
fairy = dir('O_information/SampleData_FAIRY');
load('megchanlocs.mat');
nch = 10;
fs = 610.3516;

Gage = nan(1,length(folder_OP));
Omega = nan(1,length(folder_OP));
brain = nan(length(folder_OP),nch);
JH = nan(1,length(folder_OP));
CSERbb = nan(1,length(folder_OP));
CSERnb = nan(1,length(folder_OP));

% Draw random subsamples and compare the results from neural channels to these random samplings

Nrs = 0; % N resamples
Omegars = nan(length(folder_OP),Nrs);
JHrs = nan(length(folder_OP),Nrs);
CSERrs = nan(length(folder_OP),Nrs);
limit = 100;
NDX = nan(Nrs,nch,length(folder_OP));

try
    parpool
catch
    fprintf('Seems like parallel pool is already open\n')
end

%% Now the show begins ...
start = 1;
parfor idir = start:length(folder_OP)
    if ~contains(folder_OP(idir).name,'mMCG2.ds'), continue, end
    fprintf('loading %ith file...\n',idir)
    loadme = sprintf('./O_information/SampleData_OP/%s/%s.meg4',...
        folder_OP(idir).name,folder_OP(idir).name(1:end-3));
    MEG = pop_fileio(loadme);

    % get relevant channels
    pickme = [];
    mylbl = cell(1,1);
    for ich = 1:length(MEG.chanlocs)
        if any(ismember(Tmeg.ChannelName,MEG.chanlocs(ich).labels))
            pickme = [pickme ich];
            mylbl{length(pickme)} = MEG.chanlocs(ich).labels;
        end
    end

    dat = MEG.data(pickme,:,:);
    if size(dat,3) > 1, dat = reshape(dat, size(dat,1),size(dat,2)*size(dat,3)); end
    idx = strfind(folder_OP(idir).name,'SSW');
    ga = str2double(folder_OP(idir).name(idx-2:idx-1));
    Gage(idir) = ga;
    next_idir = false;  % Flag for passing control to next iteration of idir

    % get the neural channel indices
    fidx = idir + 3;
    tmp = folder_OP(idir).name;
    assert(contains(fairy(fidx).name,tmp(1:end-15)),'Cannot find the matching fairy file')
    pth = sprintf('./O_information/SampleData_FAIRY/%s/',fairy(fidx).name);
    fry = load(sprintf('%s3_fMEG/initial_brain_sensors.mat',pth));

    if size(dat,1) == 154
        T = Tmeg;
        b = fry.initial_brain_sensors(:,1);
    elseif size(dat,1) == 153 % if this is true, then probably channels MRP3 is missing
        T = Tmeg153;
        assert(contains(setdiff(Tmeg.ChannelName,mylbl),"MRP3"),'Wrong assumption')
        % Comment this out because somehow the initial_brain_sensors.mat
        % file implicitly contains the MRJ3 channel that is missing from
        % the .meg3 file

%        C = find(Tmeg.ChannelName == "MRJ3");
%         b = nan(1,nch);
%         for ich = 1:nch
%             if fry.initial_brain_sensors(ich,1) < C
%                 b(ich) = fry.initial_brain_sensors(ich,1);
%             else 
%                 b(ich) = fry.initial_brain_sensors(ich,1)+1;
%             end
%         end
%         brain(idir,:) = b;
        b = fry.initial_brain_sensors(:,1);
    elseif size(dat,1) == 152 % if this is true, then probably channels MLF3 and MRJ3 are missing
        T = Tmeg152;
        C1 = find(Tmeg.ChannelName == "MLF3");
        C2 = find(Tmeg.ChannelName == "MRJ3");
        assert(any(ismember(setdiff(Tmeg.ChannelName,mylbl),"MLF3")),'Wrong assumption')
        assert(any(ismember(setdiff(Tmeg.ChannelName,mylbl),"MRJ3")),'Wrong assumption')
        assert(C1<C2,'Wrong assumptions about channel order')
        b = nan(1,nch);
        for ich = 1:nch
            if fry.initial_brain_sensors(ich,1) < C1 
                b(ich) = fry.initial_brain_sensors(ich,1);
            elseif fry.initial_brain_sensors(ich,1) >= C1 && fry.initial_brain_sensors(ich,1) < C2
                b(ich) = fry.initial_brain_sensors(ich,1)+1;
            elseif fry.initial_brain_sensors(ich,1) >= C2
                b(ich) = fry.initial_brain_sensors(ich,1)+2;
            else
                error('Something is coded wrong, there shouldn''t be any logical possibilities left')
            end
        end
    else
        error('Bad number of channels')
    end

    % Double check the channels
    for ich = 1:nch
        try
            assert(fry.initial_brain_sensors(ich,2) == T.X( b(ich) ),'Channel locations don''t match')
            assert(fry.initial_brain_sensors(ich,3) == T.Y( b(ich) ),'Channel locations don''t match')
            assert(fry.initial_brain_sensors(ich,4) == T.Z( b(ich) ),'Channel locations don''t match')
        catch
            try
                b(ich) = b(ich) - 1;
                assert(fry.initial_brain_sensors(ich,2) == T.X( b(ich) ),'Channel locations don''t match')
                assert(fry.initial_brain_sensors(ich,3) == T.Y( b(ich) ),'Channel locations don''t match')
                assert(fry.initial_brain_sensors(ich,4) == T.Z( b(ich) ),'Channel locations don''t match')

            catch
                b(ich) = b(ich) - 1;
                assert(fry.initial_brain_sensors(ich,2) == T.X( b(ich) ),'Channel locations don''t match')
                assert(fry.initial_brain_sensors(ich,3) == T.Y( b(ich) ),'Channel locations don''t match')
                assert(fry.initial_brain_sensors(ich,4) == T.Z( b(ich) ),'Channel locations don''t match')
            end
        end
    end
    brain(idir,:) = b;

    % get the o-information from fetal brain channels
    gehirn = dat(brain(idir,:),:);
    SIGMA = cov(zscore(gehirn'));
    JH(idir) = -log(det(SIGMA))/2; % EQ. 3 in Varley et al. 2024 (joint entropy of the system)
    Omega(idir) = soinfo_from_covmat(SIGMA,length(gehirn));

%     for irs = 1:Nrs
%         rng(irs+seed(idir),'twister')
%         checkme = 0;
%         sanity = 0;
%         while checkme <= 0 | checkme == inf | isnan(checkme)
%             resample = randperm(size(dat,1),nch);
%             datrs = dat(resample,:).*1e15;
%             NDX(irs,:,idir) = resample; % save indices
%             SIGMA = cov(zscore(datrs'));
%             checkme = det(SIGMA);
%             sanity = sanity + 1;
%             if sanity == limit
%                 fprintf('Giving up after %i attempts\n',limit)
%                 next_idir = true;  % Set flag to break out of both loops
%                 break
%             end
%         end
% 
%         if next_idir
%             % Break out of both loops
%             break
%         end
% 
%         if sanity < limit
%             JHrs(idir,irs) = -log(det(SIGMA))/2; % EQ. 3 in Varley et al. 2024
%             assert(length(datrs)==length(gehirn),'resample has different data length than original data')
%             Omegars(idir,irs) = soinfo_from_covmat(SIGMA,length(datrs));
% 
%             assert(~isnan(Omegars(idir,irs)),'Nan output for o-info')
%             assert(~isinf(Omegars(idir,irs)),'Inf output for o-info')
% 
%             if mod(irs,10)==0, fprintf('.'), end
%         end
%     end

    if next_idir
        % Continue to the next iteration of idir
        continue;
    end
end

clear dat datrs MEG
save 'Oinfo_OP_brain_chans.mat'

%%

% fix corrupted data
Omegars = real(Omegars);
dump = any(isnan(Omegars),2);
Omegars(dump,:) = [];
Gage(dump) = [];
NDX(:,:,dump) = [];

R = nan(1,Nrs);
for irs = 1:Nrs
    R(irs) = corr(Gage',Omegars(:,irs),'rows','complete');
end

% actual correlation
[r,p] = corr(Gage',Omega')
myfigure
scatter(Gage,Omega,100,'MarkerEdgeColor','k','MarkerFaceColor','c','LineWidth',2)
mylsline
xlabel('Gestational age (weeks)')
ylabel('\Omega(X) (O-information)')
title(sprintf('fMEG OP processing, r = %1.2f, P = %1.3f',r,p),'FontSize',18)
makefigpretty
print('-dpng','./Fairy_OPproc_Oinfo.png')

myfigure
histogram(R,50)
xlim([-0.5 0.5])
%title(sprintf('Actual r = %1.2f, P = %1.2f (one-tailed empirical)',r,Prs),'fontsize',18)
xlabel('Resampled Pearson coefficient (r)')
ylabel('Number of resamples')
makefigpretty
print('-dpng','./Oinfo_correlation_resampling_OP_processing.png')

%% Now see which channels were drawn that give a correlation

idx = find(R < -0.4);
squeeze(NDX(idx,:,1))


%% Helper function
% function[omega] = oinfo(dat)
% 
% nch = size(dat,1);
% SIGMA = cov(zscore(dat'));
% JH = -log(det(SIGMA))/2; % EQ. 3 in Varley et al. 2024
% tmp = (2 - nch)*JH;
% chidx = 1:nch;
% parfor ich = 1:nch
%     idx = setdiff(chidx,ich);
%     sig = cov(zscore(dat(idx,:)'));
%     % correction from Pedro
%     tmp = tmp + log(det(SIGMA(ich,ich))/2 - log(det(sig))/2);
% end
% omega = tmp;
% 
% end




