% O-information analysis of fMEG
% Joel Frohlich
% 25.03.24
%
% See Varley et al. 2024 BioRxiv for O-information background
% https://www.biorxiv.org/content/10.1101/2024.02.23.581375v1.full.pdf

%% Compute Omega (X) for MEG with Fairy processing algorithm
clearvars
close all
rng(329014901)

try
    parpool
catch
    fprintf('Seems like parallel pool is already open\n')
end

%%% For CSER
% Load Octave packages if relevant
if exist('OCTAVE_VERSION', 'builtin') && ~exist('fcdf', 'file')
  pkg('load', 'statistics');
end

% Silently add MVGC to current path
try
    p = mfilename('fullpath');
    addpath(strrep(p, 'myStateSpaceEntropyRate', 'private/mvgc_v2.0'));
    evalc('mvgc_startup;');
catch
    error('Failled to add multivariate ganger causality to path')
end
%%%

load('megchanlocs.mat');
folder_fairy = dir('O_information/Moser_data_FAIRY');

nch = 10;
fs = 610.3516;
Gage = nan(1,length(folder_fairy));
Omega = nan(1,length(folder_fairy));
JH = nan(1,length(folder_fairy));
Omega2 = nan(1,length(folder_fairy));
JH2 = nan(1,length(folder_fairy));
Omega3 = nan(1,length(folder_fairy));
JH3 = nan(1,length(folder_fairy));

CSERbb = nan(1,length(folder_fairy));
CSERnb = nan(1,length(folder_fairy));

% CSER parameters
fctr = 1170; % trial size (1.9 sec)
morder = 10;

skipme = [];

for idir = 1:length(folder_fairy)
    if length(folder_fairy(idir).name) < 15, continue, end
    fprintf('load %ith file...\n',idir)
    pth = sprintf('./O_information/Moser_data_FAIRY/%s/',folder_fairy(idir).name);
    dat = load(sprintf('%s3_fMEG/brain.mat',pth));
    dat = dat.brain.*1e15;

    assert(size(dat,1)==nch,'Wrong number of channels')
    ga = load(sprintf('%s/0_info/ssw.mat',pth));
    Gage(idir) = str2double(ga.ssw);
    SIGMA = cov(zscore(dat'));
    JH(idir) = -log(det(SIGMA))/2; % EQ. 3 in Varley et al. 2024
    Omega(idir) = soinfo_from_covmat(SIGMA,length(dat));


    % compare with data that still has fetal heart but no maternal heart
    loadme = sprintf('./O_information/Moser_data_FAIRY/%s/1_mMCG/raw_mMCG_data.ds/raw_mMCG_data.meg4',folder_fairy(idir).name);
    MCG = pop_fileio(loadme);
    % get relevant channels
    fry = load(sprintf('%s3_fMEG/initial_brain_sensors.mat',pth));
    dat2 = MCG.data(fry.initial_brain_sensors(:,1),:);
    assert(size(dat2,1)==nch,'Wrong number of channels')
    SIGMA2 = cov(zscore(dat2'));
    JH2(idir) = -log(det(SIGMA2))/2; % EQ. 3 in Varley et al. 2024
    Omega2(idir) = soinfo_from_covmat(SIGMA2,length(dat2));


    % pure heart gives a covariance matrix with an insanely small
    % determinant, which messes up the calculations ... thus, this block of
    % code is commented out for now. 

%     % compare with pure fetal heart
%     loadme = sprintf('./O_information/Moser_data_FAIRY/%s/2_fMCG/pure_fMCG_art.ds/pure_fMCG_art.meg4',folder_fairy(idir).name);
%     MCG = pop_fileio(loadme);
%     % get relevant channels
%     fry = load(sprintf('%s3_fMEG/initial_brain_sensors.mat',pth));
%     dat3 = MCG.data(fry.initial_brain_sensors(:,1),:);
%     assert(size(dat3,1)==nch,'Wrong number of channels')
%     SIGMA3 = cov(zscore(dat3'));
%     JH3(idir) = -log(det(SIGMA3))/2; % EQ. 3 in Varley et al. 2024
%     Omega3(idir) = soinfo_from_covmat(SIGMA3,length(dat3));  
%     for ich = 1:nch
%         for jch = 1:nch
%             MI3(ich,jch,idir) = mi(dat3(ich,:)',dat3(jch,:)');
%         end
%     end


    try
        assert(mod(size(dat,2),fctr)==0)
    catch
        fprintf('Wrong data size, skipping this one for CSER ...\n')
        skipme = [skipme idir];
        continue
    end

    layers = size(dat,2)/fctr;
    d = reshape(dat,10,fctr,layers);

    % random subsampling of data trials to reduce computational time
    [H1,H2] = myStateSpaceEntropyRate(d(:,:,randperm(layers,60)), fs,'yes',[1 15],morder);
    CSERbb(idir) = H1;
    CSERnb(idir) = H2;
    fprintf('     %1.2f%% complete\n\n',idir/length(folder_fairy)*100)
end

%% plotting

% Uncomment below to trim outliers
% [~,rx] = rmoutliers([Gage' CSERnb' Omega'],'ThresholdFactor',4);
% Gage(rx) = [];
% CSERnb(rx) = [];
% Omega(rx) = [];

[r,p] = corr(Gage',Omega','rows','complete');
myfigure
scatter(Gage,Omega,100,'MarkerEdgeColor','k','MarkerFaceColor','c','LineWidth',2)
mylsline
xlabel('Gestational age (weeks)')
ylabel('\Omega(X) (O-information)')
title(sprintf('fMEG Fairy processing, r = %1.2f, P = %1.5f',r,p),'FontSize',18)
makefigpretty
print('-dpng','./Fairy_Oinfo.png')

[r,p] = corr(Gage',Omega2','rows','complete');
myfigure
scatter(Gage,Omega2,100,'MarkerEdgeColor','k','MarkerFaceColor','c','LineWidth',2)
mylsline
xlabel('Gestational age (weeks)')
ylabel('\Omega(X) (O-information)')
title(sprintf('fMEG Fairy processing with fetal heart, r = %1.2f, P = %1.5f',r,p),'FontSize',18)
makefigpretty
print('-dpng','./Fairy_Oinfo_with_heart.png')
 


%% plotting

[r,p] = corr(Gage',CSERnb','rows','complete');
myfigure
scatter(Gage,CSERnb,100,'MarkerEdgeColor','k','MarkerFaceColor','c','LineWidth',2)
mylsline
xlabel('Gestational age (weeks)')
ylabel('CSER (state space entropy rate, 1 - 15 Hz)')
title(sprintf('fMEG Fairy processing, r = %1.2f, P = %1.3f',r,p),'FontSize',18)
makefigpretty
print('-dpng','./Fairy_CSER.png')

% jack-knief to find the min correlation with leave-one-out approach

useme = 1:length(CSERnb);
R = nan(1,length(CSERnb));
P = nan(1,length(CSERnb));

for i = 1:length(CSERnb)
    idx = setdiff(useme,i);
    [r,p] = corr(Gage(idx)',CSERnb(idx)','rows','complete');
    R(i) = r; P(i) = p;
end



[a,b] = min(abs(R));

fprintf('Minimum jackknief correlation: r = %1.2f, p = %1.3f\n',R(b),P(b))


%% Draw random subsamples and compare the results from neural channels to these random samplings
Nrs = 1000; % N resamples

%%% For CSER
% Load Octave packages if relevant
if exist('OCTAVE_VERSION', 'builtin') && ~exist('fcdf', 'file')
  pkg('load', 'statistics');
end

% Silently add MVGC to current path
try
    p = mfilename('fullpath');
    addpath(strrep(p, 'myStateSpaceEntropyRate', 'private/mvgc_v2.0'));
    evalc('mvgc_startup;');
catch
    % error('I can''t find MVGC')
end
%%%

folder_fairy = dir('O_information/Moser_data_FAIRY');
Omegars = nan(length(folder_fairy),Nrs);
%Omega_Cardiacrs = nan(length(folder_fairy),Nrs);
JHrs = nan(length(folder_fairy),Nrs);
%JHrs2 = nan(length(folder_fairy),Nrs);
% MI = nan(10,10,length(folder_fairy));
% MI2 = nan(10,10,length(folder_fairy),Nrs);
% MI3 = nan(10,10,length(folder_fairy),Nrs);

CSERrs = nan(length(folder_fairy),Nrs);
seed = round(rand(1,length(folder_fairy))*1000);
folder_fairy = dir('O_information/Moser_data_FAIRY');

start_ndx = 1;
parfor idir = start_ndx:length(folder_fairy)
    if length(folder_fairy(idir).name) < 15, continue, end
    fprintf('load %ith file...\n',idir)
    % Resampling from all channels
    loadme = sprintf('./O_information/Moser_data_FAIRY/%s/3_fMEG/fMEG_data.ds/fMEG_data.meg4',folder_fairy(idir).name);
    MEG = pop_fileio(loadme);
    dat = MEG.data;

%     % Resampling from pure (back-projected?) fetal MCG (cardiac) 
%     loadme = sprintf('./O_information/Moser_data_FAIRY/%s/2_fMCG/pure_fMCG_art.ds/pure_fMCG_art.meg4',folder_fairy(idir).name);
%     MCG = pop_fileio(loadme);
%     dat2 = MCG.data;
%     pth = sprintf('./O_information/Moser_data_FAIRY/%s/',folder_fairy(idir).name);
%     ga = load(sprintf('%s/0_info/ssw.mat',pth));
%     assert(Gage(idir) == str2double(ga.ssw),'Somehow these files are indexed differently?!')


   next_idir = false;  % Flag for passing control to next iteration of idir

   for irs = 1:Nrs
        rng(irs+seed(idir),'twister')
        checkme = 0;
        sanity = 0;
        limit = 100;
        while checkme <= 0 || checkme == inf || isnan(checkme)
            resample = randperm(size(dat,1),nch);
            datrs = dat(resample,:).*1e15;
            %datrs2 = dat2(resample,:).*1e15;
            local_SIGMA = cov(zscore(datrs'));
            checkme = det(local_SIGMA);
            sanity = sanity + 1;
            if sanity == limit
                fprintf('Giving up after %i attempts\n',limit)
                next_idir = true;  % Set flag to break out of both loops
                break
            end
        end

        if next_idir
            % Break out of both loops
            break
        end

        if sanity < limit
            % Resampling from fairy-cleaned data
            SIGMA = cov(zscore(datrs'));
            JHrs(idir,irs) = -log(det(SIGMA))/2; % EQ. 3 in Varley et al. 2024
            assert(length(datrs)==length(dat),'resample has different size than original data')
            Omegars(idir,irs) = soinfo_from_covmat(SIGMA,length(datrs));

            assert(~isnan(Omegars(idir,irs)),'Nan output for o-info')
            assert(~isinf(Omegars(idir,irs)),'Inf output for o-info')


%             % Resampling from cardiac activity 
% 
%             SIGMA2 = cov(zscore(datrs2'));
%             JHrs2(idir,irs) = -log(det(SIGMA2))/2; % EQ. 3 in Varley et al. 2024
%             assert(length(datrs2)==length(dat2),'resample has different size than original data')
%             Omega_Cardiacrs(idir,irs) = soinfo_from_covmat(SIGMA2,length(datrs2));
% 
%             assert(~isnan(Omega_Cardiacrs(idir,irs)),'Nan output for o-info cardiac')
%             assert(~isinf(Omega_Cardiacrs(idir,irs)),'Inf output for o-info cardiac')
            
            try
                assert(mod(size(datrs,2),fctr)==0)
            catch
                continue
            end

            layers = size(datrs,2)/fctr; 
            d = reshape(datrs,10,fctr,layers);

            % random subsampling of data trials to reduce computational time
            [H1,H2] = myStateSpaceEntropyRate(d(:,:,randperm(layers,60)), fs,'yes',[1 15],morder);    
            CSERrs_bb(idir,irs) = H1;
            CSERrs_nb(idir,irs) = H2;
            if mod(irs,10)==0, fprintf('.'), end
        end
    end

    if next_idir
        % Continue to the next iteration of idir
        continue;
    end
end

clear dat datrs MEG
save 'Oinfo_Fairy_brain_chans.mat'

%% Empirical P-value for correlation

good = any(~isnan(Omegars),2);

[r,p] = corr(Gage(good)',Omega(good)');

%myfigure

R = nan(1,Nrs);
for irs = 1:Nrs
    R(irs) = corr(Gage(good)',Omegars(good,irs),'rows','complete');
%     if irs <= 25
%         subplot(5,5,irs)
%         scatter(Gage(good),Omegars(good,irs))
%         mylsline
%     end
end

Prs1 = sum(r >= R)/Nrs*2;
Prs2 = sum(r <= R)/Nrs*2;
Prs = min([Prs1 Prs2]);

if Prs == 0
    Prs = 1/(Nrs+1);
end

fprintf('Nonparametric P = %1.6f\n',Prs)

myfigure2
histogram(R,33)
title(sprintf('FAIRY processing: r = %1.2f, P = %1.3f (two-tailed empirical, N = %i resamples)',r,Prs,Nrs),'fontsize',18)
xlabel('Resampled Pearson coefficient (r)')
ylabel('Number of resamples')
makefigpretty
print('-dpng','./Oinfo_correlation_resampling.png')

%% 

Rcser = nan(1,Nrs);

for irs = 1:Nrs
    Rcser(irs) = corr(CSERrs_nb(:,irs),Gage2','rows','complete');
end


%% Plotting

P = sum(Omega' <= Omegars,2)/Nrs;
P(P==0) = 1/(Nrs+1);

myfigure
scatter(Gage,-log10(P))
mylsline
[r,p] = corr(Gage',-log10(P),'rows','complete')




% OSize = nan(1,length(folder_fairy));
% OMax  = nan(1,length(folder_fairy));
%
% for idir = 1:length(folder_fairy)
%     allO = [];
%     allS = [];
%     if length(folder_fairy(idir).name) < 15, continue, end
%     fprintf('loading %ith file...\n',idir)
%     pth = sprintf('./O_information/Moser_data_FAIRY/%s/',folder_fairy(idir).name);
%     dat = load(sprintf('%s3_fMEG/brain.mat',pth));
%     dat = dat.brain;
%     assert(size(dat,1)==nch,'Wrong number of channels')
%     % annealing analysis
%     parfor k = 3:nch
%         tmp = nchoosek(1:nch, k);
%         NDX = unique(sort(tmp),'rows');
%         for indx = 1:size(NDX,1)
%             fprintf('.')
%             allO = [allO soinfo_from_covmat(dat(NDX(indx,:),:))];
%             allS = [allS length(NDX(indx,:))];
%         end
%     end
%     [a,b] = min(allO);
%     OMax(idir) = a;
%     Osize(idir) = allS(b);
%
%     fprintf('\n     %1.2f%% complete\n\n',idir/length(folder_fairy)*100)
% end
%
% save Oinfo_Fairy_brain_chans
%
% return






