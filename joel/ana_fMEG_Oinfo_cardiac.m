% O-information analysis of fMEG
% Joel Frohlich
% 25.03.24
%
% See Varley et al. 2024 BioRxiv for O-information background
% https://www.biorxiv.org/content/10.1101/2024.02.23.581375v1.full.pdf

%% Compute Omega (X) for MEG with Fairy processing algorithm
clearvars
close all

folder_fairy = dir('O_information/SampleData_FAIRY');
nch = 10;
fs = 610.3516;
Gage = nan(1,length(folder_fairy));
Omega = nan(1,length(folder_fairy));
JH = nan(1,length(folder_fairy));
CSERbb = nan(1,length(folder_fairy));
CSERnb = nan(1,length(folder_fairy));


for idir = 1:length(folder_fairy)
    if length(folder_fairy(idir).name) < 15, continue, end
    fprintf('load %ith file...\n',idir)
    pth = sprintf('./O_information/SampleData_FAIRY/%s/',folder_fairy(idir).name);
    dat = load(sprintf('%s2_fMCG/initial_cluster.mat',pth));
    dat = dat.initial_cluster.*1e15;
    assert(size(dat,1)==nch,'Wrong number of channels')
    ga = load(sprintf('%s/0_info/ssw.mat',pth));
    Gage(idir) = str2double(ga.ssw);
    SIGMA = cov(zscore(dat'));
    JH(idir) = -log(det(SIGMA))/2; % EQ. 3 in Varley et al. 2024
    Omega(idir) = soinfo_from_covmat(SIGMA,length(dat));
    [H1,H2] = myStateSpaceEntropyRate(dat, fs,'yes',[1 15],20);
    CSERbb(idir) = H1;
    CSERnb(idir) = H2;
    fprintf('     %1.2f%% complete\n\n',idir/length(folder_fairy)*100)
end

%% plotting

[~,rx] = rmoutliers([Gage' CSERnb' Omega'],'ThresholdFactor',4);
Gage(rx) = [];
CSERnb(rx) = [];
Omega(rx) = [];

[r,p] = corr(Gage',Omega','rows','complete');
myfigure
scatter(Gage,Omega,100,'MarkerEdgeColor','k','MarkerFaceColor','c','LineWidth',2)
mylsline
xlabel('Gestational age (weeks)')
ylabel('\Omega(X) (O-information)')
title(sprintf('fMEG Fairy processing, r = %1.2f, P = %1.3f',r,p),'FontSize',18)
makefigpretty
print('-dpng','./Fairy_Oinfo.png')

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

try
    parpool
catch
    fprintf('Seems like parallel pool is already open\n')
end

Nrs = 1000; % N resamples

folder_fairy = dir('O_information/SampleData_FAIRY');
Omegars = nan(length(folder_fairy),Nrs);
JHrs = nan(length(folder_fairy),Nrs);
CSERrs = nan(length(folder_fairy),Nrs);

folder_fairy = dir('O_information/SampleData_FAIRY');

%%
start_ndx = 1;
for idir = start_ndx:length(folder_fairy)
    if length(folder_fairy(idir).name) < 15, continue, end
    fprintf('load %ith file...\n',idir)
    loadme = sprintf('./O_information/SampleData_FAIRY/%s/2_fMCG/pure_fMCG_art.ds/pure_fMCG_art.meg4',folder_fairy(idir).name);
    MEG = pop_fileio(loadme);
    dat = MEG.data;
    pth = sprintf('./O_information/SampleData_FAIRY/%s/',folder_fairy(idir).name);
    ga = load(sprintf('%s/0_info/ssw.mat',pth));
    assert(Gage(idir) == str2double(ga.ssw),'Somehow these files are indexed differently?!')
    next_idir = false;  % Flag for passing control to next iteration of idir

    for iperm = 1:Nrs
        rng(iperm,'twister')
        checkme = 0;
        sanity = 0;
        limit = 1000;
        while checkme <= 0 || checkme == inf || isnan(checkme)
            resample = randperm(size(dat,1),nch);
            datrs = dat(resample,:).*1e15;
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
            SIGMA = cov(zscore(datrs'));
            JHrs(idir,iperm) = -log(det(SIGMA))/2; % EQ. 3 in Varley et al. 2024
            assert(length(datrs)==length(dat),'resample has different size than original data')
            Omegars(idir,iperm) = soinfo_from_covmat(SIGMA,length(datrs));

            assert(~isnan(Omegars(idir,iperm)),'Nan output for o-info')
            assert(~isinf(Omegars(idir,iperm)),'Inf output for o-info')

            % uncomment below to also do CSER again with resampling
            %[H1,H2] = myStateSpaceEntropyRate(dat, fs,'yes',[1 15],20);
            %CSERrs_bb(idir,iperm) = H1;
            %CSERrs_nb(idir,iperm) = H2;
            if mod(iperm,10)==0, fprintf('.'), end
        end
    end

    if next_idir
        % Continue to the next iteration of idir
        continue;
    end
end

clear dat datrs MEG
save 'Oinfo_Fairy_cardiac_chans.mat'

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

fprintf('Nonparametric P = %1.3f\n',Prs)

myfigure2
histogram(R,33)
title(sprintf('FAIRY processing: r = %1.2f, P = %1.2f (two-tailed empirical, N = %i resamples)',r,Prs,Nrs),'fontsize',18)
xlabel('Resampled Pearson coefficient (r)')
ylabel('Number of resamples')
makefigpretty
print('-dpng','./Oinfo_correlation_resampling_cardiac.png')

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
%     pth = sprintf('./O_information/SampleData_FAIRY/%s/',folder_fairy(idir).name);
%     dat = load(sprintf('%s3_fMEG/initial_cluster.mat',pth));
%     dat = dat.initial_cluster;
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


% %% Helper function
% function[omega] = soinfo_from_covmat(dat)
% 
% nch = size(dat,1);
% SIGMA = cov(dat');
% JH = -log(det(SIGMA))/2; % EQ. 3 in Varley et al. 2024
% tmp = (2 - nch)*JH;
% chidx = 1:nch;
% for ich = 1:nch
%     idx = setdiff(chidx,ich);
%     sig = cov(dat(idx,:)');
%     tmp = tmp + -log(det(sig))/2;
% end
% omega = tmp;
% 
% end




