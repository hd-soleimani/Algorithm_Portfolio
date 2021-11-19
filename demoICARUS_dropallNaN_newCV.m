% ==============================================================================
% demoICARUS
% ==============================================================================
%
% By: Mario Andrés Muñoz Acosta
%     The University of Melbourne
%     Australia
%

% This is to set up the file reading
rootdir = './';
filelist = struct2cell(dir([rootdir 'metadata_*.csv']));
filelist = filelist(1,:)';
nfiles = length(filelist);
nfolds  = 10;
nrandreps = 30;
a = 1;

for ii=1:nfiles
    % Reading the crossvalidation file
    CV = readtable([rootdir 'CV_' filelist{ii}],'Delimiter',',');
    CV = CV.CV;
    
    % Reading the file
    Xbar = readtable(filelist{ii});
    varlabels = Xbar.Properties.VariableNames;

    instlabels = Xbar{:,strcmpi(varlabels,'instances')};
    if isnumeric(instlabels)
        instlabels = num2cell(instlabels);
        instlabels = cellfun(@(x) num2str(x),instlabels,'UniformOutput',false);
    end

    issource = strcmpi(varlabels,'source');
    if any(issource)
        S = categorical(Xbar{:,issource});
    end

    emptystr = '';
    featind = 'feature_';
    algoind = 'algo_';

    isfeat = strncmpi(varlabels,featind,length(featind));
    featlabels = strrep(varlabels(isfeat),featind,emptystr);
    X = Xbar{:,isfeat};

    isalgo = strncmpi(varlabels,algoind,length(algoind));
    algolabels = strrep(varlabels(isalgo),algoind,emptystr);
    Y = Xbar{:,isalgo};
    

    % Only use those problems for which we do have algorithm results. Ignore
    % otherwise.
    nosolution = all(isnan(Y),2);
    X(nosolution,:) = [];
    Y(nosolution,:) = [];
    instlabels(nosolution) = [];
    if any(issource)
        S(nosolution) = [];
    end

    % -------------------------------------------------------------------------
    % Preprocessing the data 
    epsilon = [0.00 0.05 0.10 0.20 0.50 1.00 2.00]; % Determines the acceptable underperformance of an algorithm from the best one
    MaxPerf = false; % True for maximizing performance, false otherwise
    AbsPerf = false; % True for absolute performance, false for relative
    dataRND = [];
    dataIC = [];
    dataSFS = [];
    
    for jj=1:length(epsilon)
        % Use the definition of 'good' performance to determine if the algorithm
        % was succesful or not.
        if MaxPerf
            Yaux = Y;
            Yaux(isnan(Yaux)) = -1e9;
            [rankPerf,rankAlgo] = sort(Yaux,2,'descend');
            Ybest = rankPerf(:,1);
            P = rankAlgo(:,1);
            if AbsPerf
                Ybin = Yaux>=epsilon(jj);
            else
                Yaux(Yaux==0) = eps;
                Ybest(Ybest==0) = eps;
                Ybin = (1-bsxfun(@rdivide,Yaux,Ybest))<=epsilon(jj);
            end
        else
            Yaux = Y;
            Yaux(isnan(Yaux)) = 1e9;
            [rankPerf,rankAlgo] = sort(Yaux,2,'ascend');
            Ybest = rankPerf(:,1);
            P = rankAlgo(:,1);
            if AbsPerf
                Ybin = Yaux<=epsilon(jj);
            else
                Yaux(Yaux<=eps) = eps;
                Ybest(Ybest<=eps) = eps;
                Ybin = (bsxfun(@rdivide,Yaux,Ybest)-1)<=epsilon(jj);
            end
        end

        % Create the structures that contain the preferences
        [ninsts,nalgos] = size(Y);
        inputRanks = zeros(nalgos,4,ninsts);
        algoranks = zeros(ninsts,nalgos);
        votes = algoranks;

        % INPUT RANKS IS ORGANIZED AS FOLLOW:
        % (1 COL) RANKING OF THE ALGORITHM ON THE CURRENT SET
        % (2 COL) INDEX OF THE ALGORITHM ON THE COMPLETE SET
        % (3 COL) BINARY RESULT, I.E., IS THE ALGORITHM "SUCCESFUL" (0) OR NOT (1)
        % (4 COL ONWARDS) PERFORMANCE CRITERIA. THE RANKINGS ARE ORGANIZED WITH
        % INCREASINGLY IMPORTANT CRITERIA AS TO BREAK TIES. THEREFORE, THERE IS
        % SCOPE TO INCREASE THE CRITERIA BY ADDING MORE COLUMNS.
        inputRanks(:,3,:) = double(~Ybin'); % In MATILDA, (1) are good, we need to swap them around
        inputRanks(:,4,:) = (Y./Ybest)';

        % Estimating the preferences
        for kk=1:ninsts
            if MaxPerf
                [aux,algoranks(kk,:)] = sortrows(inputRanks(:,3:end,kk),'descend','MissingPlacement','last');
            else
                [aux,algoranks(kk,:)] = sortrows(inputRanks(:,3:end,kk),'ascend','MissingPlacement','last');
            end
            votes(kk,:) = algoranks(kk,:);
            votes(kk,aux(:,1)==1) = NaN;
            [inputRanks(:,2,kk),inputRanks(:,1,kk)] = sort(algoranks(kk,:));
        end

        % Setting up ICARUS to run
        portfolioSizeICARUS = zeros(nalgos,nfolds);
        portfolioRegretICARUS = portfolioSizeICARUS;
        portfolioICARUS = cell(nalgos,nfolds);
        outputTrainRanksICARUS = portfolioICARUS;
        outputTestRanksICARUS = portfolioICARUS;
        
        portfolioRegretSFS = portfolioSizeICARUS;
        portfolioSFS = portfolioICARUS;
        outputTrainRanksSFS = portfolioICARUS;
        outputTestRanksSFS = portfolioICARUS;
        
        portfolioRegretRND = portfolioSizeICARUS;
        portfolioRND = portfolioICARUS;
        outputRanksRND = portfolioICARUS;
        
        portfolioVectorICARUS = zeros(nalgos,nalgos);
        portfolioVectorSFS = zeros(nalgos,nalgos);
        portfolioVectorRND = zeros(nalgos,nalgos);
        
        for kk=1:nalgos
            
            for ll=1:nrandreps
                % RANDOM SELECTION
                [portfolioRND{kk,ll},outputRanksRND{kk,ll}] = randps(inputRanks,kk);
                [~,portfolioRegretRND(kk,ll)] = getRegret(inputRanks,portfolioRND{kk,ll},true(ninsts,1),Yaux,Ybest);
                
                aux = false(1,nalgos);
                aux(portfolioRND{kk,ll}) = true;
                portfolioVectorRND(kk,:) = portfolioVectorRND(kk,:) + aux;
            end
            
            for ll=1:nfolds
                [portfolioICARUS{kk,ll},outputTrainRanksICARUS{kk,ll}] = icarus(votes(CV~=ll,:),...
                                                                           inputRanks(:,:,CV~=ll),kk);
                portfolioSizeICARUS(kk,ll) = length(portfolioICARUS{kk,ll});
                [outputTestRanksICARUS{kk,ll},...
                    portfolioRegretICARUS(kk,ll)] = getRegret(inputRanks,...
                                                              portfolioICARUS{kk,ll},...
                                                              CV==ll,Yaux,Ybest);
                                                          
                aux = false(1,nalgos);
                aux(portfolioICARUS{kk,ll}) = true;
                portfolioVectorICARUS(kk,:) = portfolioVectorICARUS(kk,:) + aux;
                
                % SEQUENTIAL FORWARD SELECTION
                [portfolioSFS{kk,ll},outputTrainRanksSFS{kk,ll}] = psfs(inputRanks(:,:,CV~=ll),kk);
                [outputTestRanksSFS{kk,ll},...
                    portfolioRegretSFS(kk,ll)] = getRegret(inputRanks,...
                                                              portfolioSFS{kk,ll},...
                                                              CV==ll,Yaux,Ybest);
                aux = false(1,nalgos);
                aux(portfolioSFS{kk,ll}) = true;
                portfolioVectorSFS(kk,:) = portfolioVectorSFS(kk,:) + aux;
            end
            
        end
        
        epscol = epsilon(jj).*ones(nalgos,1);
        
        dataRND = [dataRND; epscol (1:nalgos)' mean(portfolioRegretRND,2) std(portfolioRegretRND,[],2) ...
                   portfolioVectorRND./nrandreps];
        
        
        dataSFS = [dataSFS; epscol (1:nalgos)' mean(portfolioRegretSFS,2) std(portfolioRegretSFS,[],2) ...
                   portfolioVectorSFS./nfolds];
        
        aux = [epscol (1:nalgos)' mode(portfolioSizeICARUS,2) mean(portfolioRegretICARUS,2) ...
               std(portfolioRegretICARUS,[],2) portfolioVectorICARUS./nfolds];
        uniqueCardinal = unique(aux(:,3))';
        idx = 0.*uniqueCardinal;
        for kk=1:length(uniqueCardinal)
            idx(kk) = find(ismember(aux(:,3),uniqueCardinal(kk)),1);
        end
        dataIC = [dataIC; aux(idx,:)];

%         save([rootdir filelist{ii}(1:end-4) '_EPS' num2str(jj) '.mat'],...
%             'portfolioSizeICARUS','portfolioRegretICARUS','portfolioICARUS',...
%             'outputTrainRanksICARUS','outputTestRanksICARUS','portfolioRegretSFS',...
%             'portfolioSFS','outputTrainRanksSFS','outputTestRanksSFS',...
%             'portfolioRegretRND','portfolioRND','outputRanksRND');
    end
    
    lbls = cell(1,4+length(algolabels));
    lbls(1:4) = {'Epsilon','Cardinality','CV_prediction_mean','CV_prediction_std'};
    lbls(5:end) = algolabels;
    writetable(array2table(dataRND,'VariableNames',lbls),...
               'RND.xls',...
               'WriteVariableNames',true,...
               'sheet',filelist{ii}(10:end-4));
    
    lbls = cell(1,4+length(algolabels));
    lbls(1:4) = {'Epsilon','Cardinality','CV_prediction_mean','CV_prediction_std'};
    lbls(5:end) = algolabels;
    writetable(array2table(dataSFS,'VariableNames',lbls),...
               'PSFS.xls',...
               'WriteVariableNames',true,...
               'sheet',filelist{ii}(10:end-4));
    
    lbls = cell(1,5+length(algolabels));
    lbls(1:5) = {'Epsilon','N_Preferences','Cardinality','CV_prediction_mean','CV_prediction_std'};
    lbls(6:end) = algolabels;
    writetable(array2table(dataIC,'VariableNames',lbls),...
              'ICARUS.xls',...
              'WriteVariableNames',true,...
              'sheet',filelist{ii}(10:end-4));
end


%         % Post-process the results from ICARUS, find the smallest oracle given the
%         % preferences.
%         smallestPortfolioSize = find(portfolioSizeICARUS==min(portfolioSizeICARUS),1,'first');
%         smallestPorftolioRanks = outputRanksICARUS{smallestPortfolioSize};
%         smallestPorfolio = portfolioICARUS{smallestPortfolioSize};
% 
%         disp(algolabels(smallestPorfolio));
% 
%         for kk=1:ninsts
%             smallestPorftolioRanks(:,:,kk) = sortrows(smallestPorftolioRanks(:,:,kk),1);
%         end
% 
%         % Estimate the expected performance of this portfolio
%         selectedAlgorithm = squeeze(smallestPorftolioRanks(1,2,:));
%         Yportfolio = Yaux(sub2ind(size(Yaux),(1:ninsts)',selectedAlgorithm));
%         % Calculate the max 'regret' (As I understood)
%         maxRegret = max(Yportfolio./Ybest);
% ==============================================================================
% Subfunctions
% ==============================================================================
function [portfolio,outputRanks] = icarus(votes,inputRanks,nValidPrefs)
    % Automated Oracle Construction using ICARUS (Identification of
    % Complementary Algorithms By Uncovered Sets), an algorithm based on
    % election systems to identify potentially robust algorithm portafolios.
    %
    % Note: Uses code from 'election.m' by Ben Petschel.
    %

    [~,ncrit,ninsts] = size(inputRanks);
    if nargin~=3
        error('AUTOAlgoSel:ICABUS:notEnoughInputs','Not enough input arguments.');
    end
    processedVotes = votes(:,1:nValidPrefs);
    portfolio = landau(processedVotes); %Using Landau voting system
    while ~isempty(processedVotes)
        outputRanks = zeros(length(portfolio),ncrit,ninsts);
        store = outputRanks;
        for ii=1:ninsts
            outputRanks(:,:,ii) = inputRanks(portfolio,:,ii);
            [~,ordering] = sortrows(outputRanks(:,3:end,ii),'MissingPlacement','last');
            [~,outputRanks(:,1,ii)] = sort(ordering);
            store(:,:,ii) = sortrows(outputRanks(:,:,ii),'MissingPlacement','last');
        end
        processedVotes = votes((squeeze(store(1,3,:)==1)),1:nValidPrefs);
        portfolio = unique([portfolio landau(processedVotes)]);
    end
end
% -------------------------------------------------------------------------
function [portfolio,outputRanks] = psfs(inputRanks,maxPortfolioSize)
    outputRanks = inputRanks;
    [nalgos,~,ninst] = size(inputRanks);
    algorithmSet = 1:nalgos;
    [~,averageRank] = sortrows(nanmean(inputRanks(:,3:end,:),3));
    current = averageRank(1);
    portfolio = zeros(1,maxPortfolioSize);
    portfolio(1) = algorithmSet(current);
    for ii=2:maxPortfolioSize
        algorithmSet = setdiff(algorithmSet,portfolio);
        inputRanks(current,:,:) = [];
        [~,averageRank] = sortrows(nanmean(inputRanks(:,3:end,:),3));
        current = averageRank(1);
        portfolio(ii) = algorithmSet(current);
    end
    outputRanks(setdiff(1:nalgos,portfolio),:,:) = [];
    for ii=1:ninst
        [~,ordering]           = sortrows(outputRanks(:,3:end,ii),'MissingPlacement','last');
        [~,outputRanks(:,1,ii)] = sort(ordering);
    end
end
% -------------------------------------------------------------------------
function [portfolio,outputRanks] = randps(inputRanks,maxPortfolioSize)
    [nalgos,~,ninst] = size(inputRanks);
    portfolio = randperm(nalgos,maxPortfolioSize);
    outputRanks = inputRanks(portfolio,:,:);
    for ii=1:ninst
        [~,ordering]           = sortrows(outputRanks(:,3:end,ii),'MissingPlacement','last');
        [~,outputRanks(:,1,ii)] = sort(ordering);
    end
end
% -------------------------------------------------------------------------
function [outputRanks,maxRegret] = getRegret(inputRanks,portfolio,CVidx,Y,Ybest)

    outputRanks = inputRanks(portfolio,:,CVidx);
    for ii=1:size(outputRanks,3)
        [~,ordering] = sortrows(outputRanks(:,3:end,ii),'MissingPlacement','last');
        [~,outputRanks(:,1,ii)] = sort(ordering);
    end
    aux = squeeze(outputRanks(:,1,:)==1);
    if size(aux,1)~=length(portfolio)
        aux = aux';
    end
    [aux,~] =  find(aux);
    Yportfolio = Y(sub2ind(size(Y),find(CVidx)', portfolio(aux)))';
    maxRegret = max(Yportfolio./Ybest(CVidx));

end
% -------------------------------------------------------------------------
% From 'election.m' by Ben Petschel.
% -------------------------------------------------------------------------
function winner = landau(pref)
    % 'Landau' - (also known as uncovered set or Fishburn set): a set of candidates
    % that beats or ties all other candidates, or for any they don't, beats/ties
    % with a third candidate that beats/ties with the other.  In other words, I is a
    % winner if for all J, C(I,J)>=C(J,I) or C(I,K)>=C(K,I) and C(K,J)>=C(J,K) for 
    % some K.
    B = beattiematrix(pref);
    if isempty(B)
      winner = []; % no valid votes
      return
    end
    % matrix B + B^2 is positive if B(i,j) or B(i,k)&B(k,j) is true for some k
    winner = find(all((B+B^2)>0, 2)');
end
% -------------------------------------------------------------------------
function B=beattiematrix(pref)
    % Calculates the matrix B where B(i,j) is true if C(i,j)>=C(j,i) where C is the
    % matrix with the number of voters that prefer i to j (see Cmatrix for details),
    % i.e. B(i,j) is true if candidate i beats or ties with candidate j.
    C = Cmatrix(pref);
    B = C>=C';
end
% -------------------------------------------------------------------------
function C=Cmatrix(pref)
    % Determines the matrix where C(i,j) is the number of voters that prefer i to j.
    % Any candidates not listed by a voter are assumed to be equal least preferred.
    %
    % Any vote with repeats or non-candidates out of order (non-positive-integers 
    % listed before the last listed positive integer) are ignored.
    N = max(pref(isposint(pref))); % number of candidates
    C = zeros(N,N);
    for k=1:size(pref,1)
      thisvote = pref(k,:);
      iscand = isposint(thisvote);
      lastcand = find(iscand,1,'last');
      firstnon = find(~iscand,1,'first');
      if ~isempty(lastcand) && (isempty(firstnon) || lastcand>firstnon)
        % valid vote must have candidates and no non-candidates out of order
        thisvote = thisvote(1:lastcand);
        nopref = setdiff(1:N,thisvote); % candidates not expressed in preferences
        for i=1:numel(thisvote)
          % add vote for every candidate less prefered than the i'th
          C(thisvote(i),thisvote(i+1:end))=C(thisvote(i),thisvote(i+1:end))+1;
          % add vote for every candidate not explicitly listed
          C(thisvote(i),nopref) = C(thisvote(i),nopref)+1;
        end % for i=1:numel(thisvote)
      end % if ~isempty(...)
    end % for k=1:size(pref,1)
end
% -------------------------------------------------------------------------
function y=isposint(x)
% returns logical array with y(i)=true if x(i) is a positive integer
    y = isfinite(x) & (x>0) & (x==floor(x));
end
% -------------------------------------------------------------------------