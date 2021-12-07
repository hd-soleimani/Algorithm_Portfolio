% ==============================================================================
% experimentPortfolio.m
% ==============================================================================
%
% By: Mario Andrés Muñoz Acosta
%     The University of Melbourne
%     Australia
%

% This is to set up the file reading
tid = getenv('SLURM_ARRAY_TASK_ID');

if ~isempty(tid)
    disp(['Trial number: ' tid]);
    tid = str2double(tid);
    rootdir = '~/punim0320/portfolios/';
else
    rootdir = './';
end

filelist = struct2cell(dir([rootdir 'aslib/metadata_*.csv']));
filelist = filelist(1,:)';
nfiles = length(filelist);
nfolds  = 10;
nrandreps = 30;

for ii=1:nfiles
    if ~isempty(tid) && ii~=tid
        % If running on the cluster, then it skips all other jobs until it
        % finds the right one.
        continue
    end
    
    % Reading the crossvalidation file
    CV = readtable([rootdir 'cv/CV_' filelist{ii}],'Delimiter',',');
    CV = CV.CV;
    
    % Reading the file
    Xbar = readtable([rootdir '/aslib/' filelist{ii}]);
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
    dataBP = [];
    dataRND = [];
    dataSF = [];
    dataSB = [];
    dataIC = [];
    
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

        % Setting up the methods to run
        cardinalityBP = zeros(1,nfolds);
        regretBP = cardinalityBP;
        entropyBP = cardinalityBP;
        portfolioBP = cell(1,nfolds);
        trainRanksBP = portfolioBP;
        testRanksBP = portfolioBP;
        vectorBP = zeros(1,nalgos);
        
        cardinalityICARUS = zeros(nalgos,nfolds);
        regretICARUS = cardinalityICARUS;
        entropyICARUS = NaN.*cardinalityICARUS;
        portfolioICARUS = cell(nalgos,nfolds);
        trainRanksICARUS = portfolioICARUS;
        testRanksICARUS = portfolioICARUS;
        vectorICARUS = zeros(nalgos,nalgos);
        
        regretSF = cardinalityICARUS;
        entropySF = NaN.*cardinalityICARUS;
        portfolioSF = portfolioICARUS;
        trainRanksSF = portfolioICARUS;
        testRanksSF = portfolioICARUS;
        vectorSF = vectorICARUS;
        
        regretSB = cardinalityICARUS;
        entropySB = NaN.*cardinalityICARUS;
        portfolioSB = portfolioICARUS;
        trainRanksSB = portfolioICARUS;
        testRanksSB = portfolioICARUS;
        vectorSB = vectorICARUS;
        
        regretRND = zeros(nfolds,nrandreps);
        entropyRND = NaN.*regretRND;
        portfolioRND = portfolioICARUS;
        ranksRND = portfolioICARUS;
        vectorRND = vectorICARUS;
        
        for ll=1:nfolds
            % BEST PERFORMANCE
            [portfolioBP{ll},trainRanksBP{ll}] = BP(votes(CV~=ll,1),inputRanks(:,:,CV~=ll));
            cardinalityBP(ll) = length(portfolioBP{ll});
            [testRanksBP{ll},regretBP(ll)] = getRegret(inputRanks,portfolioBP{ll},...
                                                       CV==ll,Yaux,Ybest);
            vectorBP(portfolioBP{ll}) = vectorBP(portfolioBP{ll}) + 1;
            entropyBP(ll) = -getEntropy(Ybin(CV==ll,portfolioBP{ll}))./log(1./cardinalityBP(ll)); % Hportfolio./Hbest;
        end
        
        for kk=1:nalgos
            Hbest = -log(1./kk);
            
            for ll=1:nrandreps
                % RANDOM SELECTION
                [portfolioRND{kk,ll},ranksRND{kk,ll}] = randps(inputRanks,kk);
                [~,regretRND(kk,ll)] = getRegret(inputRanks,portfolioRND{kk,ll},...
                                                 true(ninsts,1),Yaux,Ybest);
                vectorRND(kk,portfolioRND{kk,ll}) = vectorRND(kk,portfolioRND{kk,ll}) + 1;
                if Hbest~=0
                    entropyRND(kk,ll) = getEntropy(Ybin(:,portfolioRND{kk,ll}))./Hbest;
                end
            end
            
            for ll=1:nfolds
                % SEQUENTIAL FORWARD SELECTION
                [portfolioSF{kk,ll},trainRanksSF{kk,ll}] = psf(inputRanks(:,:,CV~=ll),kk);
                [testRanksSF{kk,ll},regretSF(kk,ll)] = getRegret(inputRanks,...
                                                                   portfolioSF{kk,ll},...
                                                                   CV==ll,Yaux,Ybest);
                vectorSF(kk,portfolioSF{kk,ll}) = vectorSF(kk,portfolioSF{kk,ll}) + 1;
                if Hbest~=0
                    entropySF(kk,ll) = getEntropy(Ybin(CV==ll,portfolioSF{kk,ll}))./Hbest;
                end
                
                % SEQUENTIAL BACKWARD SELECTION
                [portfolioSB{kk,ll},trainRanksSB{kk,ll}] = psb(inputRanks(:,:,CV~=ll),kk);
                [testRanksSB{kk,ll},regretSB(kk,ll)] = getRegret(inputRanks,...
                                                                 portfolioSB{kk,ll},...
                                                                 CV==ll,Yaux,Ybest);
                vectorSB(kk,portfolioSB{kk,ll}) = vectorSB(kk,portfolioSB{kk,ll}) + 1;
                if Hbest~=0
                    entropySB(kk,ll) = getEntropy(Ybin(CV==ll,portfolioSB{kk,ll}))./Hbest;
                end
                
                % ICARUS
                [portfolioICARUS{kk,ll},trainRanksICARUS{kk,ll}] = icarus(votes(CV~=ll,:),...
                                                                          inputRanks(:,:,CV~=ll),kk);
                cardinalityICARUS(kk,ll) = length(portfolioICARUS{kk,ll});
                [testRanksICARUS{kk,ll},regretICARUS(kk,ll)] = getRegret(inputRanks,...
                                                                         portfolioICARUS{kk,ll},...
                                                                         CV==ll,Yaux,Ybest);
                vectorICARUS(kk,portfolioICARUS{kk,ll}) = vectorICARUS(kk,portfolioICARUS{kk,ll}) + 1;
                entropyICARUS(kk,ll) = -getEntropy(Ybin(CV==ll,portfolioICARUS{kk,ll}))./log(1./cardinalityICARUS(kk,ll)); % Hportfolio./Hbest;
            end
            
        end
        
        epscol = epsilon(jj).*ones(nalgos,1);
        
        dataRND = [dataRND; epscol (1:nalgos)' ...
                   mean(regretRND,2) std(regretRND,[],2) ...
                   mean(entropyRND,2) std(entropyRND,[],2) ...
                   vectorRND./nrandreps]; %#ok<*AGROW>
        
        dataSF = [dataSF; epscol (1:nalgos)' ...
                   mean(regretSF,2) std(regretSF,[],2) ...
                   mean(entropySF,2) std(entropySF,[],2) ...
                   vectorSF./nfolds];
        
        dataSB = [dataSB; epscol (1:nalgos)' ...
                   mean(regretSB,2) std(regretSB,[],2) ...
                   mean(entropySB,2) std(entropySB,[],2) ...
                   vectorSB./nfolds];
        
        aux = [epscol (1:nalgos)' mode(cardinalityICARUS,2) ...
               mean(regretICARUS,2) std(regretICARUS,[],2) ...
               mean(entropyICARUS,2) std(entropyICARUS,[],2) ...
               vectorICARUS./nfolds];
        uniqueCardinal = unique(aux(:,3))';
        idx = 0.*uniqueCardinal;
        for kk=1:length(uniqueCardinal)
            idx(kk) = find(ismember(aux(:,3),uniqueCardinal(kk)),1);
        end
        dataIC = [dataIC; aux(idx,:)];
        
        dataBP = [dataBP; epsilon(jj) mode(cardinalityBP,2) ...
                  mean(regretBP,2) std(regretBP,[],2) ...
                  mean(entropyBP,2) std(entropyBP,[],2) ...
                  vectorBP./nfolds];
    end
    
    lbls = cell(1,6+length(algolabels));
    lbls(1:6) = {'Epsilon','Cardinality',...
                 'CV_mean_regret','CV_std_regret',...
                 'CV_mean_entropy','CV_std_entropy'};
    lbls(7:end) = algolabels;
    writetable(array2table(dataBP,'VariableNames',lbls),...
               [rootdir 'output/BP.xls'],...
               'WriteVariableNames',true,...
               'sheet',filelist{ii}(10:end-4));
    
    writetable(array2table(dataRND,'VariableNames',lbls),...
               [rootdir 'output/RND.xls'],...
               'WriteVariableNames',true,...
               'sheet',filelist{ii}(10:end-4));
    
    writetable(array2table(dataSF,'VariableNames',lbls),...
               [rootdir 'output/SF.xls'],...
               'WriteVariableNames',true,...
               'sheet',filelist{ii}(10:end-4));
    
    writetable(array2table(dataSB,'VariableNames',lbls),...
               [rootdir 'output/SB.xls'],...
               'WriteVariableNames',true,...
               'sheet',filelist{ii}(10:end-4));
    
    lbls = cell(1,7+length(algolabels));
    lbls(1:7) = {'Epsilon','N_Preferences','Cardinality',...
                 'CV_mean_regret','CV_std_regret',...
                 'CV_mean_entropy','CV_std_entropy'};
    lbls(8:end) = algolabels;
    writetable(array2table(dataIC,'VariableNames',lbls),...
              [rootdir 'output/ICARUS.xls'],...
              'WriteVariableNames',true,...
              'sheet',filelist{ii}(10:end-4));
          
    if ~isempty(tid)
        disp(['Completed ' num2str(tid) '. Exiting...']);
        exit;   % No need to run more than once if in the cluster
    end
end
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
    portfolio = sort(portfolio);
end
% -------------------------------------------------------------------------
function [portfolio,outputRanks] = BP(votes,inputRanks)
    ninst = size(inputRanks,3);
    portfolio = unique(votes)';
    outputRanks = inputRanks(portfolio,:,:);
    for ii=1:ninst
        [~,ordering] = sortrows(outputRanks(:,3:end,ii),'MissingPlacement','last');
        [~,outputRanks(:,1,ii)] = sort(ordering);
    end
end
% -------------------------------------------------------------------------
function [portfolio,outputRanks] = psf(inputRanks,maxCardinality)
    outputRanks = inputRanks;
    [nalgos,~,ninst] = size(inputRanks);
    algorithmSet = 1:nalgos;
    [~,averageRank] = sortrows(nanmean(inputRanks(:,3:end,:),3),'ascend');
    current = averageRank(1);
    portfolio = zeros(1,maxCardinality);
    portfolio(1) = algorithmSet(current);
%     for ii=2:maxCardinality
%         algorithmSet = setdiff(algorithmSet,portfolio);
%         inputRanks(current,:,:) = [];
%         [~,averageRank] = sortrows(nanmean(inputRanks(:,3:end,:),3),'ascend');
%         current = averageRank(1);
%         portfolio(ii) = algorithmSet(current);
%     end
    for ii=2:maxCardinality
        algorithmSet = setdiff(algorithmSet,portfolio);
        costs = zeros(length(algorithmSet),2);
        for jj=1:length(algorithmSet)
            aux = [algorithmSet(jj) setdiff(portfolio,0)];
            currentRanks = outputRanks(aux,3:end,:);
            for kk=1:ninst
                currentRanks(:,:,kk) = sortrows(currentRanks(:,:,kk),'MissingPlacement','last');
            end
            costs(jj,:) = nanmean(squeeze(currentRanks(1,:,:)),2);
        end
        [~,averageRank] = sortrows(costs,'ascend');
        current = averageRank(1);
        portfolio(ii) = algorithmSet(current);
    end
    portfolio = sort(portfolio);
    outputRanks(setdiff(1:nalgos,portfolio),:,:) = [];
    for ii=1:ninst
        [~,ordering] = sortrows(outputRanks(:,3:end,ii),'MissingPlacement','last');
        [~,outputRanks(:,1,ii)] = sort(ordering);
    end
end
% -------------------------------------------------------------------------
function [portfolio,outputRanks] = psb(inputRanks,maxCardinality)
    outputRanks = inputRanks;
    [nalgos,~,ninst] = size(inputRanks);
    portfolio = 1:nalgos;
%     [~,averageRank] = sortrows(nanmean(inputRanks(:,3:end,:),3),'descend');
%     current = averageRank(1);
%     portfolio = setdiff(portfolio,portfolio(current));
%     while length(portfolio)>maxCardinality
%         inputRanks(current,:,:) = [];
%         [~,averageRank] = sortrows(nanmean(inputRanks(:,3:end,:),3),'descend');
%         current = averageRank(1);
%         portfolio = setdiff(portfolio,portfolio(current));
%     end
    while length(portfolio)>maxCardinality
        costs = zeros(length(portfolio),2);
        for jj=1:length(portfolio)
            aux = setdiff(portfolio,jj);
            currentRanks = outputRanks(aux,3:end,:);
            for kk=1:ninst
                currentRanks(:,:,kk) = sortrows(currentRanks(:,:,kk),'MissingPlacement','last');
            end
            costs(jj,:) = nanmean(squeeze(currentRanks(1,:,:)),2);
        end
        [~,averageRank] = sortrows(costs,'descend');
        current = averageRank(1);
        portfolio = setdiff(portfolio,portfolio(current));
    end
    % disp(['They are equal if 1 = ' num2str(all(portfolioA==portfolio))]);
    outputRanks(setdiff(1:nalgos,portfolio),:,:) = [];
    for ii=1:ninst
        [~,ordering] = sortrows(outputRanks(:,3:end,ii),'MissingPlacement','last');
        [~,outputRanks(:,1,ii)] = sort(ordering);
    end
end
% -------------------------------------------------------------------------
function [portfolio,outputRanks] = randps(inputRanks,maxCardinality)
    [nalgos,~,ninst] = size(inputRanks);
    portfolio = sort(randperm(nalgos,maxCardinality));
    outputRanks = inputRanks(portfolio,:,:);
    for ii=1:ninst
        [~,ordering] = sortrows(outputRanks(:,3:end,ii),'MissingPlacement','last');
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
function H = getEntropy(Ybin)
    Pgood = mean(Ybin,1);
    aux = Pgood.*log(Pgood);
    aux(isnan(aux)) = 0;
    H = -sum(aux);
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