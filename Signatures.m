%% User Parameters: 
dataPath = "./netmigration.csv";    % file path to data file
exportFigures = true;               % save figures?
showTitles = false;

%% Read Data

% read in netmigration data
netMigrationData = readtable(dataPath);
stateNames = netMigrationData{:,1};
countyNames = netMigrationData{:,2};
ages = {"0","5","10","15","20","25","30","35","40","45","50","55","60",...
    "65","70","75"};
fipsNum = netMigrationData{:,3};

% convert table to matrix form and remove 1950's
% county x age group x decade x variable
% variable: {rate, net migration, expected population}
dataMatrix = zeros(height(netMigrationData),16,6,3);
for type = 1:3
    for decade = 6:-1:1
        start = 4+16*(type-1)+3*16*(decade-1);
        dataMatrix(:,:,6-decade+1,type) = netMigrationData{:,start:start+16-1};
    end
end

clear type decade start netMigrationData dataPath


%% Consistency Over Time

% remove counties with any missing rates in 1960-2010
full = dataMatrix(:,:,:,1);
full(any(isnan(full),[2,3]),:,:) = [];

% normalize with respect to magnitude
full = normalize(full,2,'norm');

% compute the maximum possible mean cosine similarity for each county
% across its 6 decades
% compute the mean cosine similarity to the mean normalized migration
% vector
mmcs = zeros(size(full,1),1);
for i = 1:size(full,1)
    vectors = permute(full(i,:,:),[3,2,1]);
    mmcs(i) = mean(vectors * normalize(mean(vectors),2,'norm')');
end

% generate empirical estimate of the theoretical ditribution of singular
% values assuming the migration vectors have a uniform spherical
% distribution

% generate 10,000 random samples, each of 6 unit vectors
% We use a multivariate normal distribution and normalize it to achive the
% desired uniform spherical distrubtion
n = 10000;
sample = mvnrnd(zeros(1,16),eye(16),6*n);
sample = normalize(sample,2,'norm');
sample = permute(reshape(sample,[n,6,16]),[2,3,1]);

% compute their maximum mean cosine similarity
empirical = zeros(n,1);
for i = 1:n
    empirical(i) = mean(sample(:,:,i)*normalize(mean(sample(:,:,i)),2,'norm')');
end

% apply a kernel density estimate
[empiricalDist, evalPoints] = kde(empirical,Bandwidth=0.02);
evalPoints = [0;evalPoints;1];
empiricalDist = [0;empiricalDist;0];


% plot the distribution
fig = figure;
colororder(["#D95319", "#0072BD"])

yyaxis left
p1 = plot(evalPoints,empiricalDist,LineWidth=3);
l = ylim;
ylabel("Probability Density")

yyaxis right
p2 = histogram(mmcs,linspace(0,1,26));
ylabel("Number of Counties")
ylim([0,l(2)*0.04*size(mmcs,1)])

xlim([1/6,1])
set(gca,'fontname','SansSerif')
xlabel("Mean Cosine Similarity")

if showTitles
    title("Continuity of Migration Signatures")
end
fontsize(fig, 15, "points")
legend([p2,p1], ["Counties";"Theoretical Random"])


if exportFigures
    exportgraphics(fig,"./Export/ConsistencyOverTime.png",Resolution=150)
    exportgraphics(fig,"./Export/ConsistencyOverTime.eps")
end

med = median(mmcs)
avg = mean(mmcs)

clear l full fig mmcs vectors evalPoints empiricalDist n sample empirical med avg



%% Data Cleaning
cleanedData = dataMatrix();

minPopulation = 30;
magnitudeThreshold = 0.05;  % lower quantile to discard

% remove state totals
cleanedData(contains(countyNames,"Total"),:,:,:) = [];


% reshape
cleanedData = reshape(permute(cleanedData,[1,3,2,4]),[],16,3);

% remove rows with missing entries
hasNan = any(isnan(cleanedData(:,:,[1,3])),[2,3]);
cleanedData = cleanedData(~hasNan,:,:);

% remove rows if any age group's expected Pop is small
lowPop = any(cleanedData(:,:,3)<=minPopulation,2);

cleanedData = cleanedData(~lowPop,:,:);

% remove the lower percent of magnitudes
magnitude = vecnorm(cleanedData(:,:,1),2,2);
lowMag = magnitude<=quantile(magnitude,magnitudeThreshold);
cleanedData = cleanedData(~lowMag,:,:);

clear minPopulation magnitudeThreshold hasNan lowPop magnitude lowMag


% We only use the migration rates
% data structure is a list of all (post cleaning) migration rate vectors
% in a list without county labels attached
rateVectors = cleanedData(:,:,1);
unitRateVectors = normalize(rateVectors,2,'norm');  % magnitude normalized



%% Clustering

% hierarchical clustering using cosine metric and average center
Z = linkage(rateVectors,"average","cosine");

% minimum number of counties needed in a cluster for it to define a 
% signature
clusterMinSize = 350;

% threshold value of inconsistency value to define clusters from the
% linkage (prefered value achieved through trials not shown here)
inconsistencyThreshold = 8;


% cluster
labels = cluster(Z,"cutoff",inconsistencyThreshold,"depth", ...
    size(rateVectors,1));

% remove small clusters and relabel
l = unique(labels);
counts = arrayfun(@(cl) sum(labels==cl),l);
[~,i] = sort(counts,'descend');
mainClusters = i(1:sum(counts>=clusterMinSize));

clear Z clusterMinSize inconsistencyThreshold l counts i


% we determined that cluster 9 (early career) is important to include
% because it fits counties with a high population density however clusters 
% 6,7,8 while statistically relavent, conflate our signatures and lower 
% the utility of our results. therefore here we remove them while keeping 
% cluster 9
mainClusters(6:8) = [];

%% Compute Signatures

% compute the singular vector for each cluster. The singular vector as
% compared to the mean normalized vector, maximizes the total squared cosine
% similarity rather than just the total cosine similarity. Since there are
% no pairwise negative similarities in each cluster, both methods will give
% similar signature vectors. We choose the singular vector as maximizing
% the squared cosine similarity is less affected by outliers.
signatures = zeros(length(mainClusters),16);
for i = 1:length(mainClusters)  % for each cluster above the min size
    vectors = unitRateVectors(labels==mainClusters(i),:);  % take the unit rate vectors
    [a,~,s] = svds(vectors,1);
    signatures(i,:) = s*sign(mean(a));
end
clear i a s vectors

clusterNames = ["Youth Outmigration";"Exodus";"Midcarrer/Family";"Retirement"; ...
    "Institutional";"Early Career"];

colors = [[.85,.325,.098];[.494,.184,.556];[.466,.674,.188];[0,.447,.741];
    [.929,.694,.125];[.635,.078,.184]];

for i = 1:size(signatures,1)
    fig = figure;
    set(gca,'fontname','SansSerif')
    plot([0,16],[0,0],'--',Linewidth=3, Color='black')
    hold on
    plot(signatures(i,:), LineWidth=4, Color=colors(i,:))
    xlim([1,16])
    xticks(1:2:16)
    xticklabels(ages(1:2:16))
    ylim([-1,1])
    yticks(linspace(-1,1,5))
    xlabel("Age Groups")
    ylabel("Normalized NMR")
    if showTitles
        title(clusterNames(i))
    end
    fontsize(fig, 20, "points")

    if exportFigures
        exportgraphics(fig,['./Export/',erase(clusterNames{i},"/"),'.png'],resolution=150)
        exportgraphics(fig,['./Export/',erase(clusterNames{i},"/"),'.eps'])
    end
end
clear colors i fig

%% Type Scores

% compute the cosine similarity of all county decades to all signature
% types and output the scores in a csv file
scores = normalize(reshape(permute(dataMatrix(:,:,:,1),[3 1 2]), [], 16),2,'norm') * signatures';

% initialize the table
t = cell2table(num2cell(scores));
t.Properties.VariableNames = clusterNames;

% construct the labeling variables
stateName = convertCharsToStrings(reshape(repmat(stateNames,1,size(dataMatrix,3))',[],1));
countyName = convertCharsToStrings(reshape(repmat(countyNames,1,size(dataMatrix,3))',[],1));
fips = reshape(repmat(fipsNum,1,size(dataMatrix,3))',[],1);
decades = ["1960";"1970";"1980";"1990";"2000";"2010"];
decade = reshape(repmat(decades,1,size(dataMatrix,1)),[],1);

% populate the table and export it
scoresTable = [table(stateName,countyName,decade,fips), t];
writetable(scoresTable,"./Export/SignatureTypeCountyScores.csv")

clear scores t decades stateName countyName fips decade

