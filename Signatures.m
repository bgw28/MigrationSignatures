%% User Parameters: 
dataPath = "./netmigration.csv";    % file path to data file
exportFigures = true;               % save figures?


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

% compute singular values for each county on its 6 decades
% divide the square value by 6 (the number of decades)
singularValue = zeros(size(full,1),1);
for i = 1:size(full,1)
    vectors = permute(full(i,:,:),[3,2,1]);
    singularValue(i) = svds(vectors,1)^2 /6;
end

% plot the distribution
fig = figure;
histogram(singularValue,linspace(0,1,26))
xlim([0,1])
set(gca,'fontname','SansSerif')
xlabel("Mean Squared Cosine Similarity")
ylabel("Number of Counties")
fontsize(fig, 15, "points")

if exportFigures
    exportgraphics(fig,"./Figures/ConsistencyOverTime.png",Resolution=300)
    exportgraphics(fig,"./Figures/ConsistencyOverTime.eps")
end

clear full fig singularValue vectors



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

signatures = zeros(length(mainClusters),16);
for i = 1:length(mainClusters)  % for each cluster above the min size
    vectors = unitRateVectors(labels==mainClusters(i),:);  % take the unit rate vectors
    [a,~,s] = svds(vectors,1);  % find their singular vector
    signatures(i,:) = s .* sign(mean(a));  % fix the sign for consistency
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
    % title(names(i))
    fontsize(fig, 20, "points")

    if exportFigures
        exportgraphics(fig,['./Figures/',clusterNames{i},'.png'],resolution=300)
        exportgraphics(fig,['./Figures/',clusterNames{i},'.eps'],resolution=300)
    end
end
clear colors i fig

%% Type Scores

scores = normalize(reshape(permute(dataMatrix(:,:,:,1),[3 1 2]), [], 16),2,'norm') * signatures';

t = cell2table(num2cell(scores));
t.Properties.VariableNames = [clusterNames{:}];

stateName = convertCharsToStrings(reshape(repmat(stateNames,1,size(dataMatrix,3))',[],1));

countyName = convertCharsToStrings(reshape(repmat(countyNames,1,size(dataMatrix,3))',[],1));

fips = reshape(repmat(fipsNum,1,size(dataMatrix,3))',[],1);

decades = ["1960";"1970";"1980";"1990";"2000";"2010"];
decade = reshape(repmat(decades,1,size(dataMatrix,1)),[],1);

scoresTable = [table(stateName,countyName,decade,fips), t];
clear scores t decades stateName countyName fips decade