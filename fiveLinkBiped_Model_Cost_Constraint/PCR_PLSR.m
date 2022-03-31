% Load data and read different vectors and matrices
Data= load('Fukuchi_Features_DTW_RMS.mat');

Cost_Components_All = Data.Cost_Components;
RMS_Performance     = Data.RMSE_Sum;


% CV for PCR MSEP
X = Cost_Components_All;
y = RMS_Performance;

[Xl,Yl,Xs,Ys,beta,pctVar,PLSmsep] = plsregress(X,y,10,'CV',10);

n = size(X,1);
PCRmsep = sum(crossval(@pcrsse,X,y,'KFold',10),1) / n;

plot(0:10,PLSmsep(2,:),'b-o',0:10,PCRmsep,'r-^');
xlabel('Number of components');
ylabel('Estimated Mean Squared Prediction Error');
legend({'PLSR' 'PCR'},'location','NE');


% PCR = PCA first + regression
[PCALoadings,PCAScores,PCAVar] = pca(normalize(X));
betaPCR = regress(y-mean(y), PCAScores(:,1:2));
%To make the PCR results easier to interpret in terms of the original spectral data, transform to regression coefficients for the original, uncentered variables.


figure(2)
plot(100*cumsum(PCAVar)/sum(PCAVar),'r-^');
xlabel('Number of Principal Components');
ylabel('Percent Variance Explained in X');
legend({'PLSR' 'PCR'},'location','SE');


betaPCR = PCALoadings(:,1:2)*betaPCR;
betaPCR = [mean(y) - mean(X)*betaPCR; betaPCR];
yfitPCR = [ones(n,1) X]*betaPCR;
% Plot fitted vs. observed response for the PLSR and PCR fits.

% plot(y,yfitPLS,'bo',y,yfitPCR,'r^');
plot(y,yfitPCR,'bo',y,yfitPCR,'r^');
xlabel('Observed Response');
ylabel('Fitted Response');
legend({'PLSR with 2 Components' 'PCR with 2 Components'},  ...
	'location','NW');
