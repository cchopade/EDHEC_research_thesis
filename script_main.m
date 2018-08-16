%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for Thesis
% Title: Implications of Stock-Bond time varying correlationon portfolio diversification
% CHOPADE Chinmay (52330)
% MSc Financial Markets(2016 - 17), EDHEC Business School - Nice, France
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clean up and adding path

clear all % delete all variables
close all % close any charts
clc % take command window cursor to top
cd('/home/chinmay/SpiderOak Hive/Archives/MSc/SRT/Thesis/data/new')
addpath(genpath('/home/chinmay/SpiderOak Hive/Archives/MSc/SRT/Thesis/matlab'))

%addpath('C:\Users\Chinmay\Google Drive\Archives\MSc\SRT\Thesis\matlab')
%addpath(genpath('C:\Users\Chinmay\Google Drive\Archives\MSc\SRT\Thesis\matlab\MFE'))
%addpath('C:\Users\Chinmay\Google Drive\Archives\MSc\SRT\Thesis\matlab\svar')
%cd('C:\Users\Chinmay\Google Drive\Archives\MSc\SRT\Thesis\data\new')
%addpath('E:\SRT\thesis\matlab')
%addpath('E:\SRT\thesis\matlab\MFE')
%cd('E:\SRT\thesis\data\new')

%% Import data

% data for 36 month rolling window
[prices,~,everything]=xlsread ('newdata1.xlsx','data36');
dates=datenum(everything(2:end,1),'dd-mm-yyyy'); % convert dates to number format
colnames=everything(1,2:3); % extract column names
monthlydata36=[dates prices(:,1:2)]; %monthly data matrix to be used
rf36=prices(2:end,3); % risk free rate
monthlyret36=log(monthlydata36(2:end,2:end)./monthlydata36(1:end-1,2:end)); % monthly log returns of time series

% data for 60 month rolling window
[prices,~,everything]=xlsread ('newdata1.xlsx','data60');
dates=datenum(everything(2:end,1),'dd-mm-yyyy'); % convert dates to number format
colnames=everything(1,2:3); % extract column names
monthlydata60=[dates prices(:,1:2)]; %monthly data matrix to be used
rf60=prices(2:end,3); % risk free rate
monthlyret60=log(monthlydata60(2:end,2:end)./monthlydata60(1:end-1,2:end)); % monthly log returns of time series

clear prices;clear everything;clear dates; % clear variables no longer reauired

%% Summary statistics and graphs of complete data

[~,~,summary]=FNsumstatsNEW(monthlyret60,colnames);
summary
figure();
plot(monthlydata60(2:end,1),monthlyret60);
datetick('x','yyyy');
title('Graph of Monthly Log Returns over the period of sample data');
xlabel('Year');
ylabel('Monthly Log Returns');
legend(' S&P 500',' 10 year US treasury bond');

% plot prices data for stocks and bonds

figure();
x=monthlydata60(2:end,1);
y1=cumprod(1.+ monthlyret60(:,1));
y2=cumprod(1.+ monthlyret60(:,2));
%y1=monthlydata(2:end,2);
%y2=monthlydata(2:end,3);
%yyaxis left
subplot(2,1,1);
plot(x,y1);
datetick('x','yyyy');
%yyaxis right
subplot(2,1,2);
plot(x,y2);
datetick('x','yyyy');

%% Correlaion analysis

constantcorr36=corrcoef(monthlyret36(:,1),monthlyret36(:,2)); % correlation over entier period

% rolling window 36 months
for i=37:size(monthlyret36,1)
    
    r=corrcoef(monthlyret36(i-36:i-1,1),monthlyret36(i-36:i-1,2));
    m36correlation(i-36,1)=r(1,2);
    clear r;
end

constantcorr60=corrcoef(monthlyret60(:,1),monthlyret60(:,2)); % correlation over entier period

% rolling window 60 months
for i=61:size(monthlyret60,1)
    
    r=corrcoef(monthlyret60(i-60:i-1,1),monthlyret60(i-60:i-1,2));
    m60correlation(i-60,1)=r(1,2);
    clear r;
end

% make matrices for plotting the correlations
constantcorr36m=repmat(constantcorr36(1,2),size(m36correlation,1),1);
plotcorr36=[m36correlation constantcorr36m];

constantcorr60m=repmat(constantcorr60(1,2),size(m60correlation,1),1);
plotcorr60=[m60correlation constantcorr60m];

% plot the correlations
figure();
subplot(2,1,1);
plot(monthlydata36(38:end,1),plotcorr36);

indexmin=find(min(plotcorr36(:,1))==plotcorr36(:,1));
xmin=monthlydata36(indexmin,1);
ymin=plotcorr36(indexmin,1);
strmin=['Min= ',num2str(ymin)];
text(xmin,ymin,strmin,'HorizontalAlignment','left');

indexmax=find(max(plotcorr36(:,1))==plotcorr36(:,1));
xmax=monthlydata36(indexmax,1);
ymax=plotcorr36(indexmax,1);
strmax=['Max= ',num2str(ymax)];
text(xmax,ymax,strmax,'HorizontalAlignment','left');

datetick('x','yyyy','keeplimits')
title('Stock Bond Correlation');
xlabel('Year');
ylabel('Correlation Coefficient');
legend(' 36 month rolling window',' Correlation over entire period');

subplot(2,1,2);
plot(monthlydata60(62:end,1),plotcorr60);

indexmin=find(min(plotcorr60(:,1))==plotcorr60(:,1));
xmin=monthlydata60(indexmin,1);
ymin=plotcorr60(indexmin,1);
strmin=['Min= ',num2str(ymin)];
text(xmin,ymin,strmin,'HorizontalAlignment','left');

indexmax=find(max(plotcorr60(:,1))==plotcorr60(:,1));
xmax=monthlydata60(indexmax,1);
ymax=plotcorr60(indexmax,1);
strmax=['Max= ',num2str(ymax)];
text(xmax,ymax,strmax,'HorizontalAlignment','left');

datetick('x','yyyy','keeplimits')
title('Stock Bond Correlation');
xlabel('Year');
ylabel('Correlation Coefficient');
legend(' 60 month rolling window',' Correlation over entire period');

%% Testing stationarity
% PP test for trend stationarity
[pph,pppvalue,ppstat,ppcvalue,ppreg]=pptest(monthlyret60(:,1));
[pphb,pppvalueb,ppstatb,ppcvalueb,ppregb]=pptest(monthlyret60(:,2));

% ADF test for unit root
[adfh,adfpvalue,adfstat,adfcvalue,adfreg]=adftest(monthlyret60(:,1),'model','TS');
[adfhb,adfpvalueb,adfstatb,adfcvalueb,adfregb]=adftest(monthlyret60(:,2),'model','TS');

col={'Stationarity Test';'PP';'ADF'};
col1={'SP500';pph;adfh};
col2={'USG 10yr bond';pphb;adfhb};
[col col1 col2]

%% Testing for auto correlatin visually
figure('Name','Autocorrelation test for Stock Returns');
subplot(2,1,1)
autocorr(monthlyret60(:,1));
title('S&P500 Returns');
subplot(2,1,2)
autocorr(monthlyret60(:,1).^2);
title('S&P500 Squared Returns');


figure('Name','Autocorrelation test for Bond Returns');
subplot(2,1,1)
autocorr(monthlyret60(:,2));
title('US Govt. 10yr bond Returns');
subplot(2,1,2)
autocorr(monthlyret60(:,2).^2);
title('US Govt. 10yr bond Squared Returns');

%% Testing for ARCH Effect

residuals=monthlyret60(:,1)-mean(monthlyret60(:,1));
harchsp500=archtest(residuals.^2);
[hlbqsp500,pvals,stats,cvals]=lbqtest(residuals.^2);

clear residuals;
residuals=monthlyret60(:,2)-mean(monthlyret60(:,2));
harchus10=archtest(residuals.^2);
[hlbqus10,pvalb,statsb,cvalb]=lbqtest(residuals.^2);

col={'ARCH Effect Test';'ARCH LM';'Ljung Box Q'};
col1={'SP500';harchsp500;hlbqsp500};
col2={'USG 10yr bond';harchus10;hlbqus10};
[col col1 col2]

%% 36 month rolling window portfolios

pwtgmv36=zeros(size(monthlyret36,1)-36,2);
pmeangmv36=zeros(size(monthlyret36,1)-36,1);

pwtmsr36=zeros(size(monthlyret36,1)-36,2);
pmeanmsr36=zeros(size(monthlyret36,1)-36,1);

pwtew36=ones(size(monthlyret36,1)-36,2).*0.5;
pmeanew36=zeros(size(monthlyret36,1)-36,1);

for i=37:size(monthlyret36,1)
    %set up portfolio object
    rf=mean(rf36(i-36:i-1,:));
    p=Portfolio('Name','SPX','RiskFreeRate',rf);
    p=p.setAssetList(colnames);
    p=p.estimateAssetMoments(monthlyret36(i-36:i-1,:));
    p=p.setDefaultConstraints;
    

    % estimate GMV portfolio
    [pwvt,pbuy,psell]=p.estimateFrontier();
    %[prsk,pret]=p.estimatePortMoments(pwvt); 
    pwtgmv36(i-36,:)=pwvt(:,1)'; %weights of stock and bond
    w=pwvt(:,1);
    R=monthlyret36(i,:)';
    pmeangmv36(i-36,1)=w'*R; % out of sample portfolio return
    clear prsk;clear pret;clear w;clear R;
    
    % estimate MSR portfolio
    pwt = estimateMaxSharpeRatio(p);
    %[prsk, pret] = estimatePortMoments(p, pwt);
    pwtmsr36(i-36,:)=pwt'; %weights of stock and bond
    R=monthlyret36(i,:)';
    pmeanmsr36(i-36,1)=pwt'*R; % out of sample portfolio return
    
    % return for EW portoflio
    R=monthlyret36(i,:)';
    pwt=[0.5;0.5];
    pmeanew36(i-36,1)=pwt'*R;;
    clear p;clear w;clear R;clear pwt;
        
end

plot(monthlydata36(38:end,1),[pmeangmv36 pmeanmsr36 pmeanew36]);
datetick('x','yyyy');
title('Portfolio Returns');
xlabel('Year');
ylabel('Returns');
legend(' Minimum Variance Portfolio',' Maximum Sharpe Ratio Portfolio',' Equally Weighted Portfolio');

% calculate volatility and sharpe ratio for the entire period
pvolgmv36full=std(pmeangmv36);
pvolmsr36full=std(pmeanmsr36);
pvolew36full=std(pmeanew36);

psharpegmv36full=mean(pmeangmv36-rf36(37:end,1))/pvolgmv36full;
psharpemsr36full=mean(pmeanmsr36-rf36(37:end,1))/pvolmsr36full;
psharpeew36full=mean(pmeanew36-rf36(37:end,1))/pvolew36full;

col1={'Portfolios';'GMV';'MSR';'EW'};
col2={'Sharpe Ratio';psharpegmv36full;psharpemsr36full;psharpeew36full};
[col1 col2]

% calculate evolution of sharpe ratio and compare with correlation
pcorr36roll=zeros(size(pmeangmv36,1)-60,1); % variable for correlatio on sub periods

pvolgmv36roll=zeros(size(pmeangmv36,1)-60,1);
psharpegmv36roll = zeros(size(pmeangmv36,1)-60,1);

pvolmsr36roll=zeros(size(pmeanmsr36,1)-60,1);
psharpemsr36roll = zeros(size(pmeanmsr36,1)-60,1);

pvolew36roll=zeros(size(pmeanew36,1)-60,1);
psharpeew36roll = zeros(size(pmeanew36,1)-60,1);


% calculate volatility of out of sample returns in blocks of 36 months
for i=1:size(pmeangmv36,1)-60
    
    % Correlation
    r=corrcoef(monthlyret36(36+i:95+i,1),monthlyret36(36+i:95+i,2));
    pcorr36roll(i,1)=r(1,2);
    
    % GMV portfolios
    pvolgmv36roll(i,1)=std(pmeangmv36(i:i+59,:));
    pexcess=pmeangmv36(i:i+59,:)-rf36(36+i:95+i,:);
    psharpegmv36roll(i,1)=mean(pexcess)/pvolgmv36roll(i,1);
    clear pexcess;
    
    % MSR portfolios
    pvolmsr36roll(i,1)=std(pmeanmsr36(i:i+59,:));
    pexcess=pmeanmsr36(i:i+59,:)-rf36(36+i:95+i,:);
    psharpemsr36roll(i,1)=mean(pexcess)/pvolmsr36roll(i,1);
    clear pexcess;
    
    % EW portfolios
    pvolew36roll(i,1)=std(pmeanew36(i:i+59,:));
    pexcess=pmeanew36(i:i+59,:)-rf36(36+i:95+i,:);
    psharpeew36roll(i,1)=mean(pexcess)/pvolew36roll(i,1);
end


figure();
plot(monthlydata36(38:end-60,1),[pcorr36roll psharpegmv36roll psharpemsr36roll psharpeew36roll]);
datetick('x','yyyy');
title('Sharpe Ratio v/s Correlation');
xlabel('Year');
ylabel('Shaarpe Ratio / Correlation');
legend(' Correlation',' Minimum Variance Portfolio',' Maximum Sharpe Ratio Portfolio', ' Equally Weighted Portfolio');


%%

plot(monthlydata36(98:end,1),[pcorr36roll psharpegmv36roll psharpemsr36roll psharpeew36roll]);
datetick('x','yyyy');
title('Sharpe Ratio v/s Correlation');
xlabel('Year');
ylabel('Shaarpe Ratio / Correlation');
legend(' Correlation',' Minimum Variance Portfolio',' Maximum Sharpe Ratio Portfolio', ' Equally Weighted Portfolio');


%% 60 month rolling window portfolios

pwtgmv60=zeros(size(monthlyret60,1)-60,2);
pmeangmv60=zeros(size(monthlyret60,1)-60,1);

pwtmsr60=zeros(size(monthlyret60,1)-60,2);
pmeanmsr60=zeros(size(monthlyret60,1)-60,1);

pwtew60=ones(size(monthlyret60,1)-60,2).*0.5;
pmeanew60=zeros(size(monthlyret60,1)-60,1);

for i=61:size(monthlyret60,1)
    %set up portfolio object
    rf=mean(rf60(i-60:i-1,:));
    p=Portfolio('Name','SPX','RiskFreeRate',rf);
    p=p.setAssetList(colnames);
    p=p.estimateAssetMoments(monthlyret60(i-60:i-1,:));
    p=p.setDefaultConstraints;
    

    % estimate GMV portfolio
    [pwvt,pbuy,psell]=p.estimateFrontier();
    %[prsk,pret]=p.estimatePortMoments(pwvt); 
    pwtgmv60(i-60,:)=pwvt(:,1)'; %weights of stock and bond
    w=pwvt(:,1);
    R=monthlyret60(i,:)';
    pmeangmv60(i-60,1)=w'*R; % out of sample portfolio return
    clear prsk;clear pret;clear w;clear R;
    
    % estimate MSR portfolio
    pwt = estimateMaxSharpeRatio(p);
    %[prsk, pret] = estimatePortMoments(p, pwt);
    pwtmsr60(i-60,:)=pwt'; %weights of stock and bond
    R=monthlyret60(i,:)';
    pmeanmsr60(i-60,1)=pwt'*R; % out of sample portfolio return
    
    % return for EW portoflio
    R=monthlyret60(i,:)';
    pwt=[0.5;0.5];
    pmeanew60(i-60,1)=pwt'*R;;
    clear p;clear w;clear R;clear pwt;
        
end

plot(monthlydata60(62:end,1),[pmeangmv60 pmeanmsr60 pmeanew60]);
datetick('x','yyyy');
title('Portfolio Returns');
xlabel('Year');
ylabel('Returns');
legend(' Minimum Variance Portfolio',' Maximum Sharpe Ratio Portfolio',' Equally Weighted Portfolio');

% calculate volatility and sharpe ratio for the entire period
pvolgmv60full=std(pmeangmv60);
pvolmsr60full=std(pmeanmsr60);
pvolew60full=std(pmeanew60);

psharpegmv60full=mean(pmeangmv60-rf60(61:end,1))/pvolgmv60full;
psharpemsr60full=mean(pmeanmsr60-rf60(61:end,1))/pvolmsr60full;
psharpeew60full=mean(pmeanew60-rf60(61:end,1))/pvolew60full;

col1={'Portfolios';'GMV';'MSR';'EW'};
col2={'Sharpe Ratio';psharpegmv60full;psharpemsr60full;psharpeew60full};
[col1 col2]

% calculate evolution of sharpe ratio and compare with correlation
pcorr60roll=zeros(size(pmeangmv60,1)-60,1); % variable for correlatio on sub periods

pvolgmv60roll=zeros(size(pmeangmv60,1)-60,1);
psharpegmv60roll = zeros(size(pmeangmv60,1)-60,1);

pvolmsr60roll=zeros(size(pmeanmsr60,1)-60,1);
psharpemsr60roll = zeros(size(pmeanmsr60,1)-60,1);

pvolew60roll=zeros(size(pmeanew60,1)-60,1);
psharpeew60roll = zeros(size(pmeanew60,1)-60,1);


% calculate volatility of out of sample returns in blocks of 36 months
for i=1:size(pmeangmv60,1)-60
    
    % Correlation
    r=corrcoef(monthlyret60(60+i:119+i,1),monthlyret60(60+i:119+i,2));
    pcorr60roll(i,1)=r(1,2);
    
    % GMV portfolios
    pvolgmv60roll(i,1)=std(pmeangmv60(i:i+59,:));
    pexcess=pmeangmv60(i:i+59,:)-rf60(60+i:119+i,:);
    psharpegmv60roll(i,1)=mean(pexcess)/pvolgmv60roll(i,1);
    clear pexcess;
    
    % MSR portfolios
    pvolmsr60roll(i,1)=std(pmeanmsr60(i:i+59,:));
    pexcess=pmeanmsr60(i:i+59,:)-rf60(60+i:119+i,:);
    psharpemsr60roll(i,1)=mean(pexcess)/pvolmsr60roll(i,1);
    clear pexcess;
    
    % EW portfolios
    pvolew60roll(i,1)=std(pmeanew60(i:i+59,:));
    pexcess=pmeanew60(i:i+59,:)-rf60(60+i:119+i,:);
    psharpeew60roll(i,1)=mean(pexcess)/pvolew60roll(i,1);
end

figure();
plot(monthlydata60(122:end,1),[pcorr60roll psharpegmv60roll psharpemsr60roll psharpeew60roll]);
datetick('x','yyyy');
title('Sharpe Ratio v/s Correlation');
xlabel('Year');
ylabel('Shaarpe Ratio / Correlation');
legend(' Correlation',' Minimum Variance Portfolio',' Maximum Sharpe Ratio Portfolio', ' Equally Weighted Portfolio');


%% extras

%summary statistics and histograms of correlations
correlations=[m36correlation m60correlation];
colnames={'36 month correlation','60 month correlation'};
[~,~,summarycorr]=FNsumstatsNEW(correlations,colnames)

figure();
subplot(2,1,1);
hist(m36correlation);
xlabel('36 month rolling window correlation');
ylabel('Frequency');
title('Distribution of 36 month rolling window correlation');

subplot(2,1,2);
hist(m60correlation);
xlabel('60 month rolling window correlation');
ylabel('Frequency');
title('Distribution of 60 month rolling window correlation');

% stationarity of correlation time series

[pph36,pppvalue36,ppstat36,ppcvalue36,ppreg36]=pptest(m36correlation);
[pph60,pppvalue60,ppstat60,ppcvalue60,ppreg60]=pptest(m60correlation);

% ADF test for unit root
[adfh36,adfpvalue36,adfstat36,adfcvalue36,adfreg36]=adftest(m36correlation,'model','TS');
[adfh60,adfpvalue60,adfstat60,adfcvalue60,adfreg60]=adftest(m60correlation,'model','TS');

col={'Stationarity Test (Correlation)';'36 month';'60 month'};
col1={'PP Statistic';ppstat36;ppstat60};
col2={'p-value';pppvalue36;pppvalue60};
col3={'ADF Statistic';adfstat36;adfstat60};
col4={'p-value';adfpvalue36;adfpvalue60};

[col col1 col2 col3 col4]

%% Autocorrelation for correlation series
figure('Name','Autocorrelation test for 36 month correlation');
subplot(2,1,1)
autocorr(m36correlation);
title('36 month window correlation');
subplot(2,1,2)
autocorr(m60correlation);
title('60 month window correlation');

%% scatter plots for 36 month correlations and sharpe

figure();
scatter(m36correlation(61:end,1),psharpegmv36roll)
title('Sharpe ratio (Min. Variance) v/s Correlation')
xlabel('Correlation');
ylabel('Sharpe Ratio (Minimum Variance)');

figure();
scatter(m36correlation(61:end,1),psharpemsr36roll)
title('Sharpe ratio (Max. Sharpe) v/s Correlation')
xlabel('Correlation');
ylabel('Sharpe Ratio (Maximum Sharpe ratio)');

figure();
scatter(m36correlation(61:end,1),psharpeew36roll)
title('Sharpe ratio (Equally Weighted) v/s Correlation')
xlabel('Correlation');
ylabel('Sharpe Ratio (Equally Weighted)');

%% scatter plots for 60 month correlation and sharpe

figure();
scatter(m60correlation(61:end,1),psharpegmv60roll)
title('Sharpe ratio (Min. Variance) v/s Correlation')
xlabel('Correlation');
ylabel('Sharpe Ratio (Minimum Variance)');

figure();
scatter(m60correlation(61:end,1),psharpemsr60roll)
title('Sharpe ratio (Max. Sharpe) v/s Correlation')
xlabel('Correlation');
ylabel('Sharpe Ratio (Maximum Sharpe ratio)');

figure();
scatter(m60correlation(61:end,1),psharpeew60roll)
title('Sharpe ratio (Equally Weighted) v/s Correlation')
xlabel('Correlation');
ylabel('Sharpe Ratio (Equally Weighted)');

%% parametric testing of sharpe ratio times series

[h1,p1,ci1,stats1]=ttest(psharpegmv36roll,psharpemsr36roll);
[h2,p2,ci2,stats2]=ttest(psharpegmv36roll,psharpeew36roll);
[h3,p3,ci3,stats3]=ttest(psharpemsr36roll,psharpeew36roll);

[h4,p4,ci4,stats4]=ttest(psharpegmv60roll,psharpemsr60roll);
[h5,p5,ci5,stats5]=ttest(psharpegmv60roll,psharpeew60roll);
[h6,p6,ci6,stats6]=ttest(psharpemsr60roll,psharpeew60roll);

%% regressions with sarpe ratio

%36 month 
xcorr36 = m36correlation(61:end,1);
ygmv36 = psharpegmv36roll;
ymsr36=psharpemsr36roll;
yew36=psharpeew36roll;

%60 month
xcorr60 = m60correlation(61:end,1);
ygmv60 = psharpegmv60roll;
ymsr60=psharpemsr60roll;
yew60=psharpeew60roll;

% linear model
mdlsharpegmv36=LinearModel.fit(xcorr36,ygmv36);
[mdlsharpegmv36r,stats]=robustfit(xcorr36,ygmv36);
mdlsharpemsr36=LinearModel.fit(xcorr36,ymsr36);
mdlsharpeew36=LinearModel.fit(xcorr36,yew36);

mdlsharpegmv60=LinearModel.fit(xcorr60,ygmv60);
mdlsharpemsr60=LinearModel.fit(xcorr60,ymsr60);
mdlsharpeew60=LinearModel.fit(xcorr60,yew60);

% non linear model
modelfun = @(b,x)b(1) + b(2)*x(:,1) + b(3)*x(:,1).^2; 
beta0=[2 2 2];

mdlsharpegmvsq36 = NonLinearModel.fit(xcorr36,ygmv36,modelfun,beta0);
mdlsharpemsrsq36 = NonLinearModel.fit(xcorr36,ymsr36,modelfun,beta0);
mdlsharpeewsq36 = NonLinearModel.fit(xcorr36,yew36,modelfun,beta0);

mdlsharpegmvsq60 = NonLinearModel.fit(xcorr60,ygmv60,modelfun,beta0);
mdlsharpemsrsq60 = NonLinearModel.fit(xcorr60,ymsr60,modelfun,beta0);
mdlsharpeewsq60 = NonLinearModel.fit(xcorr60,yew60,modelfun,beta0);

%% regression for stock weights

xwcorr36 = m36correlation;
ywgmv36=pwtgmv36(:,1);
ywmsr36=pwtmsr36(:,1);
ywew36=pwtew36(:,1);

xwcorr60 = m60correlation;
ywgmv60=pwtgmv60(:,1);
ywmsr60=pwtmsr60(:,1);
ywew60=pwtew60(:,1);


mdlswtgmv36=LinearModel.fit(xwcorr36,ywgmv36);
mdlswtmsr36=LinearModel.fit(xwcorr36,ywmsr36);
mdlswtew36=LinearModel.fit(xwcorr36,ywew36);

mdlswtgmv60=LinearModel.fit(xwcorr60,ywgmv60);
mdlswtmsr60=LinearModel.fit(xwcorr60,ywmsr60);
mdlswtew60=LinearModel.fit(xwcorr60,ywew60);

% non linear model
modelfun = @(b,x)b(1) + b(2)*x(:,1) + b(3)*x(:,1).^2; 
beta0=[2 2 2];

mdlswtgmvsq36 = NonLinearModel.fit(xwcorr36,ywgmv36,modelfun,beta0);
mdlswtmsrsq36 = NonLinearModel.fit(xwcorr36,ywmsr36,modelfun,beta0);
mdlswtewsq36 = NonLinearModel.fit(xwcorr36,ywew36,modelfun,beta0);

mdlswtgmvsq60 = NonLinearModel.fit(xwcorr60,ywgmv60,modelfun,beta0);
mdlswtmsrsq60 = NonLinearModel.fit(xwcorr60,ywmsr60,modelfun,beta0);
mdlswtewsq60 = NonLinearModel.fit(xwcorr60,ywew60,modelfun,beta0);

%% other risk measures

% Max Drawdown
gmv36price=ret2price(pmeangmv36,100);
msr36price=ret2price(pmeanmsr36,100);
ew36price=ret2price(pmeanew36,100);

gmv60price=ret2price(pmeangmv60,100);
msr60price=ret2price(pmeanmsr60,100);
ew60price=ret2price(pmeanew60,100);

gmv36dd=maxdrawdown(gmv36price);
msr36dd=maxdrawdown(msr36price);
ew36dd=maxdrawdown(ew36price);

gmv60dd=maxdrawdown(gmv60price);
msr60dd=maxdrawdown(msr60price);
ew60dd=maxdrawdown(ew60price);

col1={'Portfolios';'GMV 36';'MSR 36';'EW 36';'GMV 60';'MSR 60';'EW 60'};
col2={'Max DD';gmv36dd;msr36dd;ew36dd;gmv60dd;msr60dd;ew60dd};
[col1 col2]


% Tracking error

[irgmv36, tegmv36] = inforatio(pmeangmv36, pmeanew36);
[irmsr36, temsr36] = inforatio(pmeanmsr36, pmeanew36);

[irgmv60, tegmv60] = inforatio(pmeangmv60, pmeanew60);
[irmsr60, temsr60] = inforatio(pmeanmsr60, pmeanew60);

col1={'Portfolios';'GMV 36';'MSR 36';'GMV 60';'MSR 60'};
col2={'Tracking Error';tegmv36;temsr36;tegmv60;temsr60};
col3={'Information Ratio';irgmv36;irmsr36;irgmv60;irmsr60};
[col1 col2 col3]

   
 %% mfe toolbox
 
 data = monthlyret36;
 [PARAMETERS,LL,HT,VCV,SCORES,DIAGNOSTICS] = dcc(data,[],1,0,1);
 
 cc=zeros(size(HT,3),1);
 for i=1:size(HT,3)
     
    sigma=HT(:,:,i);
    dsigma=diag(sigma);
    sdsigma=sqrt([dsigma(1),0;0,dsigma(2)]);
    cd=inv(sdsigma)*sigma*inv(sdsigma);
    conditionalcorr36(i,1)=cd(1,2);
    clear emp; clear temp;
 end
 
plot(monthlydata36(38:end,1),[m36correlation m60correlation conditionalcorr36(37:end,1)])
datetick('x','yyyy','keeplimits')
legend('36 month RW Correlation','60 month RW Correlation','DCC GARCH(1,1) conditional correlation');
title('Comparison of Rolling Window and Conditional Correlation');
xlabel('Year');
ylabel('Correlation');

[~,~,summary]=FNsumstatsNEW(conditionalcorr36);

%%
figure();
subplot(3,1,1);
hist(m36correlation);
xlabel('36 month rolling window correlation');
ylabel('Frequency');
title('Distribution of 36 month rolling window correlation');

subplot(3,1,2);
hist(m60correlation);
xlabel('60 month rolling window correlation');
ylabel('Frequency');
title('Distribution of 60 month rolling window correlation');

subplot(3,1,3);
hist(conditionalcorr36);
ylim([0 150]);
xlabel('Conditional Correlation');
ylabel('Frequency');
title('Distribution of conditional correlation');

%%

mdl=garch(1,1);
[EstMdls,EstParamCovs,logLs,infos]=estimate(mdl,monthlyret36(:,1));
[EstMdlb,EstParamCovb,logLb,infob]=estimate(mdl,monthlyret36(:,2));
