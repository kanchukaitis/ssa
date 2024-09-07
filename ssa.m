function [RC,A,evec,eval,varexp,mctest,orderedPeriods,noisePeriods,D,C] = ssa(x,M,fspec,itime,iunit,mct,graphics)

% SSA Singular Spectrum Analysis of Time Series
%
%   [RC,A,evec,eval,varexp,mctest,D,C] = ssa(x,M,fspec,itime,iunit,mct,graphics)
%
% Perform SSA (Allen and Smith, 1996, Ghil et al. 2022) and develop confidence intervals
% for SSA selection using Monte Carlo simulations; call with 2, 3, 4, or 5 arguments.
%
% Inputs:
%  x = data (1 x N or N x 1)
%  M = embedding dimension (1 x 1)
%  fspec = Do spectral analysis of eigenvectors? 1=yes, 0=no; Default is no (optional)
%  itime = time vector for plotting the reconstructed components; Default is dimensionless (optional)
%  iunit = time units (e.g. 'years', which is the default) (optional)
%  mct = switch for Monte Carlo testing. 0 = No [default]; 1 = Yes (potentially slow)
%  graphics = switch for graphical diagnosis. 0 = No [default]; 1 = Yes
%
% Outputs:
%  RC = reconstructed components
%  A = Principal Components
%  evec = eigenvectors
%  eval = eigenvalues
%  varexp = fraction variance explained for each eigenvalue
%  mctest = structure containing the c95_val_w,c99_val_ar,c95_val_ar,c90_val_ar, and c50_val_ar confidence intervals from 
%          white and red noise Monte Carlo simulation 
%
% References:
% Allen M., Smith L.A., 1996: Monte Carlo SSA: Detecting irregular 
% oscillations in the presence of coloured noise,  J. Clim.,  9 , 3373-3404. 
%
% Ghil M., R. M. Allen, M. D. Dettinger, K. Ide, D. Kondrashov, M. E. Mann, 
% A. Robertson, A. Saunders, Y. Tian, F. Varadi, and P. Yiou, 2002: 
% Advanced spectral methods for climatic time series, Rev. Geophys.
% 40(1), pp. 3.1-3.41, 10.1029/2000GR000092.
%
% Tsonis, A.A., and J.B. Elsner, 1992: Oscillating global temperature, Nature, 356, 751.
%
% 11/04 - original version (v1)
% 08/08 - modified to run faster in batch mode
% 
% Optionally calls Dave Meko's (dmeko@arizona.edu) specbt1.m routine
% to determine the leading periodicity of the first 4 eigenvectors
% However, this routine requires the Signal Processing Toolbox from Mathworks

if nargin == 2; fspec = 0; end
if nargin <= 3; itime = [1:length(x)]; end 
if isempty(itime); itime = [1:length(x)]; end 
if nargin <= 4; iunit = 'years'; end
if nargin <6; mct = 0; end
if nargin <7; graphics = 0; end
    
[m,N] = size(x);

% make sure its a single row vector
if m > N; x = x';[ m,N] = size(x); end

% get NaN information for calculating and removing the mean
innt = length(x) - sum(isnan(x));

% Calculate the vector means
xtemp = x; xtemp(find(isnan(x))) = 0; Xbar = sum(xtemp)./innt;

% Remove the vector means 
xc = x - Xbar;

% Now get simple first order autocorrelation and variance for Monte Carlo tests
foac = corrcoef(xc(1,1:end-1),xc(1,2:end)); foac = foac(1,2);

foac = 0.67858

disp(['First Order Autocorrelation calculated: ', num2str(foac)])
datavar = var(xc);
disp(['Variance of original datset calculated: ', num2str(datavar)])

% pause

% Prepare the trajectory matrix D
if M > N; error(['Window length greater than series length! M should be less than N']); return; end;
D = NaN .* ones(M,N-M+1);
for i=1:M
  D(i,:)=xc(1,i:N-M+i); 
end                  

% Now, do the covariance matrix calculation;
% First, find NaNs and calculate sample size N (not unbiased following Ghil et al. 2002)
nnan = ~isnan(D); xnan = ind2sub(size(nnan),nnan); xsize = ((xnan*xnan')); 

% Set NaNs in trajectory matrix to zero if there are any
if sum(isnan(D)) > 1; D(find(isnan(D))) = 0; end
% now we can calculate the covariance matrix, again _not_ unbiased following Ghil et al. 2002)
C = (D*D').*(1./xsize);

% Find Eigenvalues and Eigenvectors using Singular Value Decomposition of
% the covariance matrix
[U,S,V] = svd(C);  

% get eigenvalues and fraction variance explained
evec = U;
eval = diag(S); 
varexp = diag(S)/trace(S);

% Make Principal Component time series, A 
A = U'*D;

% Now make the reconstructed components
RC = zeros(m,N);
for i=1:M
    RC(i,:) = conv(U(:,i),A(i,:));
end

% Now, normalize the RCs appropriately 
% (Ghil et al. 2002 with error corrected)
for i=1:N
    if i >= 1 && i <= M-1
        RC(:,i) = RC(:,i)*(1/i);
    elseif i <= N-M+1 && i >= M
        RC(:,i) = RC(:,i)./M;
    else
        RC(:,i) = RC(:,i)/(N-i+1);
    end
end

if mct>0

% disp(['Beginning Monte Carlo Significance Test (White)...'])
%% Monte Carlo based significance test, both white and red noise spectrums
% Gaussian white noise spectrum test
R=randn(m,N);
for i=1:mct

  for j=1:N
      G(1,j) = R(1,j) .* sqrt(datavar);
  end

  for j=1:M
      Dr(j,:)=G(:,j:N-M+j); 
  end

 Cr = (Dr*Dr').*(1./(N-M+1));
 lambda_white(i,:) = diag(U'*Cr*U)';
end

slambda_white = lambda_white'


% clean things up some
clear Dr Cr R G 

% disp(['Beginning Monte Carlo Significance Test (AR1)...'])
% AR(1) model, (naive) red spectrum test
for i=1:mct

 G(1,1) = randn(1,1);

for k=1:N
  G(1,k+1) = (foac*G(1,k)) + randn(1,1);
end

G = G - repmat(mean(G),m,1);
G = G/std(G);

% G = standardize(G)

% for j=1:N
%   G(1,j) = G(1,j) .* sqrt(datavar);
% end


for j=1:M
  Dr(j,:)=G(:,j:N-M+j); 
end

 Cr = (Dr*Dr').*(1./(N-M+1)); % the noise covarince matrix
 [Ur,~,~] = svd(Cr);  
evecNoise = Ur;

if fspec==1
    for gg=1:20
        [G,~,f] = quickmtm(evecNoise(:,gg),3,0,0,1);
        Noiseper(gg,1) = 1/f((G==max(G)));
        % disp(['TEOF (SSA) ',num2str(gg), ': Variance = ', num2str(varexp(gg)*100), ' Period = ', num2str(Dper(gg,1))])
        clear G
    end
        noisePeriods(:,i) = Noiseper;
        clear Noiseper
end



 lambda_red(i,:) = diag(U'*Cr*U)';
 lambda_red_norm(i,:) = diag(U'*Cr*U)'/trace(U'*Cr*U)';
end

%slambda_red = (lambda_red);

slambda_red =fliplr(flipud(sort(sort(lambda_red))));
slambda_red_norm =fliplr(flipud(sort(sort(lambda_red_norm)')));

% get some critical values for plotting and returning
mctest.c95_val_w=slambda_white(:,5);
mctest.c99_val_ar=slambda_red(:,2);
mctest.c95_val_ar=slambda_red(:,5);
mctest.c90_val_ar=slambda_red(:,10);
mctest.c50_val_ar=slambda_red(:,50);

mctest.c95_norm_w=slambda_white(:,5);
mctest.c99_norm_ar=slambda_red_norm(:,2);
mctest.c95_norm_ar=slambda_red_norm(:,5);
mctest.c90_norm_ar=slambda_red_norm(:,10);
mctest.c50_norm_ar=slambda_red_norm(:,50);

mctest.ar = lambda_red;

else

mctest = [];

end

%% Calculate the primary periodicity of the Reconstructed Components (eigenvectors)
warning off;
% if fspec==1
%     for gg=1:20
%         [G,~,f] = quickmtm(evec(:,gg),3,0,0,1);
%         Dper(gg,1) = 1/f((G==max(G)));
%         disp(['TEOF (SSA) ',num2str(gg), ': Variance = ', num2str(varexp(gg)*100), ' Period = ', num2str(Dper(gg,1))])
%         clear G
%     end
% end
if fspec==1
    for gg=1:20
        [G,~,f] = quickmtm(RC(gg,:),3,0,0,1);
        Dper(gg,1) = 1/f((G==max(G)));
        disp(['TEOF (SSA) ',num2str(gg), ': Variance = ', num2str(varexp(gg)*100), ' Period = ', num2str(Dper(gg,1))])
        clear G
    end
    orderedPeriods = Dper;
end


warning on;

%% Summary Figure, showing significance levels
% if graphics
if graphics
figure
orient tall
subplot(1,2,1)
semilogy([1:length(eval)],eval,'bd');
hold on;
if mct > 0
z = semilogy([1:length(eval)],mctest.c95_val_w,'b--',[1:length(eval)],mctest.c99_val_ar, 'r:',[1:length(eval)],mctest.c95_val_ar,'r-.',[1:length(eval)],mctest.c90_val_ar,'r--',[1:length(eval)],mctest.c50_val_ar,'r');
% legend(z,'95% confidence level (Guassian)','99% confidence level (AR1)','95% confidence level (AR1)','90% confidence level (AR1)','Mean Red Noise (50% confidence level)',4);
end
ylabel('Eigenvalue Amplitudes');
xlabel('Eigenvalue Rank');
hold off;

subplot(1,2,2);
semilogy([1:length(eval)],varexp,'bd');
hold on;
if mct>0
z = semilogy([1:length(eval)],mctest.c99_norm_ar, 'r:',[1:length(eval)],mctest.c95_norm_ar,'r-.',[1:length(eval)],mctest.c90_norm_ar,'r--',[1:length(eval)],mctest.c50_norm_ar,'r');
% legend(z,'99% confidence level (AR1)','95% confidence level (AR1)','90% confidence level (AR1)','Mean Red Noise (50% confidence level)',4);
end
ylabel('Normalized Eigenvalue Amplitudes');
xlabel('Eigenvalue Rank');
hold off;

% Figure showing Eigenvector waveform
figure
for ff=1:4;
    subplot(2,2,ff);
    plot([1:M],evec(:,ff));
    grid;
    title(['EOF ', int2str(ff)]);
    yl = ylim; ypos = ((yl(1,2) - yl(1,1))/20);
    if fspec==1;text(2,(yl(1,1)+ypos),['Dominant Period = ', num2str(Dper(ff,1),3), ' ',iunit]);end
end

clear ff

% Figure showing reconstructed components
figure
for ff=1:4;
    subplot(2,2,ff);
    plot(itime,RC(ff,:));
    grid;
    title(['Reconstructed Component ', int2str(ff)]);
    yl = ylim; ypos = ((yl(1,2) - yl(1,1))/20);
    xl = xlim; xpos = ((xl(1,2) - xl(1,1))/50);
    if fspec==1;text((xl(1,1)+ypos),(yl(1,1)+ypos),['Dominant Period = ', num2str(Dper(ff,1),3), ' ',iunit]);end
end


% tile the figure windows, if tilefigs.m exists -- available from Matlab Central: 
 if exist('tilefigs') >= 1; tilefigs; end
end