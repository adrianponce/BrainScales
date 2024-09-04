function [FC,Cov,Ct,lags,PSD,freqs,Csp,Remax]=LangevinSystem(C,tau,sigma,varargin)
%
% linear stochastic network of N nodes
% each node evolves as follows:
%
% dx/dt = -x/tau + C*x + noise
%
% where C are the network interactions.
% 
% Calculates the FC, the lagged-covariance,the power spectral density, 
% and the cross-spectrum analytically by means of the Jacobian matrix
%
% Inputs:
%   - C   : connectivity matrix (N-by-N)
%   - tau : characteristic time
%   - sigma : noise amplitude (scalar)
%   - varargin (optiona) : lags of the lagged-covariance (default : 0:10 seconds)
%
% Outputs:
%   - FC  : correlation matrix of real(z)
%   - Cov : covariance matrix of real(z)
%   - Ct  : lagged covariance of real(z)
%   - PSD : power spectral density of real(z)
%   - lags: lags used for Ct
%   - freqs: frequencies used for PSD
%   - Csp : cross-spectrum
%   - remax : max. real part eigenvalue
%
% Adrián Ponce-Alvarez 06-07-2022
%--------------------------------------------------------------------------


N = size(C,1);

% Jacobian matrix:
A = (C - eye(N)*tau);
    
% input noise covariance:
Qn = sigma^2*eye(N);
    
% Check stability:
[~,d] = eig(A);
d = diag(d);
Remax = max(real(d));
if Remax >= 0
   disp('Warning: not stable') 
   FC = []; Cov=[]; Ct=[]; lags=[]; PSD=[]; freqs=[]; Csp=[];
   return
end
    
% Covariance equation:
Cov = lyap(A,Qn); % Solves the Lyapunov equation: A*Cv + Cv*A' + Qn = 0
FC=corrcov(Cov);

% Lagged covariance:
%--------------------------------------------------------------------------
if nargin < 6
L=10;
lags=0:.1:L; % default lags to evaluate the lagged-covariances
else
lags = varargin{1}; % if asked, use these lags    
end
num_lags = length(lags);
Ct=zeros(N,N,num_lags);

 for n=1:length(lags)
     t = lags(n);
     Y=expm(A*t)*Cov;
     Ct(:,:,n) = Y;
 end


% Power spectrum:
% wo needs to be in radians per second!
%--------------------------------------------------------------------------
if nargin < 7
dfr=0.005;
freqs=0.05:dfr:6; % default frequencies to evaluate the PSDs and cross-spectrum
else
freqs = varargin{2}; % if asked, use these frequencies   
end
numF=length(freqs);

S   = zeros(numF,N);
Csp = zeros(N,N,numF,'single');

     for k=1:numF
            f   = freqs(k);
            J   = A+1i*f*(2*pi)*eye(N);
            Q   = inv(J);
            X   = Q*(Q'); 
            Y   = X(1:N,1:N);
            Csp(:,:,k) = single(Y);  %cross-spectrum    
            S(k,:) = diag(Y); % spectral density
            % this would be equivalent:
            %for i=1:N
            %    S(k,i)=sum(abs(Q(i,:)).^2);
            %end     
     end
            
% PSD:
%(the factor 2 is due to adding the contributions from positive and
%negative frequencies)
PSD = 2*sigma^2*S;

% cross-spectrum:
Csp = 2*sigma^2*Csp;
 
 return