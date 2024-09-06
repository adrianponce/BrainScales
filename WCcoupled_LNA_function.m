function [Cov,CeE,rE,rI,pw,peakfreq,freqs,Ct,lags] = WCcoupled_LNA_function(Wglobal,Io,freqs,lags)

%
% Statistics of the network of the Wilson-Cowan large-scale model using linear noise
% approximation
%
% Inputs:
% - Wglobal : coupling connectivity between brain regions
% - Io : background inputs
%
% Outputs:
% - CeE : correlation matrix E-E
% - rE  : E firing rates
% - rI  : I firing rates
% - peakfreq : E peak frequencies
% - pw : E power spectral densities
% - freqs : tested frequencies
% - Ct : lagged-covariances
% - lags : time lags for Ct
%
% A. Ponce-Alvarez, 14/08/2024
%--------------------------------------------------------------------------

% model parameters:
%--------------------------------------------------------------------------
% number of neural populations:
N = 2*size(Wglobal,1);

% Local connectivity:
wII=4;
wIE=16;
wEI=12;
wEE=12;

% Full Connectivity matrix:
W = [Wglobal+wEE*eye(N/2), -wEI*eye(N/2);...
     Wglobal+wIE*eye(N/2), -wII*eye(N/2)];


% Time constants:
tauE = 1; % (non-dimensional)
tauI = 2; 
tau = [tauE*ones(N/2,1);tauI*ones(N/2,1)];

tauEms = 0.010; % we assume tauE = 10 ms

% transfer functions:
%-------------------------------
scale = 1;
F = @(x) scale./(1 + exp(-x) );
% Derivative:
DF = @(x) scale*feval(F,x)*(1-feval(F,x));


% simulation parameters:
%--------------------------------------------------------------------------
dt=0.005;
tmax = 500; % multiples of tauE
tspan=0:dt:tmax;
L = length(tspan);

% downsampling:
ds = 5;
Tds = length(0:ds*dt:tmax)-1;

% Noise amplitude
sigma = 0.001;
Qn = sigma^2*eye(N);

% frequencies for PSD calculation:
%freqs = 0.001:.5:50;
numF = length(freqs);


% Simulation to get fixed points:
%--------------------------------------------------------------------------
r = rand(N,1);
% transient:
for t = 1:L            
     u = W*r + Io;            
     K = feval(F,u);            
     r = r + dt*(-r+ K)./tau;            
end        
% store
R = zeros(Tds,N);
tt = 0;               
for t = 1:L            
     u = W*r + Io;            
     K = feval(F,u);            
     r = r + dt*(-r+ K)./tau;            
     if mod(t,ds)==0
        tt=tt+1;
        R(tt,:)=r;
     end   
end
 
 % steady-state: r and u.
        
RstatE = R;
maxR = max(RstatE);
minR = min(RstatE);
delta = maxR-minR;         
         
        if all(delta)<10e-10 % if no self-sustained oscillations 
            
        % Linear noise approximation (LNA)    
        % Jacobian matrix:
        f1 = scale*r.*(1-r);
        Jmat = W.*repmat( f1./tau ,[1 N]) - diag(1./tau);  
        
      
        % Covariance matrix (LNA)  
        [V,D] = eig(Jmat);
        d=diag(D);
        U=V';
        Q=V\Qn/U;
        M = repmat(d,1,N);
        Y = M + M';
        P = -Q./Y;
        rho=V*P*V';
        Cov=real(rho);
            
        ss = repmat(diag(Cov),[1,N]).*repmat(diag(Cov),[1,N])';    
        % Correlation matrix:
        Corr= Cov ./ sqrt(ss);
            
        % store E correlations
        CeE = Corr(1:N/2,1:N/2);
        % store firing rates:
        rE = r(1:N/2);
        rI = r(N/2+1:end);
      
                                
         % Power Spectrum
         %-----------------------------------------------------------------              
         S   = zeros(numF,N);
         for kk=1:numF
             ff=freqs(kk);
             Amat      = Jmat/(tauEms)+1i*ff*(2*pi)*eye(N);
             AmatInv   = inv(Amat);
             X=AmatInv*AmatInv'; %cross-spectrum
             S(kk,:) = diag(X);
         end

        pw = 2*sigma^2*S;

        peakfreq = nan(N,1);
        % peak frequencies:
        for nn=1:N       
        [~,loc]=findpeaks(pw(:,nn),freqs);
            if ~isempty(loc)
            peakfreq(nn) = loc;
            else
            peakfreq(nn) = 0;
            end
        end
        
        % Autocovariance:
        num_lags = length(lags);
        Ct=zeros(N,N,num_lags);

         for n=1:length(lags)
             t = lags(n);
             Y=expm(Jmat*t)*Cov;
             Ct(:,:,n) = Y(1:N,1:N);
         end
       
       
        
        else % if self-sustained oscillations:
                  
        Cov = [];    
        CeE = [];
        rE = [];
        rI = [];
        pw = [];
        peakfreq = [];
        Ct = [];
                                 
        end
        
        
    
 
 
 

