function [rE,rI,time] = WCcoupled_StochSim(Wglobal,Io,tmax,sigma,varargin)
%
% Stochastic simulation of the large-scale Wilson-Cowan model
%
% Inputs:
% - Wglobal : coupling connectivity between brain regions
% - Io : background inputs
% - tmax : simulation time, multiples of tauE
% - varargin : transitory time till stationary regime (optional, default: tmax)
%
% Outputs:
% - rE  : E firing rates
% - rI  : I firing rates
% - time: time array
%
% A. Ponce-Alvarez, 13/08/2024
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

%tauEms = 0.010; % we assume tauE = 10 ms

% transfer functions:
%-------------------------------
scale = 1;
F = @(x) scale./(1 + exp(-x) );


% simulation parameters:
%--------------------------------------------------------------------------
dt=0.005;

tspan=0:dt:tmax;
L = length(tspan);
if nargin>4
   Ttran = varargin{1};
else
   Ttran = L;
end

% downsampling:
ds = 20;
Tds = length(0:ds*dt:tmax)-1;
time = 0:ds*dt:(tmax-ds*dt);
% % Noise amplitude
% sigma = 0.001;

% Simulation :
%--------------------------------------------------------------------------
r = rand(N,1);
% transient:
for t = 1:Ttran            
     u = W*r + Io;            
     K = feval(F,u);            
     r = r + dt*(-r+ K)./tau + sigma*sqrt(dt)*randn(N,1);            
end        
% store
R = zeros(Tds,N);
tt = 0;               
for t = 1:L            
     u = W*r + Io;            
     K = feval(F,u);            
     r = r + dt*(-r+ K)./tau + sigma*sqrt(dt)*randn(N,1);            
     if mod(t,ds)==0
        tt=tt+1;
        R(tt,:)=r;
     end   
end
 
 % steady-state: r and u.
        
rE = R(:,1:N/2);
rI = R(:,N/2+1:end);

