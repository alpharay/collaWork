%%Author: KB
%Purpose: Trying to reproduce the results of the reference paper below
%1. Adapted from ...
% Date:  2020. 02. 10. / 12:22:19 KST

% Version variation:  
% 1.

% Reference paper: 
%[1] Cooperative jamming and Power Allocation for Wireless
% Relay Networks in Presence of Eavesdropper

clear all;
close all;
clc;


%% SETTINGS and PREALLOCATIONS
tx=1;%number of transmitters
rx=1;%number of receiver_poss
mu = 0; % mean of Guassian distribution


%--settings based on ... ----
iter=10^3; % numbeer of Monte Carlo simulations
%d=30;%receiver's distance; destination fixed at 30m
d=1;%receiver's distance; destination fixed at 30m % (2020/02/10)

%eavesdynPos=(-16:2:32);%eavesdroppers dynamnic position; necessary to plot the graph of secrecy gap at with eavesdropper at different distances from the sources
eavesdynPos=(-4:0.05:3);%eavesdroppers dynamnic position; necessary to plot the graph of secrecy gap at with eavesdropper at different distances from the sources
R_sec_vec = zeros(size(eavesdynPos));%(2020/02/10)

tau=4; %path loss exponent
phi=0; %phase offset from source to destination
pathloss_comp=-tau/2;%pathloss component
h= d.^(pathloss_comp)*exp(1j*phi); %line of sight model
%------------------------------------------------------------------------



%% POWER in dB

%------------------------------------------------------------------------
%Total Power budget =20dBm,
%(a) Phase 1(Source to destination; active jammers)=>P_S+P_J1
%(b)Phase 2(Relays to destination; jammers active again) =>P_R+P_J2.
% Therefore P=P_S+P_J1+P_J2+P_R
%------------------------------------------------------------------------

num_relays=1; %number of relays
num_jammers=1;%number of jammers
l=num_jammers;%maximum number of jammers
n=num_relays;%maximum number of relays

% When "The average signal-to-noise ratio (SNR) of the S --> R and
% R --> D links, i.e., gamma_SR_bar and gamma_RD_bar are fixed at
% 15 dB"
gamma_SR_bar = 10^(15/10);
gamma_RD_bar = gamma_SR_bar;

% ---Noise specific settings ----
sigmaSqrd = 1;
%---------------------------------


p_budget_Dbm=20; %total power budget in dBm
p_budget_inWatts=10^((p_budget_Dbm-30)/10);%total power budget in dBm
p_budget_inWatts_divided=p_budget_inWatts/(1+2*(num_jammers)+num_relays);%dividing up the power 
P_S=p_budget_inWatts_divided;%source power in watts
P_R=(p_budget_inWatts_divided)*num_relays;%source power in watts
P_J1=(p_budget_inWatts_divided)*num_jammers;%power in watts for jammers in 1st phase
P_J2=(p_budget_inWatts_divided)*num_jammers;%power in watts for jammers in 2nd phase


% P_S=15;%source power in watts
% P_R=2;%source power in watts
% P_J1=2;%power in watts for jammers in 1st phase
% P_J2=1;%power in watts for jammers in 2nd phase


P_J=P_J1+P_J2;% total jammer power


p=[P_S P_R P_J1 P_J2]'; %matrix of all the powers

% num_relays=input('Please enter the number of relays: ');%number of relays
% num_jammers=input('Please enter the number of relays: ');%number of jammers

%-------------------------------------------------------------------------------------------------------------------
%% DISTANCES

snr_gap=zeros(1,length(eavesdynPos));

snr_gap_cummulator=0; %for cummulating the sum of snr gap to be finally averaged out after  the iterations

%--source position--
%source_pos = zeros(1,1);
source_pos = [-1;0];

for eavesd_dist_iter=1:size(eavesdynPos,2),
%for iter_num=1:iter,
    %% SOURCE
    % USING DYNAMIC RELAY/JAMMER POSITIONS
    
    %% RELAY(S)
    %A1. NON-DYNAMIC RELAY POSITION    
    x_rel_dist =  0;
    y_rel_dist =  0; 
    angles_SR = atan((y_rel_dist-source_pos(2))/(x_rel_dist-source_pos(1)));
    
    %A2. DYNAMIC II (MOST DYNAMIC).
    disk_ctr=[0,0]; %source as center
    disk_rad=18;%radius of disk used to distribute relay positions

%     [x,y,angles_out]=randomNodeCoordinateGenerator(disk_ctr,disk_rad,num_relays,0,2*pi);
%     y_rel_dist=y';
%     x_rel_dist=x';    
%     angles_SR=angles_out';%angles from source to relay. Giving their phase offset component
    
    figure(1);
    plot(source_pos(1),source_pos(2),'rd');% for seeing the positioning of the source
    hold on;
    plot(d,0,'m^');% for seeing the positioning of the destination
    plot(x_rel_dist,y_rel_dist,'g*');% relay distances
    
    
    
     %% JAMMER(S)
%     %A1. NON-DYNAMIC RELAY POSITION
%     y_jam_dist = 0;
%     x_jam_dist = 0
%     angles_SJ = atan(y_jam_dist./x_jam_dist);
    
    %B2. DYNAMIC II (MOST DYNAMIC).
    [x_jammers,y_jammers,angles_jammers_out]=randomNodeCoordinateGenerator(disk_ctr,disk_rad,num_jammers,0,2*pi);
    y_jam_dist=y_jammers';%y-coordinates of jammers with source as center
    x_jam_dist=x_jammers';%y-coordinates of jammers with source as center
    angles_SJ=angles_jammers_out';%angles from source to jammers. Giving their phase offset component
    
    plot(x_jam_dist,y_jam_dist,'bo');% for seeing the positioning of the jammers
    %------------------------------------------------------------------------------------------
    
        
    R_sec_term3_summer = 0;%secrecy capacity agreggator (2020/02/10)
    %for eavesd_dist_iter=1:size(eavesdynPos,2),%(2020/02/10)
    gamma_SE_summer = 0; %(2020/02/11)
    gamma_RE_summer = 0; %(2020/02/11)
    
    for iter_num=1:iter,%(2020/02/10)
        
        %         try
        %% DISTANCES
            
        relay_pos=[x_rel_dist;y_rel_dist];%PNB: relay positions
        jammer_pos=[x_jam_dist;y_jam_dist];%PNB: jammer positions
        eavesd_pos=[eavesdynPos(eavesd_dist_iter);0];% eavesdroppers position
        plot(eavesd_pos(1),eavesd_pos(2),'rs');% for seeing the positioning of the eavesdropper
        dest_pos=[d;0];% destination position
        
        %---source originated distances---        
        
        source2relay_vec = relay_pos - repmat(source_pos,1,size(relay_pos,2));
        source2relay=sqrt(sum((source2relay_vec).^2)).^pathloss_comp;%PBN: source(zero origin) to relay.
        %receiver_pos=[30*ones(1,num_relays);zeros(1,num_relays)];%PNB: destination(receiver_pos's) position from origin
        
        source2dest_vec=dest_pos-source_pos;%vector from source to destination.
        source2dest=sqrt(sum((source2dest_vec).^2)).^pathloss_comp;%PNB: source to destination distance.
        
        source2eavesd_vec=eavesd_pos-source_pos;%vector from source to eavesdropper.
        source2eavesd=sqrt(sum((source2eavesd_vec).^2)).^pathloss_comp;%PNB: source to eavesdropper distance.
        
        %---relay originated distances---
        relay2dest_vect=repmat(dest_pos,1,size(relay_pos,2))-relay_pos;%matrix from relays to destination.
        relay2dest=sqrt(sum((relay2dest_vect).^2)).^pathloss_comp;%PNB: relays to destination distance.
        %% finding phase angles from relay to destinations
        xS1=relay2dest_vect(1,:);
        yS1=relay2dest_vect(2,:);
        phi_relay2dest=atan(yS1./xS1);%(angles in radians)
        %phi_relay2dest=rad2deg(atan(yS1./xS1));%(angles in degrees)
        
        relay2eavesd_vect=repmat(eavesd_pos,1,size(relay_pos,2))-relay_pos;%matrix from relays to eavesdropper.
        relay2eavesd=sqrt(sum((relay2eavesd_vect).^2)).^pathloss_comp;%PNB: relays to eavesdropper distances.
        %% finding phase angles from relay to destinations
        xS2=relay2eavesd_vect(1,:);
        yS2=relay2eavesd_vect(2,:);
        phi_relay2eavesd=atan(yS2./xS2);%(angles in radians)
        %phi_relay2eavesd=rad2deg(atan(yS2./xS2));%(angles in degrees)
        
        
        %---jammer originated distances and angles---
        
        if num_jammers || num_relays > 1,
            for i=1:num_relays,
                for j=1:num_jammers,
                    jam2relay_vect(:,j)=relay_pos(:,i)-jammer_pos(:,j)
                    jam2relay_vect2dist(j)=sqrt(sum(jam2relay_vect(:,j).^2)).^pathloss_comp%PNB: distance of jammers to specific relay
                end
                jam2relay_vec2distvec(:,i)=jam2relay_vect2dist;%finished computing distances from jammers to specific relay added to  final distance matrix
                %% finding the phase angles from jammers to relays
                x_cord=jam2relay_vect(1,:);
                y_cord=jam2relay_vect(2,:);
                phi_relay2jammer=atan(y_cord./x_cord);% angles from jammers to specific relay calculated (in radians)
                %phi_relay2jammer=rad2deg(atan(y_cord./x_cord));%angles from jammers to specific relay calculated (in degrees)
                phi_relay2jammer_complete(:,i)=phi_relay2jammer;% angles from jammers to specific relay added to final angle matrix
            end
            jam2relay=jam2relay_vec2distvec';
            jam2relay_angles=phi_relay2jammer_complete';
        else disp('Error in specifing the relay size');
            
        end
        
        
        jam2dest_vect=repmat(dest_pos,1,size(jammer_pos,2))-jammer_pos;%matrix from jammers to destination.
        jam2dest=sqrt(sum((jam2dest_vect).^2)).^pathloss_comp;%PNB: jammers to destination distances.
        %% finding phase angles from jammer to destinations
        xS3=jam2dest_vect(1,:);
        yS3=jam2dest_vect(2,:);
        phi_jam2dest=atan(yS3./xS3);%angle (in radians)
        %phi_jam2dest=rad2deg(atan(yS3./xS3));%angle (in degrees)
        
        jam2eavesd_vect=repmat(eavesd_pos,1,size(jammer_pos,2))-jammer_pos;%matrix from jammers to eavesdropper.
        jam2eavesd=sqrt(sum((jam2eavesd_vect).^2)).^pathloss_comp;%PNB: jammers to eavesdropper distances.
        xS4=jam2eavesd_vect(1,:);
        yS4=jam2eavesd_vect(2,:);
        phi_jam2eavesd=atan(yS4./xS4);%angle (in radians)
        %phi_jam2eavesd=rad2deg(atan(yS4./xS4));%angle (in degrees)
        
        
        %% CHANNELS
               
        %---CHANNEL ORIGIN: FROM SOURCE-----
        %source to relay
        phi_source2relay=angles_SR;
        h_SR1=exp(1i*phi_source2relay);
        h_SR=source2relay.*h_SR1; % Source to Relay channel
        
        %source to destination/receiver
        phi_source2dest=phi;
        h_SD1=exp(1i*phi_source2dest);
        h_SD=source2dest.*h_SD1;
        
        %source to eavesdropper
        phi_source2eavesd=atan(source2eavesd);%angle (in radians)
        h_SE1=exp(1i*phi_source2eavesd);
        h_SE=source2eavesd.*h_SE1;
        
        
        %---CHANNEL ORIGIN: RELAY-----
        %relay to destination
        h_RD1=exp(1i*phi_relay2dest);
        h_RD=relay2dest.*h_RD1; % Source to Relay channel
        h_RD=h_RD';
        
        %relay to eavesdropper
        %h_RE1=(randn(1,num_relays)+randn(1,num_relays)*1i)/sqrt(2);
        h_RE1=exp(1i*phi_relay2eavesd);
        h_RE=relay2eavesd.*h_RE1; % Source to Relay channel
        h_RE=h_RE';
        
        %---CHANNEL ORIGIN: JAMMER-----
        %jammer to relay
        
        h_JR1=exp(1i*jam2relay_angles);
        h_JR_vec=jam2relay.*h_JR1;
        h_JR=h_JR_vec;
        h_JR=h_JR';
        
        %jammer to destination/receiver
        h_JD1=exp(1i*phi_jam2dest);
        h_JD=jam2dest.*h_JD1;
        h_JD=h_JD';
        
        %jammer to eavesdropper
        h_JE1=exp(1i*phi_jam2eavesd);
        h_JE=jam2eavesd.*h_JE1;
        h_JE=h_JE';
        
       
        
        %% SECRECY CAPACITY CALCULATION
        
        
        %gamma_SR = gamma_SR_bar;% with fixed gammas        
        %gamma_RD = gamma_RD_bar;% with fixed gammas
        gamma_SR = P_S*abs(h_SR).^2./sigmaSqrd;% When gammas not fixed
        gamma_RD = P_R*abs(h_RD).^2./sigmaSqrd;% When gammas not fixed
        
        gamma_SE = P_S*abs(h_SE).^2./sigmaSqrd; 
        gamma_SE_summer = gamma_SE_summer + gamma_SE;
        
        gamma_RE = P_R*abs(h_RE).^2./sigmaSqrd;
        gamma_RE_summer = gamma_RE_summer + gamma_RE;
        
        A = 3*(gamma_RD).^2*gamma_RE*gamma_SE;
        
        B_term1 = -2*gamma_RD*gamma_RE*gamma_SR*(1+gamma_SE);
        B_term2 = -2*gamma_RD*gamma_SE*(gamma_RD*(1+gamma_RE)-gamma_RE);
        B = B_term1 + B_term2;
        
        C_term1 = gamma_SR*(1+gamma_SE)*(gamma_RD*(1+gamma_RE)-gamma_RE);
        C_term2 = -gamma_RD*gamma_SE*(1+gamma_RE);
        C = C_term1 + C_term2;
        
        betaHat_num = -B-sqrt(B^2-4*A*C);
        betaHat_denom = 2*A;
        betaHat = betaHat_num/betaHat_denom;
        
        
        
        
        R_sec_term3 = log((1+gamma_SE)*(1 + gamma_RE));        
        
        R_sec_term3_summer = R_sec_term3_summer + R_sec_term3;
    end
    gamma_SE_Bar = gamma_SE_summer/iter;
    gamma_RE_Bar = gamma_RE_summer/iter;
    R_sec_term3_avg = R_sec_term3_summer/iter;
    
    
    if betaHat < (gamma_SR/gamma_RD),
        alphaHat = gamma_RD*betaHat/gamma_SR;
    elseif betaHat > (gamma_SR/gamma_RD),
        alphaHat = 1;
    end
    
    alpha_ = alphaHat;
    beta_ = betaHat;    
    
    
    R_sec_term1 = 0.5*log(1+(1-alpha_)*gamma_SE_Bar);
    R_sec_term2 = 0.5*log((1+beta_ * gamma_RD)*(1+(1-beta_)*gamma_RE_Bar));    
    R_sec = R_sec_term1 + R_sec_term2 - 0.5*R_sec_term3_avg;
    
    if R_sec < 0,
        R_sec = 0;
        %R_sec = R_sec*-1;
    end
    
    R_sec_vec(eavesd_dist_iter) = R_sec;%secrecy capacity vector
end



figure(2)
plot(eavesdynPos,R_sec_vec,'-o');%plot of SNR gap in dB
%semilogy(eavesdynPos,R_sec_vec,'-o');%plot of SNR gap in dB

grid on;
ylabel('Secrecy capacity[dB]');
xlabel('distance from source [m]');
runtimeTimeStamp();
