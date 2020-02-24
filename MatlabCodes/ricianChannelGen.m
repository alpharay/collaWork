% Purpose: Matlab code for Generating the small scale fading model for
% Rician channel
% Author:KB
% Date: 2020. 02. 18. / 11:51:15 KST
% References:Simulation of Digital Communications Systems using Matlab
%% by Mathuranathan Viswanathan
function [h_rician] = ricianChannelGen(EbN0dBvals_in,dataGenOn_in)
dataGenOn = dataGenOn_in; %set this to one when you want to get a plot for BPSK or see the capacity plots

%EbN0dB=-5:2:20; %Eb/N0 in dB overwhich the performance has to be simulated
EbN0dB = EbN0dBvals_in; %Eb/N0 in dB overwhich the performance has to be simulated

if dataGenOn == 1
    %% 1. GENERATING RICIAN CHANNEL (SMALL SCALE FADING COMPONENT - without pathloss component multiplied) with data transmission (i.e N > 1) considered    
    
    %Eb/N0 Vs BER for BPSK over Rician Fading Channel with AWGN noiseclc;clear all;
    %--------Inputs-------------------------------------------------
    N=10^6; %Number of data samples to send across the Rician Channel
    totPower=1; %Total power of LOS path & scattered paths
    K=[1 2 5 10 20 30]; %A list of Ricial K factors to simulate
    %--------------------------------------------------------------
    d=rand(1,N)>0.5; %data generation
    x=2*d-1; %BPSK modulation
    simBER_ricean=zeros(1,length(EbN0dB));
    plotStyle={'b*-','r*-','k*-','g*-','m*-','c*-'};
    
    for index =1:length(K)
        k=K(index);
        %Derive non-centrality parameter and sigma for the underlying
        %Gaussian RVs to generate the Rician Envelope
        s=sqrt(k/(k+1)*totPower); %Non-Centrality Parameter
        sigma=totPower/sqrt(2*(k+1));
        
        for i=1:length(EbN0dB)
            noise=1/sqrt(2)*(randn(1,N)+1i*randn(1,N)); %AWGN noise with mean=0 var=1
            h=((sigma*randn(1,N)+s)+1i*(randn(1,N)*sigma+0)); %Rician Fading - single tap
            n = noise*10^(-EbN0dB(i)/20); %Scaling the noise for required Eb/N0
            y_ricean=h.*x+n; %received signal through Rician channel
            %Coherent Receiver for Rician Channel
            y_ricean_cap=y_ricean./h; %Assuming that h is known at the signal accurately
            r_ricean=real(y_ricean_cap)>0; %received symbols = 1 is real part > 0 or else it is 0
            %Receiver for AWGN channel
            simBER_ricean(i)=sum(xor(d,r_ricean));
        end
        
        simBER_ricean=simBER_ricean/N;
        %Simulated BER;
        EbN0=10.^(EbN0dB/10); %Eb/N0 in Linear Scale
        semilogy(EbN0dB,simBER_ricean,plotStyle{index},'LineWidth',2);hold on
        legendInfo{index} = ['K = ' num2str(K(index))];
    end
    
else
    %% 2. JUST GENERATING RICIAN CHANNEL (SMALL SCALE FADING COMPONENT - without pathloss component multiplied) without data transmission (i.e N = 1)   

    %Eb/N0 Vs BER for BPSK over Rician Fading Channel with AWGN noiseclc;clear all;
    %--------Inputs-------------------------------------------------
    N=1; %Number of data samples to send across the Rician Channel
    totPower=1; %Total power of LOS path & scattered paths
    %K=[1 2 5 10 20 30]; %A list of Ricial K factors to simulate
    K=5; %A list of Ricial K factors to simulate
    %--------------------------------------------------------------
    
    for index =1:length(K)
        k=K(index);
        %Derive non-centrality parameter and sigma for the underlying
        %Gaussian RVs to generate the Rician Envelope
        s=sqrt(k/(k+1)*totPower); %Non-Centrality Parameter
        sigma=totPower/sqrt(2*(k+1));
        
        for i=1:length(EbN0dB)
            noise=1/sqrt(2)*(randn(1,N)+1i*randn(1,N)); %AWGN noise with mean=0 var=1
            h=((sigma*randn(1,N)+s)+1i*(randn(1,N)*sigma+0)); %Rician Fading - single tap
            n = noise*10^(-EbN0dB(i)/20); %Scaling the noise for required Eb/N0
            y_ricean=h+n; %received signal through Rician channel
            h_rician(i) = y_ricean;
        end
        
    end
    
end