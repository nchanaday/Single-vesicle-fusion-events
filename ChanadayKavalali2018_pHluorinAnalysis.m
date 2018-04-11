% ----------------------------------------------------------------------- %
% Fluorescence analysis: single-vesicle peaks and multi-vesicular peaks   %
% Denoising the signal with a combination of FFT and non-linear filtering %
% ----------------------------------------------------------------------- %

% --> For this script to work, you need to save the data properly organized in an Excel file.

% --> Do you want to perform bleaching and backgroud correction?? 0=No; 1=Yes.
BBC=1;
% --> Do you want to perform FFT filtering?? 0=No; 1=Yes.
FFTF=0;
% --> Do you want to perform Chung-Kennedy filtering?? 0=No; 1=Yes.
CKF=1;
% --> For CKF, how many rounds of denoising do you wanna perform??
dnRounds=2;  % use 2 or 3;

for nfile=1:2
    if nfile==1
        in_filename='Data1.xlsx';
        out='Data1_';
    elseif nfile==2
        in_filename='Data2.xlsx';
        out='Data2_';
    end

[status,sheets]=xlsfinfo(in_filename);
NSheets=numel(sheets);

for sheet=1:NSheets
    out_filename=[out,num2str(sheet)];
    Data=xlsread(in_filename,sheet);
    
%% --- PARAMETERS --- %

% (units: TIME in seconds and FLUORESCENCE in a.u.)

    tini=Data(1,2);                                                         % time of first single AP stimulation
    Nstim=Data(2,2);                                                        % number of stimulus
    tinterval=Data(3,2);                                                    % interval between stimulus (inverse of frequency)
    tpeak1=Data(4,2);                                                       % time of the first peak of strong stimualtion (e.g. if 20 Hz 100 AP was applied at 100 s, then the peak's time is 105 s)
    tpeak2=Data(5,2);                                                       % sencond "big" peak
    tpeak3=Data(6,2);                                                       % third "big" peak
    tpeak4=Data(7,2);                                                       % last strong stimulation (i.e. high K of NH4)
    Total_points=Data(8,2);
    Amplitude_limit=Data(9,2);
    increment=Data(10,2);
    min_fit_points=Data(11,2);
    t=Data(:,3);
    F=Data(:,6:end);
    
% --- % 

tfin=tpeak1-tinterval;      
iini=find(t>tini,1);
time=(tini-0.5);
ifin=find(t>tfin,1);
ipeak1base=find(t>(tpeak1-5),1);
ipeak1=find(t>tpeak1,1);
ipeak2=find(t>tpeak2,1);
ipeak3=find(t>tpeak3,1);
ipeak4=find(t>tpeak4,1);


%% --- BLEACHING and BACKGROUND correction --- %
% Based on: Vicente et al. (2007) DOI: 10.1088/1742-6596/90/1/012068.

if BBC==1

expdecay2=@(y,range)y(1)*exp(-y(2)*range)+y(3);
y0=[0,0,0];
options=optimoptions(@lsqcurvefit,'MaxFunEvals',1000,'MaxIter',1000,'Algorithm','levenberg-marquardt');
lb=[];
ub=[];

for g=1:size(F,2)
    range=t(10:ifin);
    expparam=lsqcurvefit(expdecay2,y0,range,F(10:ifin,g),lb,ub,options);
    FunF=expparam(1)*exp(-expparam(2)*t)+expparam(3);
    
    for j=1:size(t)
        ratio(j)=FunF(j)/FunF(1);                                           % decay due to bleaching 
        ROI(j,g)=(F(j,g)/ratio(j))-(mean(FunF(1:iini-10)));                 % first term: bleaching correction; second term: background correction 
    end
end

elseif BBC==0
    ROI=F;
end


%% --- FFT FILTERING --- %
% Fast (Discrete) Fourier Transform of the traces with posterior
% elimination of high frequency components of the signal:

if FFTF==1

tFFT=t(1:ipeak1base);  % to avoid strong stimulation part of the signal (it generates artifacts).
xFFT=ROI(1:ipeak1base,:);
stimfreq=1/tinterval;  % to not consider the stimulation frequency (and its multiples) for the frequency distribution analysis.
Nquist=fix(size(tFFT,1)/2);
k=(0:Nquist)';
freq=k/tFFT(end);  % vector of frequencies (it ranges from 0 to the imaging speed).
y=fft(xFFT);
absy=abs(y(1:Nquist+1,:));
power=(2/tFFT(end))*(absy.^2);  % power spectra density calculation.
min=find(freq>stimfreq,1);
for i=1:size(power,2)
[pks(:,i),locs(:,i)] = findpeaks(power(min:end,i),freq(min:end,1),'SortStr','descend','NPeaks',50);
% pks = maximum power values (50 highest).
% locs = frequencies of the pks --> more commom frequencies in the sample.
end
cutoff=freq(end)/5;  % cutoff frequency for the cleaning of the trace.
clean_freq=find(freq<=cutoff);
clean_y=zeros(size(y));
for i=1:size(y,2)
clean_y(clean_freq,i)=y(clean_freq,i);
end
ROI(1:ipeak1base,:)=abs(ifft(clean_y));  % recovery of the trace, cleaned (ifft=inverse FFT).

end


%% --- DENOISING OF THE FLUORESCENCE SIGNAL --- %
% Based on: Chung and Kennedy (1991); Reuel et al. (2001).
% It's a filtering on the time domain based on signal amplitudes.

if CKF==1

% ^^^ NOTE: User can change following parameter to improve fit ^^^

for round=1:dnRounds
    
N = [2 4 8 16 32 64];                                                      % Sampling window length (size of the forward and backward mean for the predictors)
K = length(N);                                                             % Number of forward and backward predictors
M = 20;                                                                    % Analysis window to compare the predictors
P = 10;                                                                    % Weighting factor

% --- Forward-Backward non-linear Algorithm

    for y=1:size(ROI,2)
        testROI(:,1)=ROI(:,y);
        ltime=length(t);
        I_avg_f = zeros(ltime,K);
        I_avg_b = zeros(ltime,K);
        
        for g = 1:ltime
            for k = 1:K
                % Average forward predictor
                window = N(k);
                if g == 1
                    I_avg_f(g,k) = testROI(1,1);
                elseif g - window - 1 < 0
                    I_avg_f(g,k) = sum(testROI(1:g-1,1))/g;
                else
                    epoint = g - window;
                    spoint = g - 1;
                    I_avg_f(g,k) = sum(testROI(epoint:spoint,1))/window;
                end
                % Average backward predictor
                if g == ltime
                    I_avg_b(g,k) = testROI(g,1);
                elseif g + window > ltime
                    sw = ltime - g;
                    I_avg_b(g,k) = sum(testROI(g+1:ltime,1))/sw;
                else
                    epoint = g + window;
                    spoint = g + 1;
                    I_avg_b(g,k) = sum(testROI(spoint:epoint,1))/window;
                end
            end
        end
        
        % Non-normalized forward and backward weights:
        f = zeros(ltime,K);
        b = zeros(ltime,K);
        for i = 1:ltime
            for k = 1:K
                Mstore_f = zeros(M,1);
                Mstore_b = zeros(M,1);
                Pi=1/(2*K); % Natali added this (see paper by Chung and Kennedy)
                for j = 0:M-1
                    t_f = i - j;
                    t_b = i + j;
                    if t_f < 1
                        Mstore_f(j+1,1) = (testROI(i,1) - I_avg_f(i,k))^2;
                    else
                        Mstore_f(j+1,1) = (testROI(t_f,1) - I_avg_f(t_f,k))^2;
                    end
                    % eqn. (4) and (5) in paper by Chung ang Kennedy:
                    if t_b > ltime
                        Mstore_b(j+1,1) = (testROI(i,1) - I_avg_b(i,k))^2;
                    else
                        Mstore_b(j+1,1) = (testROI(t_b,1) - I_avg_b(t_b,k))^2;
                    end
                end
                f(i,k) = Pi*(sum(Mstore_f)^(-P));
                b(i,k) = Pi*(sum(Mstore_b)^(-P));
            end
        end
        
        % Vector of normalization factors for the weights:
        C = zeros(ltime,1);
        for i = 1:ltime
            Kstore = zeros(K,1);
            for k = 1:K
                Kstore(k,1) = f(i,k) + b(i,k);
            end
            C(i,1) = 1/sum(Kstore);
        end
        
        % Putting parameters together and solving for intensities:
        ROIclean = zeros(ltime,1);
        for i = 1:ltime
            TempSum = zeros(K,1);
            for k = 1:K
                TempSum(k,1) = f(i,k)*C(i,1)*I_avg_f(i,k) + b(i,k)*C(i,1)*I_avg_b(i,k);
            % summatory over K of eqn. (2) in paper by Chung and Kennedy
            end
            ROIclean(i,1) = sum(TempSum);
        end
        
        ROI(2:ltime-1,y) = ROIclean(2:ltime-1,1);
        
    end

ROI(ROI==Inf)=500;
ROI(isnan(ROI))=0;

end

end


%% --- SIGNAL-TO-NOISE RATIO --- %

SNR_ROI=zeros(size(ROI));
sROI=zeros(size(ROI));
for g=1:size(ROI,2)
    for j=11:size(ROI,1)
        SNR_ROI(j,g)=((ROI(j,g))-(mean(ROI(j-10:j-1,g))))/std(ROI(j-10:j-1,g));
    end
    sROI(:,g)=smooth(ROI(:,g),10);
end


%% --- DERIVATING THE SIGNAL (local slope) --- %

lin=@(l,w)l(1)+l(2)*w;
l0=[0,0];
options=optimoptions(@lsqcurvefit,'MaxFunEvals',1000,'MaxIter',1000,'Algorithm','levenberg-marquardt');
lb=[];
ub=[];
window=7;
w=(1:1:window)';
q=(window-1)/2;
derivROI=zeros(size(F));
s_derivROI=zeros(size(F));

for h=1:size(ROI,2)
    for i=(q+1):(size(ROI,1)-q)
        FitLin=polyfit(w,ROI((i-q):(i+q),h),1);
        derivROI(i,h)=FitLin(1);
    end
    s_derivROI(:,h)=smooth(derivROI(:,h),3);
end


%% --- SINGLE VESICLE EVENTS FINDING and FITTING --- %

AmpliReal=zeros(size(ROI));
dF_multi=zeros(Total_points,1);
dF_single=zeros(Total_points,1);
time_peak_a=0;
time_peak_b=0;
a=1;
b=1;
Tau_single=0;
Dwell_single=0;
A_single=1;
Ampli_single=1;
Tau_multi=0;
Dwell_multi=0;
A_multi=1;
Ampli_multi=1;
GofFit_a=1;
GofFit_b=1;
search_dwell=Total_points-16-min_fit_points;
expdecay=@(x,tfit)x(1)*exp(-x(2)*tfit);
x0=[0,0];
options=optimoptions(@lsqcurvefit,'MaxFunEvals',1000,'MaxIter',1000,'Algorithm','levenberg-marquardt');
lb=[];
ub=[];
p=1;
NO_Pr=0;

for k=1:size(ROI,2)
    for i=iini-16:ifin
        stop=0;
        avF=mean(ROI(i-15:i-1,k));
        sdF=std(ROI(i-15:i-1,k),1);
            if (abs(rem((t(i)-time),tinterval))<=1.5) && (t(i)>(time-3))
                if (ROI(i,k)>(avF+3*sdF)) && (mean(AmpliReal(i-16:i-1,k))==0) && ((ROI(i,k)-avF)>=15) 
                    if (AmpliReal(i,k)<Amplitude_limit)
                        for u=i:i+search_dwell
                            endlimit=u+min_fit_points;
                            if ((mean(ROI(i:u,k)))>(avF+2*sdF)) && (derivROI(u+1,k)<=(-1)) % dwell time finishes when derivetive is negative
                                for e=i+Total_points-16:-3:endlimit
                                    tfit=t(u:e)-t(u);
                                    ffit=ROI(u:e,k)-avF;
                                    Fit=lsqcurvefit(expdecay,x0,tfit,ffit,lb,ub,options);
                                    Tau1=1./(Fit(2));
                                    FitFun=Fit(1)*exp(-Fit(2)*tfit);
                                    GofFit_a(a)=100*goodnessOfFit(FitFun,ffit,'NRMSE');
                                    if (Tau1>increment) && (Tau1<tinterval) && ((mean(ROI(e-3:e,k)))<(avF+sdF)) && (GofFit_a(a)>50)
                                        Tau_single(a)=Tau1;
                                        Dwell_single(a)=t(u)-t(i);
                                        Ampli_single(a)=Fit(1);
                                        A_single(a)=ROI(i,k)-avF;
                                        dF_single(:,a)=ROI(i-15:i+Total_points-16,k)-avF;
                                        time_peak_a(a)=t(i);
                                        AmpliReal(i,k)=ROI(i,k)-avF;
                                        a=a+1;
                                        stop=1;
                                        break
                                    end
                                end
                            elseif (u==(i+search_dwell)) && (mean(ROI(i:u,k))>(avF+2*sdF)) 
                                Tau_single(a)=NaN;
                                Dwell_single(a)=Inf;
                                Ampli_single(a)=NaN;
                                A_single(:,a)=ROI(i,k)-avF;
                                dF_single(:,a)=ROI(i-15:i+Total_points-16,k)-avF;
                                time_peak_a(a)=t(i);
                                AmpliReal(i,k)=ROI(i,k)-avF;
                                a=a+1;
                                stop=1;
                            end
                            if stop==1
                                break
                            end
                        end
                    elseif (AmpliReal(i,k)>Amplitude_limit)
                        for u=i:i+search_dwell
                            endlimit=u+min_fit_points;
                            if ((mean(ROI(i:u,k)))>(avF+2*sdF)) && (derivROI(u+1,k)<=(-1)) 
                                for e=i+Total_points-16:-3:endlimit
                                    tfit=t(u:e)-t(u);
                                    ffit=ROI(u:e,k)-avF;
                                    Fit=lsqcurvefit(expdecay,x0,tfit,ffit,lb,ub,options);
                                    Tau1=1./(Fit(2));
                                    FitFun=Fit(1)*exp(-Fit(2)*tfit);
                                    GofFit_b(b)=100*goodnessOfFit(FitFun,ffit,'NRMSE');
                                    if (Tau1>increment) && (Tau1<tinterval) && ((mean(ROI(e-3:e,k)))<(avF+sdF)) && (GofFit_b(b)>50) 
                                        Tau_multi(b)=Tau1;
                                        Dwell_multi(b)=t(u)-t(i);
                                        Ampli_multi(b)=Fit(1);
                                        A_multi(b)=ROI(i,k)-avF;
                                        dF_multi(:,b)=ROI(i-15:i+Total_points-16,k)-avF;
                                        time_peak_b(b)=t(i);
                                        AmpliReal(i,k)=ROI(i,k)-avF;
                                        b=b+1;
                                        stop=1;
                                        break
                                    end
                                end
                            elseif (u==(i+search_dwell)) && (mean(ROI(i:u,k))>(avF+2*sdF)) 
                                Tau_multi(b)=NaN;
                                Dwell_multi(b)=Inf;
                                Ampli_multi(b)=NaN;
                                A_multi(:,b)=ROI(i,k)-avF;
                                dF_multi(:,b)=ROI(i-515:i+Total_points-16,k)-avF;
                                time_peak_b(b)=t(i);
                                AmpliReal(i,k)=ROI(i,k)-avF;
                                b=b+1;
                                stop=1;
                            end
                            if stop==1
                                break
                            end
                        end
                    end
                end
            end
    end
end

A_single=A_single(A_single>0);
meanA_single=mean(A_single);
medianA_single=median(A_single);
meanTau_single=nanmean(Tau_single);
medianTau_single=nanmedian(Tau_single);
meanDwell_single=mean(Dwell_single(Dwell_single<Inf));
medianDwell_single=median(Dwell_single(Dwell_single<Inf));
meanAmpli_single=nanmean(Ampli_single);
mean_dF_single=mean(dF_single,2);

A_multi=A_multi(A_multi>0);
meanA_multi=mean(A_multi);
medianA_multi=median(A_multi);
meanTau_multi=nanmean(Tau_multi);
medianTau_multi=nanmedian(Tau_multi);
meanDwell_multi=mean(Dwell_multi(Dwell_multi<Inf));
medianDwell_multi=median(Dwell_multi(Dwell_multi<Inf));
meanAmpli_multi=nanmean(Ampli_multi);
mean_dF_multi=mean(dF_multi,2);

Amplitudes=horzcat(A_single(A_single>1),A_multi(A_multi>1));
Taus=horzcat(Tau_single(Tau_single>0),Tau_multi(Tau_multi>0));
Dwells=horzcat(Dwell_single,Dwell_multi);
time_peak=horzcat(time_peak_a(time_peak_a>0),time_peak_b(time_peak_b>0));

% For release probability we consider all events:
Pr_perROI=sum(AmpliReal~=0,1)/Nstim;                                 


%% --- STRONG STIMULATION PEAKS FITTING %

% BigTau >> decay time constant
% BigA >> amplitude

expdecay3=@(z,tBig)z(1)*exp(-z(2)*tBig)+z(3); 
z0=[0,0,0];
lindecay3=@(l,tBig)l(1)*tBig+l(2);
l0=[0,0];
options=optimoptions(@lsqcurvefit,'MaxFunEvals',1000,'MaxIter',1000,'Algorithm','levenberg-marquardt');
lb=[];
ub=[];

BigTau1=0;
BigA1=0;
BigSlope1=0;

for g=1:size(ROI,2)
    if tpeak2>0
    tBig1=t(ipeak1:ipeak2);
    BigExpdecay1=lsqcurvefit(expdecay3,z0,tBig1,ROI(ipeak1:ipeak2,g),lb,ub,options);
    FitFun_Big=BigExpdecay1(1)*exp(-BigExpdecay1(2)*tBig1)+BigExpdecay1(3);
    GofFit_Big(g)=100*goodnessOfFit(FitFun_Big,ROI(ipeak1:ipeak2,g),'NRMSE');
    BigTau1(g)=1./(BigExpdecay1(2));
    BigA1(g)=mean(ROI(ipeak1:ipeak1+19,g))-mean(ROI(ipeak1-225:ipeak1-200,g));
    BigLindecay1=lsqcurvefit(lindecay3,l0,tBig1,ROI(ipeak1:ipeak2,g),lb,ub,options);
    BigSlope1(g)=BigLindecay1(1);
    end
end

BigTau1=BigTau1(BigTau1>0);
meanBigTau1=mean(BigTau1);
BigSlope1=BigSlope1(BigSlope1<0);
meanSlopeTau1=mean(BigSlope1);
meanBigA1=mean(BigA1);


%% --- finally, SAVING the data!! --- %

save(out_filename);
 
% --> and cleaning the memory to start again:
clearvars -except in_filename out status sheets NSheets sheet dnRounds RBall CKF FFTF BBC

end
end

clear