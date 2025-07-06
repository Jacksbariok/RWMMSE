function [hist] = genChannel(L,K,N,M)

NM = N*M;
%% H_lk C^{NxM} 拼起来的分布服从 CN(0,βI)，其中 β=128.1+37.6lg(d)[dB]
rng(209)
AP_pos = zeros(L,2);
for l = 1:L
    AP_pos(l,:) = [rand, rand] * 500;
end
UE_pos = zeros(K,2);
dis = zeros(L,K);

rng('shuffle')
for k = 1:K
    t = [rand,rand]*500;
    dt = AP_pos - t;
    distance = sqrt(sum(dt.^2,2));
    while (sum((distance<100)) ~= 0) || (sum((distance>300)) ~= 0)
        t = [rand,rand]*500;
        dt = AP_pos - t;
        distance = sqrt(sum(dt.^2,2));
    end
    dis(:,k) = distance;
    UE_pos(k,:) = t;
end


%% simulate channel matrix H
H = zeros(N,M,L,K);
Pathloss = 128.1 + 37.6 *log10(dis/1000);
Channelgain = 10.^(-Pathloss/10);
for l = 1:L
    for k = 1:K
        H(:,:,l,k) = reshape(sqrt(Channelgain(l,k)/2) *(randn(NM,1) + 1i*randn(NM,1)),[N,M]);
    end
end

%% history results
hist.AP_pos = AP_pos;
hist.UE_pos = UE_pos;
hist.distance = dis;
hist.H = H;
end