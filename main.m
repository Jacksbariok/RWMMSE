clc
clear variables

%% 设置仿真参数

L = 4;
M = 64;
K = 16;
N = 4;
D = 4; % the number of data stream
SNR = 0;

hist = genChannel(L,K,N,M);
H = hist.H;


noise_pow = 0;
for k = 1:K
    sum_normH = 0;
    for l = 1:L
        sum_normH = sum_normH + norm(H(:,:,l,k),'fro')^2;
    end
    noise_pow = noise_pow + log10(sum_normH/N);
end


sigma_ue = sqrt(10^(noise_pow/K) * 10^(-SNR/10));

P_max = 1;

distance = hist.distance;
serve_AP = zeros(L,K);
for k = 1:K
    [~,index_L] = sort(distance(:,k));
    serve_AP(index_L(1:2),k) = 1;
end
%% RWMMSE Algorithm

% Step 1: construct \bar{H}_l = H_l * H_l^H

H_l = zeros(K*N,M,L);
H_bar = zeros(K*N,K*N,L);
for l = 1:L
    tmp = zeros(K*N,M);
    for k = 1:K
        tmp((k-1)*N+1:k*N,:) = H(:,:,l,k);
    end
    H_l(:,:,l) = tmp;
    H_bar(:,:,l) = tmp * tmp';
end

H_tilde = zeros(N,L*K*N,K);
for k = 1:K
    for l = 1:L
        H_tilde(:,(l-1)*K*N+1:l*K*N,k) = H(:,:,l,k) * H_l(:,:,l)';
    end
end


V = zeros(M,D,L,K);

W = zeros(D,D,K);
for k = 1:K
    W(:,:,k) = eye(D);
end

X = zeros(L*K*N,K*D);
for l = 1:L
    for k = 1:K
        [~,~,Vs] = svd(H(:,:,l,k));
        V(:,:,l,k) = Vs(:,1:D);
        X((l-1)*K*N+1:l*K*N,(k-1)*D+1:k*D) = pinv(H_l(:,:,l)') * Vs(:,1:D);
    end
    power = sum(V.*conj(V), [1,2,4]);
    V(:,:,l,:) = V(:,:,l,:) * sqrt(P_max)/power(l);

    X((l-1)*K*N+1:l*K*N,:) = P_max * X((l-1)*K*N+1:l*K*N,:)/trace(H_bar(:,:,l)*X((l-1)*K*N+1:l*K*N,:) * X((l-1)*K*N+1:l*K*N,:)');
end

W_old = W;
U = zeros(N,D,K);
max_iter = 20;

sumHV = get_sumHV(H, V);
oldSumRate = get_sumrate(sumHV, sigma_ue);
sumRates = [oldSumRate];
for i = 1:max_iter 
    %% Record W_k
    W_old = W;
    
    %% Updata U_k
    XX = zeros(L*K*N);
    for m = 1:K
        XX = XX + X(:,(m-1)*D+1:m*D) * X(:,(m-1)*D+1:m*D)';
    end
    for k = 1:K
        U(:,:,k) = pinv(H_tilde(:,:,k) * XX * H_tilde(:,:,k)' + sigma_ue^2*eye(N) ) * H_tilde(:,:,k) * X(:,(k-1)*D+1:k*D);
    end

    %% Upadte W_k
    for k = 1:K
        W(:,:,k) = pinv(eye(D) - U(:,:,k)' * H_tilde(:,:,k) * X(:,(k-1)*D+1:k*D));
    end


    %% Update X_lk
    for l = 1:L
        HUWUH = zeros(K*N,K*N);
        for m = 1:K
            HUWUH = HUWUH + H_bar((m-1)*N+1:m*N,:,l)' * U(:,:,m) * W(:,:,m) * U(:,:,m)' * H_bar((m-1)*N+1:m*N,:,l);
        end
        
        B = zeros(K*N,D,K);
        for k = 1:K
            
            HUW = zeros(K*N,D);
            for m = 1:K
                HUW = HUW + H_bar((m-1)*N+1:m*N,:,l)' * U(:,:,m) * W(:,:,m) * U(:,:,m)' * (H_tilde(:,:,m) * X(:,(k-1)*D+1:k*D)-H_bar((m-1)*N+1:m*N,:,l)*X((l-1)*K*N+1:l*K*N,(k-1)*D+1:k*D));
            end
            B(:,:,k) = (H_bar((k-1)*N+1:k*N,:,l)' * U(:,:,k) * W(:,:,k) - HUW);
            X((l-1)*K*N+1:l*K*N,(k-1)*D+1:k*D) = pinv(HUWUH) * B(:,:,k);
        end
        
        powl = trace(H_bar(:,:,l) * X((l-1)*K*N+1:l*K*N,:) * X((l-1)*K*N+1:l*K*N,:)');
        if powl > P_max
            low = 0;
            high = 10;

            for k = 1:K
                X((l-1)*K*N+1:l*K*N,(k-1)*D+1:k*D) = pinv(HUWUH + high*H_bar(:,:,l)) * B(:,:,k);
            end
            powl = trace(H_bar(:,:,l)*X((l-1)*K*N+1:l*K*N,:) * X((l-1)*K*N+1:l*K*N,:)');
           
            while powl > P_max
                low = high;
                high = 10*high;

                for k = 1:K
                    X((l-1)*K*N+1:l*K*N,(k-1)*D+1:k*D) = pinv(HUWUH + high*H_bar(:,:,l)) * B(:,:,k);
                end
                powl = trace(H_bar(:,:,l)*X((l-1)*K*N+1:l*K*N,:) * X((l-1)*K*N+1:l*K*N,:)');
            end

            mid = (low + high)/2;

            for k = 1:K
                X((l-1)*K*N+1:l*K*N,(k-1)*D+1:k*D) = pinv(HUWUH + mid*H_bar(:,:,l)) * B(:,:,k);
            end
            powl = trace(H_bar(:,:,l)*X((l-1)*K*N+1:l*K*N,:) * X((l-1)*K*N+1:l*K*N,:)');

            while abs((powl - P_max)/P_max) > 1e-5 && high-low > 1e-10
                 if powl > P_max
                     low = mid;
                 else
                     high = mid;
                 end
                 mid = (low + high)/2;

                 for k = 1:K
                     X((l-1)*K*N+1:l*K*N,(k-1)*D+1:k*D) = pinv(HUWUH + mid*H_bar(:,:,l)) * B(:,:,k);
                 end
                 powl = trace(H_bar(:,:,l)*X((l-1)*K*N+1:l*K*N,:) * X((l-1)*K*N+1:l*K*N,:)');
            end

            mu_l = (low + high)/2;

            for k = 1:K
                X((l-1)*K*N+1:l*K*N,(k-1)*D+1:k*D) = pinv(HUWUH + mu_l*H_bar(:,:,l)) * B(:,:,k);
                V(:,:,l,k) = H_l(:,:,l)' * X((l-1)*K*N+1:l*K*N,(k-1)*D+1:k*D);
            end
            
        end

    sumHV = get_sumHV(H, V);
    end

    sumRate = get_sumrate(sumHV, sigma_ue);
    sumRates = [sumRates, sumRate];

    fprintf("Reduce WMMSE, iter %2d, sum_rate: %f, %s\n", i, sumRate, datestr(now));

    if abs(sumRate - oldSumRate) / oldSumRate < 1e-8
        break
    end

    oldSumRate = sumRate;

end

% figure
% scatter(hist.AP_pos(:,1),hist.AP_pos(:,2),'>')
% hold on
% scatter(hist.UE_pos(:,1),hist.UE_pos(:,2),'o')
% legend('AP Positions','UE Positions')
% grid on

% figure
% plot(sumRates,'o-b','LineWidth',2)
% grid on
% legend('RWMMSE',Location='southeast')
% xlabel('迭代次数')
% ylabel('和速率（bpcu）')
% xlim([0,max_iter])