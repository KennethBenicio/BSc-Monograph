function [NMSE,Capacity] = IRS2(snr_db,P,N)
    Mt = 2;
    Mr = 2;
    K = N*Mr;
    snr_linear = 10^(snr_db/10);   
    var_signal = 1;
    variance_noise = var_signal/snr_linear;

%-----------------------PARTE 1: Gerando sinal Y---------------------------%
    
    %Inicializando os ângulos da IRS por uma dft KxK
    zeta_nk = dftmtx(K);
    Zeta_matrices = cell(P,1);
    for pp = 1:P
        Zeta_matrices{(pp)} = zeta_nk(:,1:N)';
    endfor

    alpha_irs = cell(K,1);
    alpha_irs(:) = zeros(P,P);
    Airs_matrices = cell(P,1);
    Birs_matrices = cell(P,1);
    for pp = 1:P
        irs_pos = 0:1:(N-1);
        
        AOA_irs = pi.*rand(1,1) - pi/2;
        airs = 1 .* exp(-j*pi*(irs_pos*cos(AOA_irs)))';
        airs = airs./abs(airs);
        Airs_matrices{(pp)} = airs;
        
        AOD_irs = pi.*rand(1,1) - pi/2;
        birs = 1 .* exp(-j*pi*(irs_pos*cos(AOD_irs)))';
        birs = birs./abs(birs);
        Birs_matrices{(pp)} = birs;
        
        for kk = 1:K
            alpha_irs{(kk)}(pp,pp) = airs.'*diag(Zeta_matrices{(pp)}(:,kk))*birs;
        endfor
    endfor

    alpha = (1/sqrt(2))*(randn(P,P) + 1j*randn(P,P));
    alpha = diag(diag(alpha));
    beta = (1/sqrt(2))*(randn(P,P) + 1j*randn(P,P));
    beta = diag(diag(beta));
    Arx_matrices = zeros(Mr,P);
    Btx_matrices = zeros(Mt,P);
    for pp = 1:P
        AOA_rx = pi.*rand(1,1) - pi/2;
        arx_pos = 0:1:(Mr-1);
        arx = exp(-j*pi*(arx_pos*cos(AOA_rx)))';
        Arx_matrices(:,pp) = arx;
        
        AOD_tx = pi.*rand(1,1) - pi/2;
        brx_pos = 0:1:(Mt-1);
        btx = exp(-j*pi*(brx_pos*cos(AOD_tx)))'; 
        Btx_matrices(:,pp) = btx;
    endfor

    X = randi(4,[Mt K]) - 1;
    X = pskmod(X,4,0,"gray");
    Y = cell(P,K);
    for pp = 1:P
        for kk = 1:K
            noise = (sqrt(variance_noise)/sqrt(2))*(randn(Mr,1) + 1j*randn(Mr,1));
            aux = Arx_matrices(:,pp)*alpha(pp,pp)*alpha_irs{(kk)}(pp,pp)*beta(pp,pp)*Btx_matrices(:,pp).'*X(:,kk);
            Y(pp,kk) = aux + noise; 
        endfor
    endfor
    
    %--------------------------PARTE 2: Estimação------------------------------%
    
    %Inicializando as celulas que serão utilizadas
    C = cell(P,1);
    Z = cell(P,1);
    v_hat = cell(P,1);
    zeta_opt = cell(P,1);
    
    W = cell(P,1);
    Q = cell(P,1);
    NMSE = cell(P,1);
    Hefc = zeros(Mr,Mt);
    for pp = 1:P
        c_aux = kr(Zeta_matrices{(pp)},X);
        c_aux = kron(c_aux,eye(Mr)).';
        C{(pp)} = c_aux;

        %Obtendo o vetor y da p-esima irs
        y = (cell2mat(Y(pp,:)'));

        %Conferindo a igualdade entre pinv(C)*y = kr(btx*birs',arx*airs')
        H = Btx_matrices(:,pp)*Birs_matrices{(pp)}.';
        G = Arx_matrices(:,pp)*Airs_matrices{(pp)}.';
        v = alpha(pp,pp)*beta(pp,pp)*vec(kr(H,G));
        teste = real(y./(c_aux*v));

        Z{(pp)} = reshape(pinv(C{(pp)})*y,[Mr*Mt,N]);
        
        %Estimando as componentes de Z
        [U,S,V] =  svd(Z{(pp)},"econ");
        U_p = U(:,1)*S(1,1);
        V_p = conj(V(:,1));
        %NMSE em Z na p-ésima irs
        v_hat{(pp)} = vec(reshape(U_p, [Mt*Mr, 1])*reshape(V_p, [1, N]));
        NMSE{(pp)} = norm(v - v_hat{(pp)}).**2/norm(v).**2;
        
        zeta_opt{(pp)} = exp(-angle(V_p));
    endfor

    %NMSE na 1-irs
    NMSE = cell2mat(NMSE);
    NMSE = NMSE(1,1);
    
    X = randi(4,[P 1]) - 1;
    X = (pskmod(X,4,0,"gray")./sqrt(2*P))';
    %Rxx = X*(conj(X).')/100;
    
    %Capacity
    Hefc = zeros(P,P);
    for pp = 1:P
      Hefc(pp,pp) = Airs_matrices{(pp)}.'*diag(zeta_opt{(pp)})*Birs_matrices{(pp)};
    endfor
    Hefc = Arx_matrices*alpha*Hefc*beta*Btx_matrices.'*X;
    Capacity = abs(log2(det(eye(Mr) + (Hefc*(conj(Hefc)'))./(Mt*variance_noise))));
endfunction