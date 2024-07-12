function mainL1(training,predictions,C,N,eps,epsilon,kernel, param_kernel, lambda_0, max_iter,tol)

    Q=zeros(N,N); % inizializzazione matrice
    Q_piccola=zeros(N/2,N/2); % inizializzazione sottomatrice che definisce la mia funzione quadratica
    if strcmp(kernel, 'Polinomiale')
        for i=1:N/2 % sfrutto simmetria matrice
            for j=1:i-1
                Q_piccola(i,j)=PolyKer(training(i),training(j),param_kernel); % PolyKer funzione kernel polinomiale
            end
        end
        Q_piccola=Q_piccola+ Q_piccola';

        for i = 1 : N/2 % calcolo la diagonale
            Q_piccola(i,i)=PolyKer(training(i),training(i),param_kernel);
        end
        Q=[Q_piccola,-Q_piccola;-Q_piccola, Q_piccola];

    elseif strcmp(kernel, 'Radiale')

        for i=1:N/2 % sfrutto simmetria matrice
            for j=1:i-1
                Q_piccola(i,j)=RadialKer(training(i),training(j),param_kernel);
            end
        end
        Q_piccola=Q_piccola+ Q_piccola'; 

        for i = 1 : N/2 %calcolo la diagonale
            Q_piccola(i,i)=RadialKer(training(i),training(i),param_kernel);
        end

        Q=[Q_piccola,-Q_piccola;-Q_piccola, Q_piccola]; 

    else
        fprintf('Tipo di kernel non supportato\n');
        return;
    end 

    q=zeros(1, N);  % Inizializzazione del vettore funzione quadratica
    first_indices = 1:(N/2);
    second_indices =((N/2)+1):N;
    q(first_indices) = eps - predictions;
    q(second_indices) = eps + predictions;

    
    % chiamata algoritmo risolutivo -> resituisce due valori:
    % 1) minimo della funzione 2) iterazioni eseguite

    [sol, iter_true]= PGML1(Q,q,epsilon,C,lambda_0,N,max_iter,tol);
    
    % definizione vincoli per chiamata a Quadprog()
    lb=zeros(N,1);
    ub=C*ones(N,1);
    A_eq=[ones(1,N/2),-ones(1,N/2)];
    b_eq=0;
    A=[];
    b=[];

    options=optimoptions('quadprog','Display','off', 'OptimalityTolerance',10^-14,'ConstraintTolerance',10^-8);
    tic;
    [min2,f_min2]=quadprog(Q,q',A,b,A_eq,b_eq,lb,ub,lambda_0',options); % https://it.mathworks.com/help/optim/ug/quadprog.html#bssh6y6-1_sep_shared-Aeq
    toc;
    
    
    fprintf('|lambda_0-min*|=%d\n',norm(lambda_0-min2'));

    % Definizione grafico errore

    errore=(sol-f_min2)/max(abs(f_min2),1); % definizione errore relativo
    x=1:(iter_true);
    semilogy(x, errore, 'r', 'DisplayName', sprintf('Convergenza step size=1/L')); % errore in funzione del numero di iterazioni
    legend('show');

    % Definizione grafico errore Torema 2 convergenza teorica

    %errore_x=norm(lambda_0-min2').^2;
    %L=eig(Q);
    %L=L(end);
    %errore_x=(L./(2*x))*errore_x;
    %semilogy( x,errore_x,'r','DisplayName', sprintf('L/2k * ||(x_0-x*)||^2')); % errore usando un passo costante pari a 1/max_autovalore
    
end