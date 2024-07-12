
% Routine che risolve il  problema CQKnP (sottoproblema 1.5 report)

function d_star=AP(grad_lambda,lambda_n,N, C, tol)
    a=[ones(1,N/2),-ones(1,N/2)]; % vettore a(i)
    ind_C=C-lambda_n<=tol;% indici del punto lambda corrente che hanno componente uguale a C -> avranno vincolo (-infty,0]
    ind_0=lambda_n<=tol;% indici del punto lambda corrente che hanno componente uguale a 0 -> avranno vincolo [0,+infty)
    ind_wc=~ind_C&~ind_0;% indici del punto lambda corrente che non hanno vincoli
    indici=ind_C|ind_0;% indici dei punti che hanno un vincolo (che sia C o 0 è indifferente)
    brk=zeros(1,N);% array relativo ai break points
    
    % Definizione array breakpoints
    n_brk=0;
    for i=1:N
        if indici(i)==1 && i<=N/2
            n_brk=n_brk+1;
            brk(n_brk)=-grad_lambda(i);
        elseif indici(i)==1 && i>N/2
            n_brk=n_brk+1;
            brk(n_brk)=grad_lambda(i);
        end
    end

    brk=brk(1,1:n_brk);
    brk_sort=unique(brk,'sorted');% ordinamento e rimozione eventuali duplicati
    n=n_brk;
    last_psi_der=0; % penultima derivata
    
    for i=1:n
        m=brk_sort(i);
        d=compute_direction(m,grad_lambda,a,ind_C,ind_0,ind_wc,N);
        psi_der=sum(a.*d); % derivata funzione corrispondente al breakpoint "m" corrente 
        if abs(psi_der)<=tol
            d_star=d;
            break
        end
        if i==1 && psi_der<-tol % caso con derivata iniziale negativa
            m_left=m;
            while psi_der<0
                m_old=m_left;
                m_left=m_left-1;
                d=compute_direction(m_left,grad_lambda,a,ind_C,ind_0,ind_wc,N);
                psi_der=sum(a.*d);
            end
            % interpolazione tra i punti m_old e m_left
            m_half=(m_left+m_old)/2;
            int_half=[m_left,m_half,m_old];
            d=compute_direction(m_half,grad_lambda,a,ind_C,ind_0,ind_wc,N);
            psi_der=sum(a.*d);% derivata nel punto medio "m_half"
            while abs(psi_der)>=tol
                if psi_der<-tol
                    int_half(3)=m_half;
                    m_half=(int_half(1)+m_half)/2;
                    int_half(2)=m_half;
                elseif psi_der>tol
                    int_half(1)=m_half;
                    m_half=(m_half+int_half(3))/2;
                    int_half(2)=m_half;
                end
                d=compute_direction(m_half,grad_lambda,a,ind_C,ind_0,ind_wc,N);
                psi_der=sum(a.*d);
            end
            d_star=d;
            break
        end
        if i==n && psi_der>tol % caso con derivata finale positiva
            m_right=m;
            while psi_der>0
                m_old=m_right;
                m_right=m_right+1;
                d=compute_direction(m_right,grad_lambda,a,ind_C,ind_0,ind_wc,N);
                psi_der=sum(a.*d);
            end
            % interpolazione tra i punti m_old e m_right
            m_half=(m_old+m_right)/2;
            int_half=[m_old,m_half,m_right];
            d=compute_direction(m_half,grad_lambda,a,ind_C,ind_0,ind_wc,N);
            psi_der=sum(a.*d);% derivata nel punto medio "m_half"
            while abs(psi_der)>=tol
                if psi_der<-tol
                    int_half(3)=m_half;
                    m_half=(int_half(1)+m_half)/2;
                    int_half(2)=m_half;
                elseif psi_der>tol
                    int_half(1)=m_half;
                    m_half=(m_half+int_half(3))/2;
                    int_half(2)=m_half;
                end
                d=compute_direction(m_half,grad_lambda,a,ind_C,ind_0,ind_wc,N);
                psi_der=sum(a.*d);
            end
            d_star=d;
            break
        end
        if last_psi_der>tol && psi_der<-tol % caso con cambiamento segno derivata
            m_half=(brk_sort(i)-brk_sort(i-1))/2;
            int_half=[brk_sort(i-1),m_half, brk_sort(i)];
            d=compute_direction(m_half,grad_lambda,a,ind_C,ind_0,ind_wc,N);
            psi_der=sum(a.*d); % derivata nel punto medio "m_half"
            while abs(psi_der)>=tol
                if psi_der<-tol
                    int_half(3)=m_half;
                    m_half=(int_half(1)+m_half)/2;
                    int_half(2)=m_half;
                elseif psi_der>tol
                    int_half(1)=m_half;
                    m_half=(m_half+int_half(3))/2;
                    int_half(2)=m_half;
                end
                d=compute_direction(m_half,grad_lambda,a,ind_C,ind_0,ind_wc,N);
                psi_der=sum(a.*d);
            end
            d_star=d;
            break
        end
        last_psi_der=psi_der;
    end
end

                

  
