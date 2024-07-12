function [d] = compute_direction(m,grad_lambda,a,ind_C,ind_0,ind_wc,N)
    d(ind_wc)=-(grad_lambda(ind_wc)+m*a(ind_wc)); % no vincolo
    for j=1:N
        if ind_C(j)==1 && j<=N/2 %vincolo (-infty,0]
            if  m>-grad_lambda(j)
                d(j)=-(grad_lambda(j)+m*a(j));
            else
                d(j)=0;
            end
        elseif ind_C(j)==1 && j>N/2
            if m<grad_lambda(j)
                d(j)=-(grad_lambda(j)+m*a(j));
            else
                d(j)=0;
            end
        end
        if ind_0(j)==1 && j<=N/2 %vincolo [0,+infty)
            if m<-grad_lambda(j)
                d(j)=-(grad_lambda(j)+m*a(j));
            else
                d(j)=0;
            end
        elseif ind_0(j)==1 && j>N/2
            if m>grad_lambda(j)
                d(j)=-(grad_lambda(j)+m*a(j));
            else
                d(j)=0;
            end
        end  
    end
end