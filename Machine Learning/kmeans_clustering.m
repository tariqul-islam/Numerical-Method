function g_u = kmeans_clustering(X,K,SPACE)
%u=kmeans_clustering(X,K,SPACE)
%X
%   Each vector in row
%
%K
%   Number of clusters
%
%SPACE
%   The algorithm chooses initial values of u
%   randomly. Hence to know where to choose from
%   the input space is required. [lower_limit upper_limit]
%

%written by
%tariqul islam ponir
%ponir(dot)bd(at)hotmail(dot)com
%

    [nx,dx]=size(X);
    [ns,ds]=size(SPACE);
    [nk,dk]=size(K);
    if ~(nk==1 && dk==1)
        error('Multiple K values are not supported yet.');
    elseif K==0
        error('K has to be more than zero.');
    end
    if ~(ns==1 && ds==2)
        error('Different sample region for different dimensions is not supported yet.');
    end
    if (nx<K)
        error(['Cannot cluster ' num2str(nx) ' points in ' num2str(K) ' cluster']);
    end
    
    u=zeros(K,dx);
    SK = [];
    S=zeros(K,nx);
    
    function update_cluster()
        SK=zeros(1,K);
        for i=1:nx
            uc_min=0;
            uc_minval=10^20;
            for k=1:K
                val = norm(X(i,:)-u(k,:));
                if val<uc_minval
                    uc_minval=val;
                    uc_min=k;
                end
            end
            SK(uc_min)=SK(uc_min)+1;
            S(uc_min,SK(uc_min))=i;
        end
    end

    function update_u()
        for k=1:K
            uu_i=SK(k);
            if uu_i==0
                continue;
            end
            sum=zeros(1,dx);
            for i=1:uu_i
                uu_j=S(k,i);
                sum=sum+X(uu_j,:);
            end
            u(k,:)=sum/uu_i;
        end
    end

    function z = check_value()
        z=0;
        for k=1:K
            cv_i=SK(k);
            for i=1:cv_i
                cv_j=S(k,i);
                z=z+norm(X(cv_j,:)-u(k,:))^2;
            end
        end
    end

    g_u=[];
    g_min=10^30;
    for gi=1:100
        u=SPACE(1)+rand(K,dx)*(SPACE(2)-SPACE(1));
        update_cluster();

        chk_val=check_value();
        old_chk_val=10^30;

        iter=0;
        while abs(chk_val-old_chk_val)>10^-10
            update_u();
            update_cluster();
            old_chk_val=chk_val;
            chk_val=check_value();
            if iter==1000
                warning('maximum iterations reached');
                break;
            end
            iter=iter+1;
        end
        
        %disp(chk_val);
        if chk_val<g_min
            g_u=u;
            g_min=chk_val;
            
        end
    end
    %iter=iter+2;
end
