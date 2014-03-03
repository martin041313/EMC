%% the main program

function B = emc(I,wo,ho)
    B=1;
    %% initial begin
    
    % Convert input sRGB image to CIELab color space.
    if exist('applycform','file')
       I = applycform(I,makecform('srgb2lab'));
    else
       I = colorspace('Lab<-RGB',I);
    end
    
    [wi,hi,~] = size(I);
    rx = fix(wi/wo);
    ry = fix(hi/ho);
    all = wo*ho;
    u = zeros(2,1,all);
    v = zeros(3,1,all);
    u_sigma = zeros(2,2,all);
    v_sigma = zeros(all);
    for k = 1 : all
        xk = fix(mod(k,wo)+1/2);
        yk = fix(floor(k/wo)+1/2);
        u(:,:,k) =[xk,yk]';
        u_sigma(:,:,k) = [rx/3 0;0 ry/3];
        v(:,:,k) = [1/2 1/2 1/2]';
        v_sigma(k) = 0.0001;
    end 
    
    m_diff = 0; % 1:m阶段产生了新结果
    c_diff = 0; % 1:c阶段产生了新结果
    %% initial end
    
    w = zeros(hi,wi,all);
    while(1)
       %% e_step
       
       
       
       %% compute all kernels
       for k = 1 : all
           xk = fix(mod(k,wo)+1/2);
           yk = fix(floor(k/wo)+1/2);
           
           w_sum = 0;
           for yi = yk-2*ry+1:yk+2*ry
               for xi = xk-2*rx+1:xk+2*rx
                   if(xi>0 && xi<=wi && yi>0 && yi<=hi)
                   p_delta = ([xi,yi]'-u(k));
                   w(k,yi,xi) = exp(-1/2 * p_delta' / (u_sigma) * p_delta - norm(I(xi,yi)-v(k))/(2*v_sigma^2));
                   w_sum = w_sum+w(k,yi,xi);
                   end
               end
           end
           for yi = yk-2*ry+1:yk+2*ry
               for xi = xk-2*rx+1:xk+2*rx
                   if(xi>0 && xi<=wi && yi>0 && yi<=hi)
                    w(k,yi,xi) = w(k,yi,xi)/w_sum;
                   end
               end
           end
       end
       
       %% Normalize per pixel
       y = zeros(all,hi,wi);
       wn_sum = zeros(hi,wi);
       for yi = 1:hi
           for xi = 1:wi
               for k = 1: all
                   xk = fix(mod(k,wo)+1/2);
                   yk = fix(floor(k/wo)+1/2);
                   if( abs(xi-xk)<2*rx && abs(yi-yk) <2*ry)
                       wn_sum(yi,xi) = wn_sum(yi,xi)+w(k,yi,xi);
                   end
               end
           end
       end
       for yi = 1:hi
           for xi = 1:wi
               for k = 1: all
                   y(k,yi,xi) = w(k,yi,xi) / wn_sum(yi,xi);
               end
           end
       end
       
                       
       
               
       %% e_step end 
       %% m_step maximization step
       for k=1:all
           w_sum=0;
           for yi = yk-2*ry+1:yk+2*ry
               for xi = xk-2*rx+1:xk+2*rx
                   if(xi>0 && xi<=wi && yi>0 && yi<=hi)
                        w_sum = w_sum+y(k,yi,xi);
                   end
               end
           end
           
           u_sigma_t=0;
           for yi = yk-2*ry+1:yk+2*ry
               for xi = xk-2*rx+1:xk+2*rx
                   if(xi>0 && xi<=wi && yi>0 && yi<=hi)
                       p_delta = ([xi,yi]'-u(k));
                       u_sigma_t = u_sigma_t+y(k,yi,xi)*p_delta*p_delta';
                   end
                   
               end
           end
           u_sigma_t = u_sigma_t/w_sum;
           if(u_sigma(k)~=u_sigma_t)
               
               u_sigma(k)=u_sigma_t;
               m_diff=1;
           end
           
            u_t=0;
           for yi = yk-2*ry+1:yk+2*ry
               for xi = xk-2*rx+1:xk+2*rx
                   if(xi>0 && xi<=wi && yi>0 && yi<=hi)
                       u_t = u_t+y(k,yi,xi)*I(xi,yi)';
                   end
                   
               end
           end
           u_t = u_t/w_sum;
           if(u(k)~=u_t)
               
               u(k)=u_t;
               m_diff=1;
           end
           
           v_t=0;
           for yi = yk-2*ry+1:yk+2*ry
               for xi = xk-2*rx+1:xk+2*rx
                   if(xi>0 && xi<=wi && yi>0 && yi<=hi)
                       v_t = v_t+y(k,yi,xi)*[xi,yi]';
                   end
                   
               end
           end
           v_t = v_t/w_sum;
           if(v(k)~=v_t)
               
               v(k)=v_t;
               m_diff=1;
           end
       end
           
       %% m_step end 
       %% c_step corretion step
       
       % spatial constraints
       u_avg=zeros(2,1);
       
       for k=1:all
           xk = fix(mod(k,wo)+1/2);
           yk = fix(floor(k/wo)+1/2);
           len = [-1,1,-wo,wo];
           count=0;
           for i=1:4
               t=k+len(i);
               if(t>0 && t<=all)
                   u_avg = u_avg+u(t);
                   count=count+1;
               end
               u_avg = u_avg/count;
               u_avg = clampBox((u(k)+u_avg)/2,[xk,yk]'+[rx/4,ry/4]',[xk,yk]'-[rx/4,ry/4]');

               if(u(k)~=u_avg)
                   u(k)=u_avg
                   c_diff=1;
               end
           end
       end
       
       %constant spatial variance
       
       for k=1:all
           [U,S,V] = svd(u_sigma(k));
           S(1,1) = clamp(S(1,1),0.05,0.1);
           S(2,2) = clamp(S(2,2),0.05,0.1);
           t = U*S*V;
           if(t~=u_sigma(k))
               c_diff=1;
               u_sigma(k)=t;
           end
           
       end
       
       % shape constraints
       for k=1:all
           xk = fix(mod(k,wo)+1/2);
           yk = fix(floor(k/wo)+1/2);
           len = [-1,1,-wo,wo,-wo-1,-wo+1,wo-1,wo+1];
           for i=1:8
               t = k+len(i);
               xn = fix(mod(t,wo)+1/2);
               yn = fix(floor(t/wo)+1/2);
               d = [xn-xk,yn-yk]';
               % directional variance
               s=0;
               for yi = yk-2*ry+1:yk+2*ry
                 for xi = xk-2*rx+1:xk+2*rx
                   if(xi>0 && xi<=wi && yi>0 && yi<=hi)
                       s = s+y(k,yi,xi)*max(0,([xi,yi]'-u(k))'*d)^2;
                   end
                 end
               end
               f=0;
               for yi = yk-2*ry+1:yk+2*ry
                 for xi = xk-2*rx+1:xk+2*rx
                   if(xi>0 && xi<=wi && yi>0 && yi<=hi)
                       f = f+y(k,yi,xi)*y(t,yi,xi);
                   end
                 end
               end
               % edge orientation 
               o = 0;
               for yi = yk-2*ry+1:yk+2*ry
                 for xi = xk-2*rx+1:xk+2*rx
                   if(xi>0 && xi<=wi && yi>0 && yi<=hi)
                       
                   end
                 end
               end
               
               
              
               
           
       
       %% c_step end
       if(m_diff==0 && c_diff==0)
           break
       end
       
           end
       end
       
    end
    
    
    % Convert filtered image back to sRGB color space.
    if exist('applycform','file')
       B = applycform(B,makecform('lab2srgb'));
    else  
       B = colorspace('RGB<-Lab',B);
    end

end

