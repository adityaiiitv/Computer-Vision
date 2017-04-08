% SFS
clear all;
close all;

function val = Rp(p,q,ps,qs)
    val = (ps*(q*q+1) - p*(qs*q+1)) / (sqrt((qs*qs+ps*ps+1) * (p*p+q*q+1)) * (p*p+q*q+1));
end

function val = Rq(p,q,ps,qs)
    val = (qs*(p^2+1) - q*(ps*p+1)) / (sqrt((qs^2+ps^2+1) * (p^2+q^2+1)) * (p^2+q^2+1));
end

function val = Rv(p,q,s)
    val = (s(1)*p + s(2)*q+1) / sqrt((s(2)*s(2)+s(1)*s(1)+1) * (p*p+q*q+1));
end

M=256;
N=256;

r=100;
% RADIUS OF THE SPHERE
Depth=zeros(M,N);
E=0.2*ones(M,N);
% BACKGROUND LIGHT
s_init=[0,0];
% p,q COORDINATE OF THE SOURCE

% INITIALIZATION OF DIFFERENT MATRICES
p_init=zeros(M,N);
q_init=zeros(M,N);
f_init=zeros(M,N);
g_init=zeros(M,N);
p_orig=zeros(M,N);
q_orig=zeros(M,N);
f_orig=zeros(M,N);

mask=zeros(M,N);
boundary=zeros(M,N);
boundaryfg=zeros(M,N);
g_orig=zeros(M,N);

for i=1:M,
    for j=1:N,
        curr_r = sqrt((i-M/2)^2 + (j-N/2)^2);
        if(curr_r<r) 
			% WE ASSIGN DEPTH ONLY TO THE POINTS PART OF THE SPHERE
            Depth(i,j) = round(sqrt(r^2 - (i-M/2)^2 - (j-N/2)^2));
            % ESTIMATE DEPTH FOR EACH PIXEL
            p = (i-M/2)/Depth(i,j);
            q = (j-N/2)/Depth(i,j);
            mask(i,j)=1;
            % CALCULATE E FROM p,q & s
            t = Rv(p, q, s_init);
            % CLOSED CURVE VALUES FOR p,q
            if(round(curr_r)==round(r/2))
                p_init(i,j)=p;
                q_init(i,j)=q;
                boundary(i,j)=1;
            end
            p_orig(i,j)=p;
            q_orig(i,j)=q;
            % CLOSED CURVE VALUES FOR f,g
            if(round(curr_r)>=r -1 && round(curr_r)<r+1)
                f_init(i,j)=2*p/(1+sqrt(1+p^2+q^2));
                g_init(i,j)=2*q/(1+sqrt(1+p^2+q^2));
                boundaryfg(i,j)=1;
            end
            f_orig(i,j)=2*p/(1+sqrt(1+p^2+q^2));
            g_orig(i,j)=2*q/(1+sqrt(1+p^2+q^2));
            
            % ENSURE THAT E IS NEGATIVE
            if(t>0)
                E(i,j)=t;
            else
                E(i,j)=0;
            end
            
        end
    end
end

% THE SOURCE IS PASSED TO ALGORITHM FOR RECONSTRUCTION
s=s_init;
figure;
imshow(E);

% NOISY IMAGE WITH GAUSSIAN NOISE
E_noise=imnoise(E,'gaussian',0,5);
save('ASSIGN2.mat','E','s','r','mask','boundary','p_init','q_init','f_init','g_init','E_noise','Depth','boundaryfg');
