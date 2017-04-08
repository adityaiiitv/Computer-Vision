% Shape from Shading
% Image generator with source direction = s_orig
close all;
clear all;	

function val = Ravg(A,i,j)
    val = ( A(i,j-1) + A(i-1,j) + A(i,j+1) + A(i+1,j) )/6 +(A(i-1,j-1)+ A(i-1,j+1)+A(i+1,j-1)+A(i+1,j+1))/12;
end

function val = Rp(p,q,ps,qs)
    val = (ps*(q*q +1) - p*(qs*q + 1))/(sqrt((qs*qs + ps*ps +1) * (p*p + q*q +1)) * (p*p+q*q+1));
end


function val = Rq(p,q,ps,qs)
    val = (qs*(p^2 +1) - q*(ps*p + 1))/(sqrt((qs^2 + ps^2 +1) * (p^2 + q^2 +1)) * (p^2+q^2+1));
end

function val = Rval(p,q,s)
    val = (s(1)*p + s(2)*q + 1)/sqrt( (s(2)*s(2) + s(1)*s(1) +1) * (p*p + q*q +1));
end


%This is the size of the image(MxN)
M=256;
N=256;

% Radius of sphere
radius=100;

%	Initializing depth of the image
Depth = zeros(M,N);	 
p_map =	zeros(M,N);
q_map =	zeros(M,N);
% Background Light 
E = 0.0 * ones(M,N); 

% p,q co-ordinate of source
%s_orig=[0,1]; 
%s_orig=[0,0];
s_orig=[0.5,0.5];

for i=1:M,
    for j=1:N,
        current_radius = sqrt((i-M/2)^2 + (j-N/2)^2);
        % Assign Depth only for points that are part of sphere in image
        if(current_radius < radius) 		            
         	% estimate depth for each pixel of sphere
         	Depth(i,j) = round(sqrt(radius^2 - (i-M/2)^2 - (j-N/2)^2));

        	%Surface gradients of the image.
           	p = (i-M/2)/Depth(i,j);
            p_map(i,j)=p;
            %Rate of change of depth in x and y direction
            q = (j-N/2)/Depth(i,j);
      		q_map(i,j)=q;      
           	% Calculate E from p,q & s
            temp = Rval(p, q, s_orig);
            
            % Ensure nonegativity of E
            if(temp>0)
                E(i,j) = temp;
            
            else
                E(i,j) = 0;

            end
            
        end
    end
end

figure;
imshow(mat2gray(E));

projection=E;
for i=1:M,
    for j=1:N,
        current_radius = sqrt((i-M/2)^2 + (j-N/2)^2);
        if(current_radius < radius) 	% Assign Depth only for points that are part of sphere in image
        projection(i,j)=1;
    	endif
    endfor
endfor

figure;
imshow(mat2gray(projection));

