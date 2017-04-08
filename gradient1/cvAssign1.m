R = 255;   
for i = 1:512
	for j=1:512
		% THE CENTER OF THE IMAGE IS AT 256, 256:: WE MAKE IT 0, 0
		X = 256-i;
   		Y = 256-j;

  		if(X^2 + Y^2 < R^2)     
  			c(i,j) = sqrt(R^2-X^2 - Y^2)*100;
    
    			% FOR PS=0 AND QS=0
    			r1(i,j) = sqrt(1-((X^2 +Y^2)/R^2));
    
    			% FOR PS=0.5 AND QS=0.5
   			r2(i,j) = ((sqrt(R^2 -(X^2 +Y^2)))- 0.5*(X + Y))/(sqrt(1.5)*R);
   
   			% FOR PS=0 AND QS=1
   			r3(i,j) = (sqrt(R^2 - (X^2 +Y^2)) - Y)/(sqrt(2)*R); 
 		endif
 	endfor 
endfor
imshow(r2)

