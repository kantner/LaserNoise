function cmap = temperatureColormap(m)

%from: http://soliton.vm.bytemark.co.uk/pub/cpt-city/arendal/temperature.cpt

x = linspace(0,1,18); 

c= [  30  92 179
	  23 111 193
	  11 142 216
	   4 161 230
	  25 181 241
	  51 188 207
	 102 204 206
	 153 219 184
	 192 229 136
	 204 230  75
	 243 240  29
	 254 222  39
	 252 199   7
	 248 157  14
	 245 114  21
	 241  71  28
	 219  30  38
	 164  38  44];
 
 

y = linspace(0,1,m);

cmap = interp1(x,c,y,'pchip')/255;
