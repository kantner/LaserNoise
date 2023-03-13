function cmap = divergingColormap1(N)



rgb = [...
    94    79   162
    50   136   189
   102   194   165
   171   221   164
   230   245   152
   255   255   191
   254   224   139
   253   174    97
   244   109    67
   213    62    79
   158     1    66  ] / 255;

cmap = interp1(linspace(0,1,11),rgb,linspace(0,1,N));