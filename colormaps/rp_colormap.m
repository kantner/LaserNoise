function cmap = rp_colormap(m)

%from: http://soliton.vm.bytemark.co.uk/pub/cpt-city/arendal/temperature.cpt

x = linspace(0,1,314); 

c= [  
   152    17    52
   153    19    53
   154    22    54
   155    24    55
   157    27    57
   159    30    58
   159    32    59
   161    35    61
   162    37    62
   163    40    63
   165    43    65
   166    45    66
   168    48    67
   168    50    68
   170    53    70
   172    56    71
   172    58    72
   174    61    74
   175    63    75
   176    66    76
   178    69    78
   179    71    79
   181    74    80
   181    76    81
   183    79    83
   185    82    84
   185    84    85
   187    87    87
   188    89    88
   190    92    89
   191    95    91
   192    97    92
   194   100    93
   194   102    94
   196   105    96
   198   108    97
   199   110    98
   200   113   100
   201   115   100
   203   118   102
   204   121   104
   205   123   105
   207   126   106
   207   128   107
   209   131   109
   211   135   110
   212   136   111
   213   139   113
   214   141   113
   216   144   115
   217   148   117
   218   149   117
   220   152   119
   221   154   120
   222   157   121
   224   161   123
   225   162   124
   226   166   126
   227   167   126
   229   170   128
   230   174   130
   231   175   130
   233   179   132
   234   180   133
   235   183   134
   237   187   136
   238   188   137
   239   192   138
   240   193   139
   242   197   141
   243   200   142
   244   201   143
   246   205   145
   247   206   146
   248   210   147
   250   213   149
   251   214   150
   252   218   151
   253   219   152
   252   221   152
   251   220   152
   247   219   151
   244   218   150
   242   218   149
   239   217   148
   238   217   147
   234   216   146
   231   215   145
   229   214   145
   226   214   144
   224   213   143
   221   212   142
   218   211   141
   216   211   140
   213   210   139
   211   210   139
   208   209   138
   205   208   136
   203   207   136
   200   206   135
   198   206   134
   195   205   133
   191   204   132
   190   204   131
   186   203   130
   185   203   130
   181   202   129
   178   201   128
   177   200   127
   173   199   126
   172   199   125
   168   198   124
   165   197   123
   163   197   123
   160   196   122
   158   195   121
   155   195   120
   152   194   119
   150   193   118
   147   192   117
   145   192   117
   142   191   115
   139   190   114
   137   190   114
   134   189   113
   132   188   112
   129   188   111
   125   187   110
   124   186   109
   121   185   108
   119   185   108
   116   184   107
   112   183   106
   111   183   105
   107   182   104
   106   181   103
   102   181   102
    99   180   101
    97   179   101
    94   178    99
    92   178    99
    89   177    98
    86   176    97
    84   176    96
    81   175    95
    79   174    94
    76   173    93
    73   173    92
    71   172    92
    68   171    91
    66   171    90
    63   170    89
    60   169    88
    58   169    87
    55   168    86
    53   167    86
    50   166    85
    48   166    84
    49   166    86
    50   166    88
    50   165    89
    51   165    91
    51   165    92
    52   165    94
    53   165    96
    54   165    96
    54   164    98
    55   164    99
    56   164   101
    57   164   103
    57   164   104
    58   164   106
    58   164   107
    59   163   109
    60   163   111
    61   163   112
    61   163   114
    62   163   115
    63   162   117
    64   162   119
    64   162   120
    65   162   121
    65   162   122
    66   162   124
    67   161   126
    67   161   127
    68   161   129
    69   161   130
    70   161   132
    70   161   134
    71   160   135
    72   160   137
    72   160   138
    73   160   140
    74   160   142
    74   160   143
    75   159   144
    76   159   145
    77   159   147
    77   159   149
    78   159   150
    79   159   152
    79   159   153
    80   158   155
    81   158   157
    81   158   158
    82   158   160
    83   158   161
    83   157   163
    84   157   165
    85   157   166
    86   157   168
    86   157   168
    87   157   170
    88   156   172
    88   156   173
    89   156   175
    89   156   176
    90   156   178
    91   156   180
    92   155   181
    92   155   183
    93   155   184
    94   155   186
    95   155   188
    95   155   189
    96   154   191
    96   154   192
    97   154   193
    98   154   195
    99   154   196
    99   154   198
   100   154   199
   101   153   201
   102   153   203
   102   153   204
   104   155   205
   106   155   205
   108   157   206
   109   158   206
   112   160   207
   114   161   208
   115   162   208
   118   164   209
   119   164   210
   122   166   211
   124   168   211
   125   169   212
   128   170   213
   129   171   213
   131   173   214
   134   174   215
   135   175   215
   137   177   216
   139   177   216
   141   179   217
   144   181   218
   145   182   218
   147   183   219
   149   184   220
   151   186   220
   153   187   221
   155   188   222
   157   190   222
   158   191   223
   161   192   224
   163   194   224
   164   195   225
   167   196   226
   168   197   226
   171   199   227
   173   200   228
   174   201   228
   177   203   229
   178   204   229
   180   205   230
   183   207   231
   184   208   231
   186   209   232
   188   210   233
   190   212   233
   193   213   234
   194   214   235
   196   216   235
   197   217   236
   200   218   237
   202   220   237
   204   221   238
   206   222   239
   207   223   239
   210   225   240
   212   226   241
   213   227   241
   216   229   242
   217   230   242
   220   231   243
   222   233   244
   223   234   244
   226   235   245
   227   236   246
   229   238   246
   232   239   247
   233   240   248
   235   242   248
   237   243   249
   239   244   250
   242   246   251
   243   247   251
   245   248   252
   246   249   252
   249   251   253
   251   253   254
   253   253   254
   255   255   255];
 
 

y = linspace(0,1,m);

cmap = flipud(interp1(x,c,y,'pchip'))/255;
