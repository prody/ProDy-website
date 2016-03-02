display rendermode GLSL 
display projection orthographic
color Display Background white
display shadows on
display depthcue off
axes location off
stage location off
light 0 on
light 1 on
light 2 off
light 3 on
mol addrep 0
display resetview
mol new {./1gfl_3-7aa.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
mol modselect 0 0 protein
mol modstyle 0 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor 0 0 Structure
mol color Structure
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol selection protein
mol material Opaque
draw color red
draw line {-11.712  62.828 -14.936} {-10.489  63.603  -0.855} width 3 style solid 
draw line {-11.712  62.828 -14.936} {-12.027  66.84    0.313} width 3 style solid 
draw line {-11.712  62.828 -14.936} {  1.722  56.827 -24.305} width 3 style solid 
draw line {-11.712  62.828 -14.936} {  1.854  55.31  -20.831} width 3 style solid 
draw line {-11.712  62.828 -14.936} { 27.049  54.388  -5.042} width 3 style solid 
draw line {-11.712  62.828 -14.936} { 27.235  52.101  -8.143} width 3 style solid 
draw line {-11.712  62.828 -14.936} { 26.913  55.024 -10.577} width 3 style solid 
mol addrep 0
mol modselect 2 0 protein and name CA and resid 3 76 77 116 117 211 212 213
mol modcolor 2 0 ColorID 1
mol modstyle 2 0 VDW 0.600000 12.000000
mol color ColorID 1
mol representation VDW 1.000000 12.000000 
mol selection protein and name CA and resid 3 76 77 116 117 211 212 213
mol material Opaque 
draw color gray
draw color gray
draw color gray
draw line { -7.99   57.947 -13.219} {-10.489  63.603  -0.855} width 3 style solid 
draw line { -7.99   57.947 -13.219} {-12.027  66.84    0.313} width 3 style solid 
draw line { -7.99   57.947 -13.219} {  1.722  56.827 -24.305} width 3 style solid 
draw line { -7.99   57.947 -13.219} { 27.049  54.388  -5.042} width 3 style solid 
draw line { -7.99   57.947 -13.219} { 27.235  52.101  -8.143} width 3 style solid 
mol addrep 0
mol modselect 3 0 protein and name CA and resid 6 76 77 116 211 212
mol modcolor 3 0 ColorID 2
mol modstyle 3 0 VDW 0.600000 12.000000
mol color ColorID 2
mol representation VDW 1.000000 12.000000 
mol selection protein and name CA and resid 6 76 77 116 211 212
mol material Opaque 
draw color orange
mol addrep 0
