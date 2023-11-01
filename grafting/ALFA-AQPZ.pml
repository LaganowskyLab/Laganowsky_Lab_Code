reinitialize

# aqpz tetramer from PDB 1RC2
load 1RC2_edited.pdb 
create target, 1RC2_edited
delete 1RC2_edited

# ALFA-nB complex PDB 6I2G
load alfanb.pdb 
create nb, alfanb and chain A+B and polymer
delete alfanb

# UNCOMMENT 'fab' to try different grafts
# try c-term grafting
# fusion 0
#   -----------------**************
#fab GGIIGGLIYRTLLEKRDASRLEEELRRRLTE, hf, ss=1
# -1
fab GGIIGGLIYRTLLKRDASRLEEELRRRLTE, hf, ss=1
# -2
#fab GGIIGGLIYRTLLKRASRLEEELRRRLTE, hf, ss=1
# -3
#fab GGIIGGLIYRTLLRASRLEEELRRRLTE, hf, ss=1

align hf, target
align nb and chain B, hf

create hf2, hf
align hf2, target and chain B

create nb2, nb
align nb2 and chain B, hf2

create hf3, hf
align hf3, target and chain C

create nb3, nb
align nb3 and chain B, hf3

create hf4, hf
align hf4, target and chain D

create nb4, nb
align nb4 and chain B, hf4

color grey50, target
color chartreuse, nb
color tv_red, hf*
color violet, target and chain A
color tv_blue, target and chain B
color tv_orangee, target and chain C
color tv_yellow, target and chain D

orient
hide everything, not polymer
as cartoon, polymer
#show sticks

# bottom view
### cut below here and paste into script ###
set_view (\
     0.966829896,    0.255419463,   -0.000000068,\
    -0.255419463,    0.966829896,    0.000000076,\
     0.000000085,   -0.000000056,    1.000000000,\
     0.000000043,   -0.000000028, -307.790313721,\
   -46.775100708,   46.774955750,   21.308921814,\
   260.736022949,  343.944396973,  -20.000000000 )
### cut above here and paste into script ###

turn x, 90

