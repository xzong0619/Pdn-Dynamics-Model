#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White}
camera {orthographic
  right -52.83*x up 55.29*y
  direction 1.00*z
  location <0,0,50.00> look_at <0,0,0>}
light_source {<  2.00,   3.00,  40.00> color White
  area_light <0.70, 0, 0>, <0, 0.70, 0>, 3, 3
  adaptive 1 jitter}

#declare simple = finish {phong 0.7}
#declare pale = finish {ambient .5 diffuse .85 roughness .001 specular 0.200 }
#declare intermediate = finish {ambient 0.3 diffuse 0.6 specular 0.10 roughness 0.04 }
#declare vmd = finish {ambient .0 diffuse .65 phong 0.1 phong_size 40. specular 0.500 }
#declare jmol = finish {ambient .2 diffuse .6 specular 1 roughness .001 metallic}
#declare ase2 = finish {ambient 0.05 brilliance 3 diffuse 0.6 metallic specular 0.70 roughness 0.04 reflection 0.15}
#declare ase3 = finish {ambient .15 brilliance 2 diffuse .6 metallic specular 1. roughness .001 reflection .0}
#declare glass = finish {ambient .05 diffuse .3 specular 1. roughness .001}
#declare glass2 = finish {ambient .0 diffuse .3 specular 1. reflection .25 roughness .001}
#declare Rcell = 0.070;
#declare Rbond = 0.100;

#macro atom(LOC, R, COL, TRANS, FIN)
  sphere{LOC, R texture{pigment{color COL transmit TRANS} finish{FIN}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{torus{R, Rcell rotate 45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      torus{R, Rcell rotate -45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      translate LOC}
#end

cylinder {<-11.59,   5.51, -39.66>, <  4.76, -17.55, -26.79>, Rcell pigment {Black}}
cylinder {<  3.37,   7.65, -36.09>, < 19.72, -15.40, -23.22>, Rcell pigment {Black}}
cylinder {< -3.35,  15.86, -12.85>, < 13.00,  -7.20,   0.02>, Rcell pigment {Black}}
cylinder {<-18.31,  13.71, -16.42>, < -1.96,  -9.34,  -3.55>, Rcell pigment {Black}}
cylinder {<-11.59,   5.51, -39.66>, <  3.37,   7.65, -36.09>, Rcell pigment {Black}}
cylinder {<  4.76, -17.55, -26.79>, < 19.72, -15.40, -23.22>, Rcell pigment {Black}}
cylinder {< -1.96,  -9.34,  -3.55>, < 13.00,  -7.20,   0.02>, Rcell pigment {Black}}
cylinder {<-18.31,  13.71, -16.42>, < -3.35,  15.86, -12.85>, Rcell pigment {Black}}
cylinder {<-11.59,   5.51, -39.66>, <-18.31,  13.71, -16.42>, Rcell pigment {Black}}
cylinder {<  4.76, -17.55, -26.79>, < -1.96,  -9.34,  -3.55>, Rcell pigment {Black}}
cylinder {< 19.72, -15.40, -23.22>, < 13.00,  -7.20,   0.02>, Rcell pigment {Black}}
cylinder {<  3.37,   7.65, -36.09>, < -3.35,  15.86, -12.85>, Rcell pigment {Black}}
atom(<-14.22,   8.72, -30.57>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #0 
atom(<-10.48,   9.26, -29.67>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #1 
atom(<-12.17,   5.84, -28.96>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #2 
atom(< -8.43,   6.37, -28.06>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #3 
atom(<-10.45,   5.31, -27.40>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #4 
atom(<-11.10,  10.02, -27.51>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #5 
atom(< -6.71,   5.85, -26.51>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #6 
atom(<-14.84,   9.48, -28.40>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #7 
atom(<-12.50,   8.19, -29.01>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #8 
atom(< -9.06,   7.14, -25.90>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #9 
atom(< -8.76,   8.73, -28.12>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #10 
atom(<-12.80,   6.60, -26.79>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #11 
atom(< -7.34,   6.61, -24.35>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #12 
atom(<-11.08,   6.07, -25.24>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #13 
atom(< -9.38,   9.49, -25.96>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #14 
atom(<-13.12,   8.96, -26.85>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #15 
atom(< -7.66,   8.96, -24.40>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #16 
atom(<-11.40,   8.43, -25.29>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #17 
atom(< -5.62,   6.08, -22.79>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #18 
atom(< -9.36,   5.55, -23.69>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #19 
atom(< -6.74,   9.79, -28.78>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #20 
atom(< -3.00,  10.33, -27.89>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #21 
atom(< -4.69,   6.91, -27.17>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #22 
atom(< -0.95,   7.45, -26.28>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #23 
atom(< -2.97,   6.38, -25.62>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #24 
atom(< -3.62,  11.09, -25.73>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #25 
atom(<  0.77,   6.92, -24.72>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #26 
atom(< -7.36,  10.56, -26.62>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #27 
atom(< -5.02,   9.26, -27.23>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #28 
atom(< -1.58,   8.21, -24.12>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #29 
atom(< -1.28,   9.80, -26.33>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #30 
atom(< -5.32,   7.67, -25.01>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #31 
atom(<  0.14,   7.68, -22.56>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #32 
atom(< -3.60,   7.15, -23.45>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #33 
atom(< -1.90,  10.56, -24.17>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #34 
atom(< -5.64,  10.03, -25.06>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #35 
atom(< -0.18,  10.04, -22.62>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #36 
atom(< -3.92,   9.50, -23.51>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #37 
atom(<  1.86,   7.15, -21.01>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #38 
atom(< -1.88,   6.62, -21.90>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #39 
atom(<-10.13,   2.96, -27.35>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #40 
atom(< -6.39,   3.49, -26.46>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #41 
atom(< -8.09,   0.07, -25.74>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #42 
atom(< -4.35,   0.61, -24.85>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #43 
atom(< -6.37,  -0.45, -24.18>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #44 
atom(< -7.02,   4.26, -24.29>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #45 
atom(< -2.63,   0.08, -23.29>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #46 
atom(<-10.76,   3.72, -25.19>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #47 
atom(< -8.41,   2.43, -25.79>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #48 
atom(< -4.97,   1.37, -22.68>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #49 
atom(< -4.67,   2.96, -24.90>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #50 
atom(< -8.71,   0.84, -23.58>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #51 
atom(< -3.25,   0.85, -21.13>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #52 
atom(< -6.99,   0.31, -22.02>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #53 
atom(< -5.30,   3.73, -22.74>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #54 
atom(< -9.04,   3.19, -23.63>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #55 
atom(< -3.58,   3.20, -21.18>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #56 
atom(< -7.32,   2.66, -22.08>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #57 
atom(< -1.53,   0.32, -19.58>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #58 
atom(< -5.27,  -0.22, -20.47>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #59 
atom(< -2.65,   4.03, -25.56>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #60 
atom(<  1.09,   4.56, -24.67>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #61 
atom(< -0.61,   1.15, -23.95>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #62 
atom(<  3.13,   1.68, -23.06>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #63 
atom(<  1.11,   0.62, -22.40>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #64 
atom(<  0.46,   5.33, -22.51>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #65 
atom(<  4.85,   1.15, -21.51>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #66 
atom(< -3.28,   4.79, -23.40>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #67 
atom(< -0.93,   3.50, -24.01>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #68 
atom(<  2.51,   2.45, -20.90>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #69 
atom(<  2.81,   4.04, -23.12>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #70 
atom(< -1.23,   1.91, -21.79>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #71 
atom(<  4.23,   1.92, -19.34>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #72 
atom(<  0.49,   1.38, -20.24>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #73 
atom(<  2.18,   4.80, -20.95>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #74 
atom(< -1.56,   4.26, -21.85>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #75 
atom(<  3.90,   4.27, -19.40>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #76 
atom(<  0.16,   3.74, -20.29>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #77 
atom(<  5.95,   1.39, -17.79>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #78 
atom(<  2.21,   0.85, -18.68>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #79 
atom(< -6.04,  -2.81, -24.13>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #80 
atom(< -2.30,  -2.27, -23.24>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #81 
atom(< -4.00,  -5.69, -22.52>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #82 
atom(< -0.26,  -5.15, -21.63>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #83 
atom(< -2.28,  -6.22, -20.97>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #84 
atom(< -2.93,  -1.51, -21.08>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #85 
atom(<  1.46,  -5.68, -20.07>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #86 
atom(< -6.67,  -2.04, -21.97>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #87 
atom(< -4.32,  -3.34, -22.58>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #88 
atom(< -0.88,  -4.39, -19.47>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #89 
atom(< -0.58,  -2.80, -21.68>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #90 
atom(< -4.62,  -4.93, -20.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #91 
atom(<  1.44,  -1.74, -22.35>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #92 
atom(<  5.18,  -1.20, -21.45>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #93 
atom(<  3.48,  -4.62, -20.74>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #94 
atom(<  7.22,  -4.08, -19.84>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #95 
atom(<  5.20,  -5.15, -19.18>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #96 
atom(<  4.55,  -0.44, -19.29>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #97 
atom(<  8.94,  -4.61, -18.29>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #98 
atom(<  0.81,  -0.97, -20.18>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #99 
atom(<  3.16,  -2.26, -20.79>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #100 
atom(<  6.60,  -3.32, -17.68>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #101 
atom(<  6.90,  -1.73, -19.90>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #102 
atom(<  2.86,  -3.85, -18.57>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #103 
atom(< -1.96,  -8.57, -20.91>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #104 
atom(<  1.78,  -8.04, -20.02>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #105 
atom(<  0.09, -11.45, -19.30>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #106 
atom(<  3.83, -10.92, -18.41>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #107 
atom(<  1.81, -11.98, -17.75>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #108 
atom(<  1.16,  -7.27, -17.86>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #109 
atom(<  5.55, -11.45, -16.86>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #110 
atom(< -2.58,  -7.81, -18.75>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #111 
atom(< -0.24,  -9.10, -19.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #112 
atom(<  3.20, -10.15, -16.25>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #113 
atom(<  3.50,  -8.56, -18.47>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #114 
atom(< -0.54, -10.69, -17.14>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #115 
atom(<  5.52,  -7.50, -19.13>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #116 
atom(<  9.26,  -6.96, -18.24>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #117 
atom(<  7.57, -10.38, -17.52>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #118 
atom(< 11.31,  -9.85, -16.63>, 1.82, rgb <1.00, 1.00, 0.78>, 0.0, ase2) // #119 
atom(<  9.29, -10.91, -15.96>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #120 
atom(<  8.64,  -6.20, -16.07>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #121 
atom(< 13.03, -10.37, -15.07>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #122 
atom(<  4.90,  -6.74, -16.97>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #123 
atom(<  7.24,  -8.03, -17.57>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #124 
atom(< 10.68,  -9.08, -14.46>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #125 
atom(< 10.98,  -7.49, -16.68>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #126 
atom(<  6.94,  -9.62, -15.36>, 0.59, rgb <1.00, 0.05, 0.05>, 0.0, ase2) // #127 
