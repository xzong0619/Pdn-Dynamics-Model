#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White}
camera {orthographic
  right -17.00*x up 15.00*y
  direction 1.00*z
  location <0,0,50.00> look_at <0,0,0>}
light_source {<  2.00,   3.00,  125.00> color White
  area_light <0.95, 0, 0>, <0, 0.80, 0>, 5, 4
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
#declare Rcell = 0.100;
#declare Rbond = 0.100;

#macro atom(LOC, R, COL, TRANS, FIN)
  sphere{LOC, R texture{pigment{color COL transmit TRANS} finish{FIN}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{torus{R, Rcell rotate 45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      torus{R, Rcell rotate -45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      translate LOC}
#end

atom(<-14.50, -19.53,  -5.49>, 1.20, rgb <0.87, 0.87, 0.90>, 0.0, jmol) // #0 
atom(<-23.81, -11.44, -11.28>, 1.20, rgb <0.87, 0.87, 0.90>, 0.0, jmol) // #1 
atom(<-18.85, -15.45,  -5.49>, 1.20, rgb <0.87, 0.87, 0.90>, 0.0, jmol) // #2 
atom(<-19.47, -15.51, -11.28>, 1.20, rgb <0.87, 0.87, 0.90>, 0.0, jmol) // #3 
atom(<-16.35, -18.74,  -4.77>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #4 
atom(<-21.96, -12.22, -10.55>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #5 
atom(<-20.69, -16.24,  -4.77>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #6 
atom(<-17.62, -14.73, -10.55>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #7 
atom(<-17.87, -17.11,  -4.70>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #8 
atom(<-20.44, -13.85, -10.49>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #9 
atom(<-22.22, -17.87,  -4.70>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #10 
atom(<-16.10, -13.10, -10.49>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #11 
atom(<-15.42, -13.18,  -6.29>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #12 
atom(<-22.90, -17.78,  -0.50>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #13 
atom(<-19.76, -13.76,  -6.29>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #14 
atom(<-18.56, -17.20,  -0.50>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #15 
atom(<-17.02, -14.75,  -6.27>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #16 
atom(<-21.30, -16.22,  -0.48>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #17 
atom(<-21.36, -12.20,  -6.27>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #18 
atom(<-16.96, -18.76,  -0.48>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #19 
atom(<-23.04, -12.53,  -3.66>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #20 
atom(<-15.27, -18.43,  -9.44>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #21 
atom(<-18.70, -14.41,  -3.66>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #22 
atom(<-19.61, -16.55,  -9.44>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #23 
atom(<-23.50, -12.48,  -1.47>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #24 
atom(<-14.82, -18.48,  -7.26>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #25 
atom(<-19.16, -14.47,  -1.47>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #26 
atom(<-19.16, -16.50,  -7.26>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #27 
atom(<-18.07, -19.75,  -2.85>, 0.65, rgb <0.19, 0.31, 0.97>, 0.0, jmol) // #28 
atom(<-20.25, -11.21,  -8.64>, 0.65, rgb <0.19, 0.31, 0.97>, 0.0, jmol) // #29 
atom(<-13.73, -15.23,  -2.85>, 0.65, rgb <0.19, 0.31, 0.97>, 0.0, jmol) // #30 
atom(<-24.59, -15.73,  -8.64>, 0.65, rgb <0.19, 0.31, 0.97>, 0.0, jmol) // #31 
atom(<-18.16, -20.79,  -2.94>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #32 
atom(<-20.15, -10.18,  -8.73>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #33 
atom(<-13.82, -14.19,  -2.94>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #34 
atom(<-24.49, -16.77,  -8.73>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #35 
atom(<-17.41, -19.39,  -3.60>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #36 
atom(<-20.91, -11.58,  -9.39>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #37 
atom(<-13.07, -15.60,  -3.60>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #38 
atom(<-25.25, -15.37,  -9.39>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #39 
atom(<-17.66, -19.47,  -1.92>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #40 
atom(<-20.65, -11.50,  -7.70>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #41 
atom(<-13.32, -15.51,  -1.92>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #42 
atom(<-24.99, -15.45,  -7.70>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #43 
atom(<-16.26, -17.10,  -5.97>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #44 
atom(<-22.06, -13.86,  -0.18>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #45 
atom(<-20.60, -17.88,  -5.97>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #46 
atom(<-17.72, -13.08,  -0.18>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #47 
atom(<-16.90, -13.27,  -4.87>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #48 
atom(<-21.42, -17.70, -10.66>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #49 
atom(<-21.24, -13.68,  -4.87>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #50 
atom(<-17.08, -17.29, -10.66>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #51 
atom(<-14.95, -18.87,  -2.68>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #52 
atom(<-23.36, -12.09,  -8.47>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #53 
atom(<-19.30, -16.11,  -2.68>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #54 
atom(<-19.02, -14.85,  -8.47>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #55 
atom(<-19.81, -19.46,  -3.99>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #56 
atom(<-18.51, -11.50,  -9.77>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #57 
atom(<-15.47, -15.52,  -3.99>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #58 
atom(<-22.85, -15.44,  -9.77>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #59 
atom(<-19.25, -18.05,  -3.07>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #60 
atom(<-19.06, -12.91,  -8.85>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #61 
atom(<-14.91, -16.93,  -3.07>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #62 
atom(<-23.40, -14.03,  -8.85>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #63 
atom(<-19.88, -19.22,  -0.91>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #64 
atom(<-18.44, -11.74,  -6.70>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #65 
atom(<-15.54, -15.76,  -0.91>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #66 
atom(<-22.78, -15.20,  -6.70>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #67 
atom(<-21.27, -18.97,  -1.98>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #68 
atom(<-17.04, -11.99,  -7.77>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #69 
atom(<-16.93, -16.01,  -1.98>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #70 
atom(<-21.39, -14.95,  -7.77>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #71 
atom(<-20.53, -20.59,  -1.86>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #72 
atom(<-17.79, -10.37,  -7.64>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #73 
atom(<-16.18, -14.39,  -1.86>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #74 
atom(<-22.13, -16.58,  -7.64>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #75 
atom(<-16.82, -17.63,  -5.17>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #76 
atom(<-21.49, -13.34, -10.96>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #77 
atom(<-21.16, -17.35,  -5.17>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #78 
atom(<-17.15, -13.61, -10.96>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #79 
atom(<-16.45, -13.71,  -5.79>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #80 
atom(<-21.87, -17.25,   0.00>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #81 
atom(<-20.79, -13.23,  -5.79>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #82 
atom(<-17.53, -17.73,   0.00>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #83 
atom(<-23.40, -11.92,  -2.60>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #84 
atom(<-14.92, -19.05,  -8.39>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #85 
atom(<-19.06, -15.03,  -2.60>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #86 
atom(<-19.26, -15.94,  -8.39>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #87 
atom(<-19.41, -19.14,  -3.01>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #88 
atom(<-18.91, -11.83,  -8.80>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #89 
atom(<-15.07, -15.84,  -3.01>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #90 
atom(<-23.25, -15.12,  -8.80>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #91 
atom(<-20.32, -19.51,  -1.88>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #92 
atom(<-17.99, -11.45,  -7.67>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #93 
atom(<-15.98, -15.47,  -1.88>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #94 
atom(<-22.33, -15.50,  -7.67>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #95 
atom(<-14.50, -11.49,  -5.49>, 1.20, rgb <0.87, 0.87, 0.90>, 0.0, jmol) // #96 
atom(<-23.81,  -3.40, -11.28>, 1.20, rgb <0.87, 0.87, 0.90>, 0.0, jmol) // #97 
atom(<-18.85,  -7.42,  -5.49>, 1.20, rgb <0.87, 0.87, 0.90>, 0.0, jmol) // #98 
atom(<-19.47,  -7.47, -11.28>, 1.20, rgb <0.87, 0.87, 0.90>, 0.0, jmol) // #99 
atom(<-16.35, -10.71,  -4.77>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #100 
atom(<-21.96,  -4.18, -10.55>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #101 
atom(<-20.69,  -8.20,  -4.77>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #102 
atom(<-17.62,  -6.69, -10.55>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #103 
atom(<-17.87,  -9.08,  -4.70>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #104 
atom(<-20.44,  -5.82, -10.49>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #105 
atom(<-22.22,  -9.83,  -4.70>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #106 
atom(<-16.10,  -5.06, -10.49>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #107 
atom(<-15.42,  -5.15,  -6.29>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #108 
atom(<-22.90,  -9.74,  -0.50>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #109 
atom(<-19.76,  -5.73,  -6.29>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #110 
atom(<-18.56,  -9.17,  -0.50>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #111 
atom(<-17.02,  -6.71,  -6.27>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #112 
atom(<-21.30,  -8.18,  -0.48>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #113 
atom(<-21.36,  -4.16,  -6.27>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #114 
atom(<-16.96, -10.73,  -0.48>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #115 
atom(<-23.04,  -4.50,  -3.66>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #116 
atom(<-15.27, -10.40,  -9.44>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #117 
atom(<-18.70,  -6.38,  -3.66>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #118 
atom(<-19.61,  -8.51,  -9.44>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #119 
atom(<-23.50,  -4.45,  -1.47>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #120 
atom(<-14.82, -10.45,  -7.26>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #121 
atom(<-19.16,  -6.43,  -1.47>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #122 
atom(<-19.16,  -8.46,  -7.26>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #123 
atom(<-18.07, -11.71,  -2.85>, 0.65, rgb <0.19, 0.31, 0.97>, 0.0, jmol) // #124 
atom(<-20.25,  -3.18,  -8.64>, 0.65, rgb <0.19, 0.31, 0.97>, 0.0, jmol) // #125 
atom(<-13.73,  -7.20,  -2.85>, 0.65, rgb <0.19, 0.31, 0.97>, 0.0, jmol) // #126 
atom(<-24.59,  -7.70,  -8.64>, 0.65, rgb <0.19, 0.31, 0.97>, 0.0, jmol) // #127 
atom(<-18.16, -12.75,  -2.94>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #128 
atom(<-20.15,  -2.14,  -8.73>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #129 
atom(<-13.82,  -6.16,  -2.94>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #130 
atom(<-24.49,  -8.74,  -8.73>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #131 
atom(<-17.41, -11.35,  -3.60>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #132 
atom(<-20.91,  -3.54,  -9.39>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #133 
atom(<-13.07,  -7.56,  -3.60>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #134 
atom(<-25.25,  -7.33,  -9.39>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #135 
atom(<-17.66, -11.43,  -1.92>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #136 
atom(<-20.65,  -3.46,  -7.70>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #137 
atom(<-13.32,  -7.48,  -1.92>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #138 
atom(<-24.99,  -7.42,  -7.70>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #139 
atom(<-16.26,  -9.07,  -5.97>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #140 
atom(<-22.06,  -5.83,  -0.18>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #141 
atom(<-20.60,  -9.84,  -5.97>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #142 
atom(<-17.72,  -5.05,  -0.18>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #143 
atom(<-16.90,  -5.23,  -4.87>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #144 
atom(<-21.42,  -9.66, -10.66>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #145 
atom(<-21.24,  -5.64,  -4.87>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #146 
atom(<-17.08,  -9.25, -10.66>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #147 
atom(<-14.95, -10.84,  -2.68>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #148 
atom(<-23.36,  -4.06,  -8.47>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #149 
atom(<-19.30,  -8.07,  -2.68>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #150 
atom(<-19.02,  -6.82,  -8.47>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #151 
atom(<-19.81, -11.43,  -3.99>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #152 
atom(<-18.51,  -3.47,  -9.77>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #153 
atom(<-15.47,  -7.48,  -3.99>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #154 
atom(<-22.85,  -7.41,  -9.77>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #155 
atom(<-19.25, -10.02,  -3.07>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #156 
atom(<-19.06,  -4.88,  -8.85>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #157 
atom(<-14.91,  -8.89,  -3.07>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #158 
atom(<-23.40,  -6.00,  -8.85>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #159 
atom(<-19.88, -11.18,  -0.91>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #160 
atom(<-18.44,  -3.71,  -6.70>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #161 
atom(<-15.54,  -7.73,  -0.91>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #162 
atom(<-22.78,  -7.17,  -6.70>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #163 
atom(<-21.27, -10.94,  -1.98>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #164 
atom(<-17.04,  -3.96,  -7.77>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #165 
atom(<-16.93,  -7.98,  -1.98>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #166 
atom(<-21.39,  -6.92,  -7.77>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #167 
atom(<-20.53, -12.56,  -1.86>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #168 
atom(<-17.79,  -2.33,  -7.64>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #169 
atom(<-16.18,  -6.35,  -1.86>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #170 
atom(<-22.13,  -8.54,  -7.64>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #171 
atom(<-16.82,  -9.59,  -5.17>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #172 
atom(<-21.49,  -5.30, -10.96>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #173 
atom(<-21.16,  -9.32,  -5.17>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #174 
atom(<-17.15,  -5.57, -10.96>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #175 
atom(<-16.45,  -5.68,  -5.79>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #176 
atom(<-21.87,  -9.22,   0.00>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #177 
atom(<-20.79,  -5.20,  -5.79>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #178 
atom(<-17.53,  -9.69,   0.00>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #179 
atom(<-23.40,  -3.88,  -2.60>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #180 
atom(<-14.92, -11.01,  -8.39>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #181 
atom(<-19.06,  -6.99,  -2.60>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #182 
atom(<-19.26,  -7.90,  -8.39>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #183 
atom(<-19.41, -11.10,  -3.01>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #184 
atom(<-18.91,  -3.79,  -8.80>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #185 
atom(<-15.07,  -7.81,  -3.01>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #186 
atom(<-23.25,  -7.08,  -8.80>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #187 
atom(<-20.32, -11.48,  -1.88>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #188 
atom(<-17.99,  -3.41,  -7.67>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #189 
atom(<-15.98,  -7.43,  -1.88>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #190 
atom(<-22.33,  -7.46,  -7.67>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #191 
atom(< -5.82, -19.53,  -5.49>, 1.20, rgb <0.87, 0.87, 0.90>, 0.0, jmol) // #192 
atom(<-15.13, -11.44, -11.28>, 1.20, rgb <0.87, 0.87, 0.90>, 0.0, jmol) // #193 
atom(<-10.16, -15.45,  -5.49>, 1.20, rgb <0.87, 0.87, 0.90>, 0.0, jmol) // #194 
atom(<-10.79, -15.51, -11.28>, 1.20, rgb <0.87, 0.87, 0.90>, 0.0, jmol) // #195 
atom(< -7.67, -18.74,  -4.77>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #196 
atom(<-13.28, -12.22, -10.55>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #197 
atom(<-12.01, -16.24,  -4.77>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #198 
atom(< -8.94, -14.73, -10.55>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #199 
atom(< -9.19, -17.11,  -4.70>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #200 
atom(<-11.76, -13.85, -10.49>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #201 
atom(<-13.53, -17.87,  -4.70>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #202 
atom(< -7.42, -13.10, -10.49>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #203 
atom(< -6.74, -13.18,  -6.29>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #204 
atom(<-14.21, -17.78,  -0.50>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #205 
atom(<-11.08, -13.76,  -6.29>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #206 
atom(< -9.87, -17.20,  -0.50>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #207 
atom(< -8.33, -14.75,  -6.27>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #208 
atom(<-12.62, -16.22,  -0.48>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #209 
atom(<-12.67, -12.20,  -6.27>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #210 
atom(< -8.28, -18.76,  -0.48>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #211 
atom(<-14.36, -12.53,  -3.66>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #212 
atom(< -6.59, -18.43,  -9.44>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #213 
atom(<-10.02, -14.41,  -3.66>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #214 
atom(<-10.93, -16.55,  -9.44>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #215 
atom(<-14.81, -12.48,  -1.47>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #216 
atom(< -6.13, -18.48,  -7.26>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #217 
atom(<-10.47, -14.47,  -1.47>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #218 
atom(<-10.48, -16.50,  -7.26>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #219 
atom(< -9.39, -19.75,  -2.85>, 0.65, rgb <0.19, 0.31, 0.97>, 0.0, jmol) // #220 
atom(<-11.56, -11.21,  -8.64>, 0.65, rgb <0.19, 0.31, 0.97>, 0.0, jmol) // #221 
atom(< -5.04, -15.23,  -2.85>, 0.65, rgb <0.19, 0.31, 0.97>, 0.0, jmol) // #222 
atom(<-15.90, -15.73,  -8.64>, 0.65, rgb <0.19, 0.31, 0.97>, 0.0, jmol) // #223 
atom(< -9.48, -20.79,  -2.94>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #224 
atom(<-11.47, -10.18,  -8.73>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #225 
atom(< -5.14, -14.19,  -2.94>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #226 
atom(<-15.81, -16.77,  -8.73>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #227 
atom(< -8.73, -19.39,  -3.60>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #228 
atom(<-12.22, -11.58,  -9.39>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #229 
atom(< -4.38, -15.60,  -3.60>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #230 
atom(<-16.56, -15.37,  -9.39>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #231 
atom(< -8.98, -19.47,  -1.92>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #232 
atom(<-11.97, -11.50,  -7.70>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #233 
atom(< -4.64, -15.51,  -1.92>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #234 
atom(<-16.31, -15.45,  -7.70>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #235 
atom(< -7.57, -17.10,  -5.97>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #236 
atom(<-13.38, -13.86,  -0.18>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #237 
atom(<-11.91, -17.88,  -5.97>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #238 
atom(< -9.03, -13.08,  -0.18>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #239 
atom(< -8.21, -13.27,  -4.87>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #240 
atom(<-12.74, -17.70, -10.66>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #241 
atom(<-12.56, -13.68,  -4.87>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #242 
atom(< -8.39, -17.29, -10.66>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #243 
atom(< -6.27, -18.87,  -2.68>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #244 
atom(<-14.68, -12.09,  -8.47>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #245 
atom(<-10.61, -16.11,  -2.68>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #246 
atom(<-10.34, -14.85,  -8.47>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #247 
atom(<-11.12, -19.46,  -3.99>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #248 
atom(< -9.82, -11.50,  -9.77>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #249 
atom(< -6.78, -15.52,  -3.99>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #250 
atom(<-14.17, -15.44,  -9.77>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #251 
atom(<-10.57, -18.05,  -3.07>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #252 
atom(<-10.38, -12.91,  -8.85>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #253 
atom(< -6.23, -16.93,  -3.07>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #254 
atom(<-14.72, -14.03,  -8.85>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #255 
atom(<-11.20, -19.22,  -0.91>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #256 
atom(< -9.75, -11.74,  -6.70>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #257 
atom(< -6.85, -15.76,  -0.91>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #258 
atom(<-14.09, -15.20,  -6.70>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #259 
atom(<-12.59, -18.97,  -1.98>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #260 
atom(< -8.36, -11.99,  -7.77>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #261 
atom(< -8.25, -16.01,  -1.98>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #262 
atom(<-12.70, -14.95,  -7.77>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #263 
atom(<-11.84, -20.59,  -1.86>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #264 
atom(< -9.11, -10.37,  -7.64>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #265 
atom(< -7.50, -14.39,  -1.86>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #266 
atom(<-13.45, -16.58,  -7.64>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #267 
atom(< -8.14, -17.63,  -5.17>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #268 
atom(<-12.81, -13.34, -10.96>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #269 
atom(<-12.48, -17.35,  -5.17>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #270 
atom(< -8.47, -13.61, -10.96>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #271 
atom(< -7.76, -13.71,  -5.79>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #272 
atom(<-13.19, -17.25,   0.00>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #273 
atom(<-12.10, -13.23,  -5.79>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #274 
atom(< -8.84, -17.73,   0.00>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #275 
atom(<-14.72, -11.92,  -2.60>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #276 
atom(< -6.23, -19.05,  -8.39>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #277 
atom(<-10.37, -15.03,  -2.60>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #278 
atom(<-10.57, -15.94,  -8.39>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #279 
atom(<-10.73, -19.14,  -3.01>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #280 
atom(<-10.22, -11.83,  -8.80>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #281 
atom(< -6.38, -15.84,  -3.01>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #282 
atom(<-14.57, -15.12,  -8.80>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #283 
atom(<-11.64, -19.51,  -1.88>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #284 
atom(< -9.31, -11.45,  -7.67>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #285 
atom(< -7.30, -15.47,  -1.88>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #286 
atom(<-13.65, -15.50,  -7.67>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #287 
atom(< -5.82, -11.49,  -5.49>, 1.20, rgb <0.87, 0.87, 0.90>, 0.0, jmol) // #288 
atom(<-15.13,  -3.40, -11.28>, 1.20, rgb <0.87, 0.87, 0.90>, 0.0, jmol) // #289 
atom(<-10.16,  -7.42,  -5.49>, 1.20, rgb <0.87, 0.87, 0.90>, 0.0, jmol) // #290 
atom(<-10.79,  -7.47, -11.28>, 1.20, rgb <0.87, 0.87, 0.90>, 0.0, jmol) // #291 
atom(< -7.67, -10.71,  -4.77>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #292 
atom(<-13.28,  -4.18, -10.55>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #293 
atom(<-12.01,  -8.20,  -4.77>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #294 
atom(< -8.94,  -6.69, -10.55>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #295 
atom(< -9.19,  -9.08,  -4.70>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #296 
atom(<-11.76,  -5.82, -10.49>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #297 
atom(<-13.53,  -9.83,  -4.70>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #298 
atom(< -7.42,  -5.06, -10.49>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #299 
atom(< -6.74,  -5.15,  -6.29>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #300 
atom(<-14.21,  -9.74,  -0.50>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #301 
atom(<-11.08,  -5.73,  -6.29>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #302 
atom(< -9.87,  -9.17,  -0.50>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #303 
atom(< -8.33,  -6.71,  -6.27>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #304 
atom(<-12.62,  -8.18,  -0.48>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #305 
atom(<-12.67,  -4.16,  -6.27>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #306 
atom(< -8.28, -10.73,  -0.48>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #307 
atom(<-14.36,  -4.50,  -3.66>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #308 
atom(< -6.59, -10.40,  -9.44>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #309 
atom(<-10.02,  -6.38,  -3.66>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #310 
atom(<-10.93,  -8.51,  -9.44>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #311 
atom(<-14.81,  -4.45,  -1.47>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #312 
atom(< -6.13, -10.45,  -7.26>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #313 
atom(<-10.47,  -6.43,  -1.47>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #314 
atom(<-10.48,  -8.46,  -7.26>, 0.74, rgb <1.00, 0.00, 0.00>, 0.0, jmol) // #315 
atom(< -9.39, -11.71,  -2.85>, 0.65, rgb <0.19, 0.31, 0.97>, 0.0, jmol) // #316 
atom(<-11.56,  -3.18,  -8.64>, 0.65, rgb <0.19, 0.31, 0.97>, 0.0, jmol) // #317 
atom(< -5.04,  -7.20,  -2.85>, 0.65, rgb <0.19, 0.31, 0.97>, 0.0, jmol) // #318 
atom(<-15.90,  -7.70,  -8.64>, 0.65, rgb <0.19, 0.31, 0.97>, 0.0, jmol) // #319 
atom(< -9.48, -12.75,  -2.94>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #320 
atom(<-11.47,  -2.14,  -8.73>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #321 
atom(< -5.14,  -6.16,  -2.94>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #322 
atom(<-15.81,  -8.74,  -8.73>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #323 
atom(< -8.73, -11.35,  -3.60>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #324 
atom(<-12.22,  -3.54,  -9.39>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #325 
atom(< -4.38,  -7.56,  -3.60>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #326 
atom(<-16.56,  -7.33,  -9.39>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #327 
atom(< -8.98, -11.43,  -1.92>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #328 
atom(<-11.97,  -3.46,  -7.70>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #329 
atom(< -4.64,  -7.48,  -1.92>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #330 
atom(<-16.31,  -7.42,  -7.70>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #331 
atom(< -7.57,  -9.07,  -5.97>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #332 
atom(<-13.38,  -5.83,  -0.18>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #333 
atom(<-11.91,  -9.84,  -5.97>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #334 
atom(< -9.03,  -5.05,  -0.18>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #335 
atom(< -8.21,  -5.23,  -4.87>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #336 
atom(<-12.74,  -9.66, -10.66>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #337 
atom(<-12.56,  -5.64,  -4.87>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #338 
atom(< -8.39,  -9.25, -10.66>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #339 
atom(< -6.27, -10.84,  -2.68>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #340 
atom(<-14.68,  -4.06,  -8.47>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #341 
atom(<-10.61,  -8.07,  -2.68>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #342 
atom(<-10.34,  -6.82,  -8.47>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #343 
atom(<-11.12, -11.43,  -3.99>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #344 
atom(< -9.82,  -3.47,  -9.77>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #345 
atom(< -6.78,  -7.48,  -3.99>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #346 
atom(<-14.17,  -7.41,  -9.77>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #347 
atom(<-10.57, -10.02,  -3.07>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #348 
atom(<-10.38,  -4.88,  -8.85>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #349 
atom(< -6.23,  -8.89,  -3.07>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #350 
atom(<-14.72,  -6.00,  -8.85>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #351 
atom(<-11.20, -11.18,  -0.91>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #352 
atom(< -9.75,  -3.71,  -6.70>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #353 
atom(< -6.85,  -7.73,  -0.91>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #354 
atom(<-14.09,  -7.17,  -6.70>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #355 
atom(<-12.59, -10.94,  -1.98>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #356 
atom(< -8.36,  -3.96,  -7.77>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #357 
atom(< -8.25,  -7.98,  -1.98>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #358 
atom(<-12.70,  -6.92,  -7.77>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #359 
atom(<-11.84, -12.56,  -1.86>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #360 
atom(< -9.11,  -2.33,  -7.64>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #361 
atom(< -7.50,  -6.35,  -1.86>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #362 
atom(<-13.45,  -8.54,  -7.64>, 0.40, rgb <1.00, 1.00, 1.00>, 0.0, jmol) // #363 
atom(< -8.14,  -9.59,  -5.17>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #364 
atom(<-12.81,  -5.30, -10.96>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #365 
atom(<-12.48,  -9.32,  -5.17>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #366 
atom(< -8.47,  -5.57, -10.96>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #367 
atom(< -7.76,  -5.68,  -5.79>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #368 
atom(<-13.19,  -9.22,   0.00>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #369 
atom(<-12.10,  -5.20,  -5.79>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #370 
atom(< -8.84,  -9.69,   0.00>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #371 
atom(<-14.72,  -3.88,  -2.60>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #372 
atom(< -6.23, -11.01,  -8.39>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #373 
atom(<-10.37,  -6.99,  -2.60>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #374 
atom(<-10.57,  -7.90,  -8.39>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #375 
atom(<-10.73, -11.10,  -3.01>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #376 
atom(<-10.22,  -3.79,  -8.80>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #377 
atom(< -6.38,  -7.81,  -3.01>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #378 
atom(<-14.57,  -7.08,  -8.80>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #379 
atom(<-11.64, -11.48,  -1.88>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #380 
atom(< -9.31,  -3.41,  -7.67>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #381 
atom(< -7.30,  -7.43,  -1.88>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #382 
atom(<-13.65,  -7.46,  -7.67>, 0.68, rgb <0.25, 0.25, 0.25>, 0.0, jmol) // #383 
