#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White}
camera {orthographic
  right -26.85*x up 13.57*y
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
#declare Rcell = 0.020;
#declare Rbond = 0.100;

#macro atom(LOC, R, COL, TRANS, FIN)
  sphere{LOC, R texture{pigment{color COL transmit TRANS} finish{FIN}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{torus{R, Rcell rotate 45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      torus{R, Rcell rotate -45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      translate LOC}
#end

atom(<  7.05,   2.12, -30.65>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #0 
atom(<  8.17,   0.63, -27.24>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #1 
atom(<  4.51,   1.82, -27.73>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #2 
atom(<  5.63,   0.34, -24.33>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #3 
atom(<  4.24,   1.95, -25.38>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #4 
atom(<  8.78,   2.81, -26.50>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #5 
atom(<  5.36,   0.47, -21.97>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #6 
atom(<  7.65,   4.29, -29.90>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #7 
atom(<  6.78,   2.25, -28.29>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #8 
atom(<  6.23,   2.51, -23.58>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #9 
atom(<  7.90,   0.77, -24.89>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #10 
atom(<  5.11,   3.99, -26.99>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #11 
atom(<  5.96,   2.64, -21.22>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #12 
atom(<  4.84,   4.12, -24.63>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #13 
atom(<  8.50,   2.94, -24.14>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #14 
atom(<  7.38,   4.42, -27.55>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #15 
atom(<  8.23,   3.07, -21.78>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #16 
atom(<  7.11,   4.55, -25.19>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #17 
atom(<  5.69,   2.77, -18.86>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #18 
atom(<  4.57,   4.25, -22.27>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #19 
atom(<  9.29,  -0.85, -23.84>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #20 
atom(< 10.42,  -2.33, -20.43>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #21 
atom(<  6.75,  -1.15, -20.92>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #22 
atom(<  7.87,  -2.63, -17.51>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #23 
atom(<  6.48,  -1.02, -18.56>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #24 
atom(< 11.02,  -0.16, -19.68>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #25 
atom(<  7.60,  -2.50, -15.15>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #26 
atom(<  9.90,   1.33, -23.09>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #27 
atom(<  9.02,  -0.72, -21.48>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #28 
atom(<  8.48,  -0.46, -16.76>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #29 
atom(< 10.14,  -2.20, -18.07>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #30 
atom(<  7.35,   1.03, -20.17>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #31 
atom(<  8.20,  -0.33, -14.40>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #32 
atom(<  7.08,   1.16, -17.81>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #33 
atom(< 10.75,  -0.03, -17.32>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #34 
atom(<  9.62,   1.46, -20.73>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #35 
atom(< 10.47,   0.10, -14.96>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #36 
atom(<  9.35,   1.59, -18.37>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #37 
atom(<  7.93,  -0.19, -12.05>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #38 
atom(<  6.81,   1.29, -15.45>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #39 
atom(<  1.97,   1.52, -24.81>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #40 
atom(<  3.09,   0.04, -21.41>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #41 
atom(< -0.58,   1.22, -21.90>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #42 
atom(<  0.54,  -0.26, -18.49>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #43 
atom(< -0.85,   1.35, -19.54>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #44 
atom(<  3.69,   2.21, -20.66>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #45 
atom(<  0.27,  -0.13, -16.13>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #46 
atom(<  2.57,   3.69, -24.07>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #47 
atom(<  1.69,   1.65, -22.46>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #48 
atom(<  1.15,   1.91, -17.74>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #49 
atom(<  2.81,   0.17, -19.05>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #50 
atom(<  0.03,   3.40, -21.15>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #51 
atom(<  0.88,   2.04, -15.38>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #52 
atom(< -0.25,   3.53, -18.79>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #53 
atom(<  3.42,   2.34, -18.30>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #54 
atom(<  2.30,   3.83, -21.71>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #55 
atom(<  3.15,   2.47, -15.94>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #56 
atom(<  2.02,   3.96, -19.35>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #57 
atom(<  0.60,   2.17, -13.02>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #58 
atom(< -0.52,   3.66, -16.43>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #59 
atom(<  4.21,  -1.45, -18.00>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #60 
atom(<  5.33,  -2.93, -14.59>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #61 
atom(<  1.67,  -1.74, -15.08>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #62 
atom(<  2.79,  -3.23, -11.67>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #63 
atom(<  1.39,  -1.61, -12.72>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #64 
atom(<  5.93,  -0.75, -13.84>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #65 
atom(<  2.52,  -3.10,  -9.32>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #66 
atom(<  4.81,   0.73, -17.25>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #67 
atom(<  3.94,  -1.32, -15.64>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #68 
atom(<  3.39,  -1.05, -10.93>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #69 
atom(<  5.06,  -2.80, -12.23>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #70 
atom(<  2.27,   0.43, -14.33>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #71 
atom(<  3.12,  -0.92,  -8.57>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #72 
atom(<  2.00,   0.56, -11.98>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #73 
atom(<  5.66,  -0.62, -11.49>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #74 
atom(<  4.54,   0.86, -14.89>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #75 
atom(<  5.39,  -0.49,  -9.13>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #76 
atom(<  4.27,   0.99, -12.54>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #77 
atom(<  2.85,  -0.79,  -6.21>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #78 
atom(<  1.73,   0.69,  -9.62>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #79 
atom(< -3.12,   0.92, -18.98>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #80 
atom(< -2.00,  -0.56, -15.57>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #81 
atom(< -5.66,   0.62, -16.06>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #82 
atom(< -4.54,  -0.86, -12.65>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #83 
atom(< -5.93,   0.75, -13.70>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #84 
atom(< -1.39,   1.61, -14.82>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #85 
atom(< -4.81,  -0.73, -10.29>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #86 
atom(< -2.52,   3.10, -18.23>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #87 
atom(< -3.39,   1.05, -16.62>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #88 
atom(< -2.27,  -0.43, -13.21>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #89 
atom(< -0.88,  -2.04, -12.16>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #90 
atom(<  0.25,  -3.53,  -8.76>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #91 
atom(< -3.42,  -2.34,  -9.24>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #92 
atom(< -2.30,  -3.83,  -5.84>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #93 
atom(< -3.69,  -2.21,  -6.89>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #94 
atom(<  0.85,  -1.35,  -8.01>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #95 
atom(< -2.57,  -3.69,  -3.48>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #96 
atom(< -0.27,   0.13, -11.41>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #97 
atom(< -1.15,  -1.91,  -9.80>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #98 
atom(< -0.03,  -3.40,  -6.40>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #99 
atom(< -8.20,   0.33, -13.14>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #100 
atom(< -7.08,  -1.16,  -9.73>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #101 
atom(<-10.75,   0.03, -10.22>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #102 
atom(< -9.62,  -1.46,  -6.82>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #103 
atom(< -8.48,   0.46, -10.78>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #104 
atom(< -7.35,  -1.03,  -7.38>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #105 
atom(< -4.84,  -4.12,  -2.92>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #106 
atom(< -8.50,  -2.94,  -3.41>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #107 
atom(< -7.38,  -4.42,   0.00>, 2.04, rgb <1.00, 1.00, 0.78>, 0.0, ase3) // #108 
atom(< -6.23,  -2.51,  -3.97>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #109 
atom(< -5.11,  -3.99,  -0.56>, 0.66, rgb <1.00, 0.05, 0.05>, 0.0, ase3) // #110 
atom(< -5.88,  -2.37,  -6.23>, 1.39, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #111 
atom(< -3.69,   1.11, -14.29>, 1.39, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #112 
atom(< -0.27,   1.48, -10.49>, 1.39, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #113 
atom(<  5.35,   2.67, -16.58>, 1.39, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #114 
atom(<  9.02,   1.49, -16.09>, 1.39, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #115 
atom(<  6.47,   1.19, -13.18>, 1.39, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #116 
atom(<  7.20,   2.70, -14.97>, 1.39, rgb <0.00, 0.41, 0.52>, 0.0, ase3) // #117 
