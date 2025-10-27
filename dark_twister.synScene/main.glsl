// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
const vec4 hsv2rgb_K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
vec3 hsv2rgb(vec3 c) {
  vec3 p = abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www);
  return c.z * mix(hsv2rgb_K.xxx, clamp(p - hsv2rgb_K.xxx, 0.0, 1.0), c.y);
}
// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
//  Macro version of above to enable compile-time constants
#define HSV2RGB(c)  (c.z * mix(hsv2rgb_K.xxx, clamp(abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www) - hsv2rgb_K.xxx, 0.0, 1.0), c.y))

#define ROT(a)      mat2(cos(a), sin(a), -sin(a), cos(a))

const float
  TAU=2.*PI
, max_dist  = 10.
, max_steps = 90.
, norm_eps  = 5e-4
, tolerance = 5e-4
;

const vec2
  PathA=vec2(1,sqrt(.5))/6.
, PathB=vec2(6,4)
;

const vec3
  light_col = HSV2RGB(vec3(.58,.7,2.))
, ref_col   = HSV2RGB(vec3(.62,.7,1.))
;

float
  g_G
, g_H
;
mat2
  g_R
;

vec3 path(float z){
  return vec3(PathB*sin(PathA*z),z);
}

vec3 dpath(float z){
  return vec3(PathB*PathA*cos(PathA*z),1);
}

vec3 ddpath(float z){
  return -vec3(PathB*PathA*PathA*sin(PathA*z),0);
}

void warpWorld(inout vec3 p){
  vec3
    warp = path(p.z)
  , dwarp = normalize(dpath(p.z))
  ;
  p.xy -= warp.xy;
  p -= dwarp*dot(vec3(p.xy, 0), dwarp)*.5*vec3(1,1,-1);
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
float roundedX(vec2 p, float w, float r) {
  p = abs(p);
  return length(p-min(p.x+p.y,w)*.5) - r;
}

float hash(vec2 co) {
  return fract(sin(dot(co.xy ,vec2(12.9898,58.233))) * 13758.5453);
}

float df(vec2 p, float a) {
  const mat2
    R45=ROT(radians(45.))
  ;
  mat2
    R=ROT(a)
  ;
  float
    i
  , d=1e3
  , g
  , S=.5
  , H
  ;
  vec2
    P
  , s=vec2(0.)
  , I
  ;
  for(i=0.;i<3.;++i) {
    p*=R;
    s+=S*sign(p);
    p=abs(p);
    p-=S;
    d=min(d,roundedX(p,.25*S,.05*(S)));
    P=p;
    P*=R45;
    I=sign(P);
    P=abs(P);
    P-=.125*S;
    H=hash(s+.123*I);
    g=length(P)+.05*(smoothstep(1.,.9,sin(1.*TIME+.5*a+TAU*fract(8667.*H))));
    if (g<g_G) {
      g_G=g;
      g_H=H;
    }
    d=min(d,g);
    S*=.45;
  }
  return d;
}

float df(vec3 p)  {
  warpWorld(p);
  p.xy*=g_R;
  return df(p.xy,p.z);
}

vec3 nf(vec3 p) {
  vec2 e=vec2(norm_eps,0);
  return normalize(vec3(
    df(p+e.xyy)-df(p-e.xyy)
  , df(p+e.yxy)-df(p-e.yxy)
  , df(p+e.yyx)-df(p-e.yyx)
  ));
}


// License: MIT, author: Inigo Quilez, found: https://www.iquilezles.org/www/articles/spherefunctions/spherefunctions.htm
float raySphereDensity(vec3 ro, vec3 rd, vec4 sph, float dbuffer) {
  float ndbuffer = dbuffer/sph.w;
  vec3  rc = (ro - sph.xyz)/sph.w;
  float b = dot(rd,rc);
  float c = dot(rc,rc) - 1.0;
  float h = b*b - c;
  if(h<0.0) return 0.0;
  h = sqrt(h);
  float t1 = -b - h;
  float t2 = -b + h;
  if(t2<0.0 || t1>ndbuffer) return 0.0;
  t1 = max(t1, 0.0);
  t2 = min(t2, ndbuffer);
  float i1 = -(c*t1 + b*t1*t1 + t1*t1*t1/3.0);
  float i2 = -(c*t2 + b*t2*t2 + t2*t2*t2/3.0);
  return (i2-i1)*(3.0/4.0);
}

float raymarch(vec3 O, vec3 I, float initz, out float iter) {
  float
    i
  , z=initz
  , d
  ;
  for(i=0.;i<max_steps;++i) {
    d=df(I*z+O);
    if(z>max_dist || d < tolerance) {
      break;
    }
    z+=.9*d;
  }
  iter = i;
  return z;
}

vec3 render3D() {
  float
    i
  , B=dot(pow(vec2(syn_BassLevel, syn_BassHits),bass_pow),bass_mix)
  , D
  , G
  , F
  , H
  , L
  , z=mod(TIME,100.)
  ;
  vec2
    p=2.*_uvc;
  vec3
    O=path(z)
  , L0=path(z+flash_dist)
  , Z=normalize(dpath(z))
  , X=normalize(cross(vec3(0,1,0)-ddpath(z),Z))
  , Y=cross(X,Z)
  , I=normalize(p.x*X+p.y*Y+2.*Z)
  , col=vec3(0)
  , P
  , N
  , R
  , LD0
  ;
  vec4
    pcol
  ;
  g_R=ROT(.3*TIME);
  g_G=1e3;
  z=raymarch(O,I,.1,i);
  G=g_G;
  H=g_H;
  P=O+I*z;
  LD0=normalize(L0-P);
  N=nf(P);
  F=1.+dot(N,I);
  F=mix(.1,1.,F*F*F)*smoothstep(1.,.99,F);
  R=reflect(I,N);
  if (z<max_dist) {
    col+=F*pow(max(dot(R,LD0),0.),190.)*20./dot(O-L0,O-L0)*ref_col;
  }

  D=raySphereDensity(O,I,vec4(L0,2.),z);
  L=(1.+dot(normalize(O-L0),I));
  col+=D*D*light_col*mix(flash_mix.x, flash_mix.y, B);
  col+=1e-4/max(G*G,2e-6)*hsv2rgb(vec3(.7+.2*H,.9,.1));
  col-=1e-1*L*vec3(2,3,1);
  col=max(col,0.);
  col=tanh(col);
  col=sqrt(col);
  pcol=texture(syn_FinalPass,_uv);
  col=mix(col,pcol.xyz,.3);
  return col;
}

vec3 render2D() {
  vec2
    p =_uvc
  ;

  float
    aa=sqrt(.5)/RENDERSIZE.y
  , d=df(p,.1*TIME)
  ;

  return vec3(1)*smoothstep(aa,-aa,d);
}

vec4 renderMain() {
  vec3
    col
  ;

  col=render3D();


  return vec4(col,1);

}
