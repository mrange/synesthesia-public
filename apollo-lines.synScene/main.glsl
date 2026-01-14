// This file is released under CC0 1.0 Universal (Public Domain Dedication).
// To the extent possible under law, Mårten Rånge has waived all copyright
// and related or neighboring rights to this work.
// See <https://creativecommons.org/publicdomain/zero/1.0/> for details.

#define ROT(a)      mat2(cos(a), sin(a), -sin(a), cos(a))

const float
  TAU       = 2.*PI
;

mat2 g_rot;
float g_scale;

// License: Unknown, author: Matt Taylor (https://github.com/64), found: https://64.github.io/tonemapping/
vec3 aces_approx(vec3 v) {
  const float
    a = 2.51
  , b = 0.03
  , c = 2.43
  , d = 0.59
  , e = 0.14
  ;

  v = max(v, 0.);
  v *= d;
  return clamp((v*(a*v+b))/(v*(c*v+d)+e), .0, 1.);
}

// License: MIT, author: Inigo Quilez, found: https://www.iquilezles.org/www/articles/smin/smin.htm
float pmin(float a, float b, float k) {
  float h = clamp(0.5+0.5*(b-a)/k, 0.0, 1.0);
  return mix(b, a, h) - k*h*(1.0-h);
}

float pmax(float a, float b, float k) {
  return -pmin(-a, -b, k);
}

vec3 palette(float a) {
  return 1.+sin(vec3(0,7,2)+a);
}

float apollonian(vec4 p, float s, float w, out float off) {
  float
    z=1.
  , k
  , d
  ;

  for(int i=0; i<6;++i) {
    p-=2.*floor(p*.5+.5);
    k= s/dot(p,p);
    p*=k;
    z*=k;
  }

  vec4
    sp=p/z
  , ap=abs(sp)-w
  ;

  d=pmax(ap.w, ap.y, w*11.);
  d=min(d, pmax(ap.x, ap.z, w*11.));
  off=length(sp);
  
  return d;
}

float df(vec3 p, float w, out float off) {
  vec4
    p4 = vec4(p, woffset)
  ;

  p4.yw = p4.yw*g_rot;
  p4.zw = g_rot*p4.zw;

  return apollonian(p4, g_scale, w, off);
}

vec3 glowmarch(vec3 ro, vec3 rd, float tinit, out float tlast) {

  vec3 
    col=vec3(0)
  , p
  ;
  float 
    t=tinit
  , off
  , d
  ;
  
  for (int i=0; i<60; ++i) {
    p=ro+rd*t;
    d=df(p, 6E-5+t*t*2E-3, off);
    col+=1E-9/max(d*d, 1E-8)*(palette(log(off))+5E-2);
    t+=.5*max(d, 1E-4);
    if (t>2.) break;
  }
  
  tlast = t;
  return col;
}

// License: MIT, author: Inigo Quilez, found: https://www.iquilezles.org/www/articles/spherefunctions/spherefunctions.htm
float ray_sphere_density(vec3 ro, vec3 rd, vec4 sph, float dbuffer) {
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

vec3 render(vec3 ro, vec3 rd) {
  float 
    tlast
  , den
  ;
  vec3 
    col=.1/max(.5-rd.y+.1*rd.x*rd.x, .1)*palette(5.+.1*rd.y)
  ;
  col+=glowmarch(ro, rd, 1E-2, tlast);
  den=ray_sphere_density(ro,rd,vec4(vec3(0,0,0),glow_radii),tlast);
  col+=den*den*2.*syn_BassLevel*palette(-.5);
  return col;
}

vec3 effect(vec2 p, vec2 pp) {
  g_scale= mix(1.85, 1.5, .5-.5*cos(TAU*TIME/1600.));
  g_rot  = ROT(TIME*TAU/800.);

  vec2
    s=sin(9.*p+TIME)*length(pp)
  ;

  vec3
    ro = pos
  , ZZ = normalize(dpos)
  , XX = normalize(cross(vec3(0,1,0)-ddpos, ZZ))
  , YY = cross(ZZ, XX)
  , rd = normalize(-p.x*XX + p.y*YY + mix(2.,2.1,.5-.5*s.x*s.y)*ZZ)
  , col= render(ro, rd)
  ;

  col*=smoothstep(1.707, .707, length(pp));
  col-=3E-2*(.3+dot(pp,pp))*vec3(2,3,1);
  col= aces_approx(col);
  col= sqrt(col);

  vec4 
    m=_loadMedia()
  ;
  col=mix(col,m.xyz,(.5+p.y)*m.w*media_opacity);

  return col;
}

vec4 renderMain() {
  vec3 col = effect(2.*_uvc, -1.+2.*_uv);
  return vec4(col,1.0);
}

