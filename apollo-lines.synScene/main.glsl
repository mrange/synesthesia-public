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
    p4 = vec4(p, 0.1)
  ;

  p4.yw = p4.yw*g_rot;
  p4.zw = g_rot*p4.zw;

  return apollonian(p4, g_scale, w, off);
}

vec3 glowmarch(vec3 col, vec3 ro, vec3 rd, float tinit) {
  float
    t = tinit
  , off
  , d
  ;

  for (int i = 0; i < 60; ++i) {
    d = df(ro + rd*t, 6E-5+t*t*2E-3, off);
    col += 1E-9/max(d*d, 1E-8)*(palette(log(off))+5E-2);
    t += .5*max(d, 1E-4);
    if (t > .5) break;
  }

  return col;
}

vec3 render(vec3 col, vec3 ro, vec3 rd) {
  col = glowmarch(col, ro, rd, 1E-2);
  return col;
}


vec3 effect(vec2 p, vec2 pp, vec2 q) {
  float
    tm  = TIME
  ;

  g_scale = mix(1.85, 1.5, .5-.5*cos(TAU*tm/1600.));
  g_rot = ROT(tm*TAU/800.);

  vec2
    s=sin(9.*p+TIME)*length(p)
  ;

  vec3
    ro = pos
  , ZZ = normalize(dpos)
  , XX = normalize(cross(vec3(0,1,0), ZZ))
  , YY = cross(ZZ, XX)
  , rd = normalize(-p.x*XX + p.y*YY + mix(2.,2.05,.5-.5*s.x*s.y)*ZZ)
  , col
  ;

  col = .1/max(.5-rd.y+.1*rd.x*rd.x, .1)*palette(5.+.1*rd.y);
  col = render(col, ro, rd);
  col *= smoothstep(1.707, .707, length(pp));
  col -= 3E-2*(.3+dot(pp,pp))*vec3(2,3,1);
  col = aces_approx(col);
  col = sqrt(col);

  return col;
}

vec4 renderMain() {
  vec3 col = effect(2.*_uvc, -1.+2.*_uv, _uv);
  return vec4(col,1.0);
}

