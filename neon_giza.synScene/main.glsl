const float
  TAU                 = 2.*PI
, MAX_RAY_LENGTH_HI   = 48.0
, TOLERANCE_HI        = 0.0001
, MAX_RAY_MARCHES_HI  = 70.
, NORM_OFF            = 0.001
;

#define ROT(a)              mat2(cos(a), sin(a), -sin(a), cos(a))
#define SCA(a)              vec2(sin(a), cos(a))

// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
const vec4 hsv2rgb_K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
vec3 hsv2rgb(vec3 c) {
  vec3 p = abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www);
  return c.z * mix(hsv2rgb_K.xxx, clamp(p - hsv2rgb_K.xxx, 0.0, 1.0), c.y);
}
// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
//  Macro version of above to enable compile-time constants
#define HSV2RGB(c)  (c.z * mix(hsv2rgb_K.xxx, clamp(abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www) - hsv2rgb_K.xxx, 0.0, 1.0), c.y))

//#define PATHA (0.33*vec2(0.1147, 0.2093))
//#define PATHB (0.33*vec2(13.0, 3.0))
const vec2
  PATHA = (TAU/300.*vec2(3,2))
, PATHB = (vec2(2, 3))
;
  
vec3 path(float z) {
  return vec3(sin(z*PATHA)*PATHB, z);
}

vec3 dpath(float z) {
  return vec3(PATHA*PATHB*cos(PATHA*z), 1.0);
}

vec3 ddpath(float z) {
  return vec3(-PATHA*PATHA*PATHB*sin(PATHA*z), 0.0);
}

// License: Unknown, author: Claude Brezinski, found: https://mathr.co.uk/blog/2017-09-06_approximating_hyperbolic_tangent.html
float tanh_approx(float x) {
  //  Found this somewhere on the interwebs
  //  return tanh(x);
  float x2 = x*x;
  return clamp(x*(27.0 + x2)/(27.0+9.0*x2), -1.0, 1.0);
}

float beat() {
#ifdef KODELIFE
  return pow(1.-fract(TIME),2.);
#else
  return tanh_approx(dot(pow(vec2(syn_BassLevel,syn_BassHits), bass_pow), bass_mix));
#endif
}

float wave(float x) {
#ifdef KODELIFE
  return (.5+.25*sin(x*10.+TIME));
#else
  return (texture(syn_Spectrum,x).w-.5)*2.;
#endif
}

void rot(inout vec2 p, float a) {
  float
    c=cos(a)
  , s=sin(a)
  ;
  p=vec2(c*p.x+s*p.y,c*p.y-s*p.x);
}

// License: MIT, author: Inigo Quilez, found: https://www.iquilezles.org/www/articles/smin/smin.htm
float pmin(float a, float b, float k) {
  float h = clamp(0.5+0.5*(b-a)/k, 0.0, 1.0);
  return mix(b, a, h) - k*h*(1.0-h);
}

// License: CC0, author: Mårten Rånge, found: https://github.com/mrange/glsl-snippets
float pmax(float a, float b, float k) {
  return -pmin(-a, -b, k);
}

float rayPlane(vec3 ro, vec3 rd, vec4 p) {
  return -(dot(ro,p.xyz)+p.w)/dot(rd,p.xyz);
}

float length4(vec2 p) {
  p*=p;
  return sqrt(sqrt(dot(p,p)));
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/articles/distfunctions2d/
float flatTorus( vec3 p, vec2 t) {
  p=p.xzy;
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length4(q)-t.y;
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/articles/distfunctions2d/
float cappedTorus(vec3 p, vec2 sc, vec2 t) {
  float
    ra = t.x
  , rb = t.y
  ;
  p.x = abs(p.x);
  float k = (sc.y*p.x>sc.x*p.y) ? dot(p.xy,sc) : length(p.xy);
  return sqrt(dot(p,p) + ra*ra - 2.0*ra*k) - rb;
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/articles/distfunctions2d/
float arc(vec2 p, vec2 sc, float ra, float rb) {
  // sc is the sin/cos of the arc's aperture
  p.x = abs(p.x);
  return ((sc.y*p.x>sc.x*p.y) ? length(p-sc*ra) :
                                abs(length(p)-ra)) - rb;
}

#define ROTY(a)               \
  mat3(                       \
    +cos(a) , 0.0 , +sin(a) \
  , 0.0     , 1.0 , 0.0     \
  , -sin(a) , 0.0 , +cos(a) \
  )

#define ROTZ(a)               \
  mat3(                       \
    +cos(a) , +sin(a) , 0.0   \
  , -sin(a) , +cos(a) , 0.0   \
  , 0.0     , 0.0     , 1.0   \
  )

#define ROTX(a)               \
  mat3(                       \
    1.0 , 0.0     , 0.0       \
  , 0.0 , +cos(a) , +sin(a)   \
  , 0.0 , -sin(a) , +cos(a)   \
  )

const float
  hoff      = 0.725;
const mat3
  roty      = ROTY(radians(10.0))
;
const vec3
  sunDir     = normalize(vec3(0.0,-0.01, 1.0))*roty
, lightPos   = vec3(0.0, -60.0, -200.0)*roty
, sunColor   = HSV2RGB(vec3(hoff+0.0, 0.9, 0.0005))
, topColor   = HSV2RGB(vec3(hoff+0.0, 0.9, 0.0001))
, glowColor0 = HSV2RGB(vec3(hoff+0.0, 0.9, 0.0001))
, glowColor2 = HSV2RGB(vec3(hoff+0.3, 0.95, 0.001))
, diffColor  = HSV2RGB(vec3(hoff+0.0, 0.9, .25))
, glowCol1   = HSV2RGB(vec3(hoff+0.2, 0.85, 0.0125))
;

vec2 planeCoord(vec3 p, vec3 c, vec3 up, vec4 dim) {
  vec3
    d = p - c
  , xx = (cross(up,dim.xyz))
  , yy = (cross(xx,dim.xyz))
  ;
  return vec2(dot(d,xx), dot(d,yy));
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
float triIso(vec2 p, vec2 q) {
  p.x = abs(p.x);
  vec2 a = p - q*clamp( dot(p,q)/dot(q,q), 0.0, 1.0 );
  vec2 b = p - q*vec2( clamp( p.x/q.x, 0.0, 1.0 ), 1.0 );
  float s = -sign( q.y );
  vec2 d = min( vec2( dot(a,a), s*(p.x*q.y-p.y*q.x) ),
                vec2( dot(b,b), s*(p.y-q.y)  ));
  return -sqrt(d.x)*sign(d.y);
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/articles/intersectors/
vec2 rayBox(vec3 ro, vec3 rd, vec3 boxSize, out vec3 outNormal)  {
  vec3 m = 1.0/rd; // can precompute if traversing a set of aligned boxes
  vec3 n = m*ro;   // can precompute if traversing a set of aligned boxes
  vec3 k = abs(m)*boxSize;
  vec3 t1 = -n - k;
  vec3 t2 = -n + k;
  float tN = max( max( t1.x, t1.y ), t1.z );
  float tF = min( min( t2.x, t2.y ), t2.z );
  if( tN>tF || tF<0.0) return vec2(-1.0); // no intersection
  outNormal = (tN>0.0) ? step(vec3(tN),t1)  : // ro ouside the box
                         step(t2,vec3(tF))  ;  // ro inside the box
  outNormal *= -sign(rd);
  return vec2( tN, tF );
}

vec3 sky(vec3 ro, vec3 rd) {
  float
    hd = max(abs(rd.y+0.15), 0.00066)
  ;
  vec3
    col = sunColor/(1.0+1e-5 - dot(sunDir, rd))
  ;
  col += 100.0*glowColor0/sqrt(hd);
  col += glowColor2/(hd);
  return col;
}

vec3 glow(vec3 ro, vec3 rd, float beat) {
  vec3
    bn
  , bro = ro
  ;
  bro.y += -1000.0+70.0;
  vec2
    bi = rayBox(bro, rd, vec3(90.0, 1000.0, 90.0), bn)
  ;
  float
    lightDist = distance(lightPos, ro)
  ;
  vec3
    lightDir  = normalize(lightPos-ro)
  ;
  float
    g3        = 1.0+0.00001 - dot(lightDir, rd)
  ;

  vec3
    col = 8.*glowColor0/(g3)
  , rrd = rd*transpose(roty)*ROTX(0.027)
  ;

  if (bi != vec2(-1)) {
    float bdi = tanh_approx(0.00125*(bi.y-bi.x));
    col += 1000.0*glowColor0*(bdi/max(rrd.y, 0.005));
  }

  float sx = abs(rrd.x);
  rrd.y += mix(0.00, 0.0125, beat)*(wave(0.5*sx));
  col += 20.0*glowColor0/(abs(mix(mix(20.0, 0.25, beat)*rrd.y*rrd.y, abs(rrd.y), tanh_approx(4.0*sx)))+mix(2.0, 0.5, beat)*sx*sx*sx+0.0001);

  col *= mix(1.0, 4.0, beat);
  return col;
}

vec3 side(vec3 col, vec3 ro, vec3 rd, float t, vec4 dim, vec3 c) {
  vec3 n = dim.xyz;

  vec3 p = ro + rd*t;

  vec3 r = reflect(rd, n);
  vec3 ldiff = p - lightPos;
  vec3 ld = normalize(ldiff);
  vec3 rcol0 = sky(p, r);
  float dcol = max(dot(ld, n), 0.0);
  dcol *= dcol;
  vec2 pp = planeCoord(p, c, vec3(0.0, 1.0, 0.0), dim);
  vec2 p0 = pp;
  vec2 p1 = pp;
  const vec2 tri =vec2(485, sqrt(3.0)*356.0);
  float d0 = triIso(p0, tri);
  float d1 = triIso(p1, 0.11*tri);
  float d = d0;
  vec3 bcol = col;
  float hf = smoothstep(-600.0, -400.0, p.y);
  vec3 pcol = 3.0*diffColor*dcol;
  pcol += rcol0;
  pcol = mix(clamp(col, 0.0, 0.1), pcol, hf);
  col = mix(col, pcol, smoothstep(9., 0., d));
  col += topColor/max(0.00005*(d1-1.), 0.000025)*hf;
  return col;
}

vec3 pyramid(vec3 col, vec3 ro, vec3 rd) {
  const mat3
    rotx = ROTX(radians(-51.8))
  , rr0  = rotx
  , rr1  = rr0*ROTY(PI/2.0)
  ;

  const vec3
    n0   = normalize(vec3(.0, 0.0, 1.0))*rr0
  , c0   = vec3(0.0)
  , n1   = normalize(vec3(.0, 0.0, 1.0))*rr1
  , c1   = vec3(0.0)
  ;

  const vec4
    dim0 = vec4(n0, -dot(c0, n0))
  , dim1 = vec4(n1, -dot(c1, n1))
  ;
  float
    t0  = rayPlane(ro, rd , dim0)
  , t1  = rayPlane(ro, rd , dim1)
  ;

  if (t1 > 0.0) {
    col = side(col, ro, rd, t1, dim1, c1);
  }
  if (t0 > 0.0) {
    col = side(col, ro, rd, t0, dim0, c0);
  }


  return col;
}

float g_gd;

void warpWorld(inout vec3 p){
  vec3
    warp = path(p.z)
  , dwarp = normalize(dpath(p.z))
  , ddwarp = ddpath(p.z)
  ;
  p.xy -= warp.xy;
  p -= dwarp*dot(vec3(p.xy, 0), dwarp)*0.5*vec3(1,1,-1);
  rot(p.xy,-50.*ddwarp.x);
}

// License: Unknown, author: Unknown, found: don't remember
float hash(float co) {
  return fract(sin(co*12.9898) * 13758.5453);
}

vec2 dfArcs(vec3 p) {
  const mat2
    R = ROT(-PI/4.)
  ;
  const vec2
    sc=SCA(PI/16.)
  ;

  vec3
    p1 = p
  ;
  float
    n1 = round(.25*p1.z)
  , h1 = hash(mod(n1+37.5,75.)-37.5)
  , sh1 = -1.0+2.0*h1
  ;
  p1.z-=4.*n1;
  vec3
    p2 = p1
    ;
  rot(p1.xy,.25*TIME*sh1);
  p1.xy = abs(p1.xy);
  p1.xy *= R;
  float
    d1 = cappedTorus(p1, sc, 3.0*vec2(1.0, 0.075))
  , d2 = flatTorus(p2, 3.0*vec2(1.0, 0.03))
  ;
  d2 = max(d2, -d1);
  return vec2(d1, d2);
}

float df(vec3 p, float t) {
  warpWorld(p);
  vec3 p0 = p;
  p0.y -= 3.0;
  p0.y = -p0.y;
  vec2 dd1 = dfArcs(p);
  float
    d0 = arc(p0.xy, SCA(PI/6.0), 6.0, 0.8)
  , d = d0
  ;
  d = pmax(d, -(dd1.y-0.25), 0.125);
  d = min(d, dd1.x);
  d = min(d, dd1.y);
  float gd = dd1.x;
  t = max(t-MAX_RAY_LENGTH_HI*0.5, 0.0);
  g_gd = min(g_gd, gd+t*t*1E-4);

  return d;
}

float df(vec3 p) {
  return df(p, 0.0);
}

vec3 normal(vec3 pos) {
  vec2  eps = vec2(NORM_OFF,0.0);
  return normalize(vec3(
    df(pos+eps.xyy) - df(pos-eps.xyy)
  , df(pos+eps.yxy) - df(pos-eps.yxy)
  , df(pos+eps.yyx) - df(pos-eps.yyx)
  ));
}

float rayMarch(vec3 ro, vec3 rd, float initt, out float iter) {
  const float tol = TOLERANCE_HI;
  float
    t = initt
  , i
  ;
  for (i = 0.; i < MAX_RAY_MARCHES_HI; ++i) {
    float d = df(ro + rd*t, t);
    if (d < TOLERANCE_HI || t > MAX_RAY_LENGTH_HI) {
      break;
    }
    t += d;
  }
  iter = i;
  return t;
}


vec3 render0(vec3 ro, vec3 rd, float beat) {
  const
    vec3 ro_ = vec3(0.0, 0.0, -2700.0)*roty
  ;
  const
    mat3 rrd = ROTX(-0.2)*ROTY(0.)
  ;
  const float
    rdd=  3.0
  , mm = 4.0
  ;
  rd *= rrd;
  ro = ro_;

  vec3
    glowCol = glow(ro, rd, beat)
  , col = sky(ro, rd)
  ;

  col = pyramid(col, ro, rd);
  col += glowCol;
  col = clamp(col, 0., 10.);

  return col;
}

vec3 render1(vec3 ro, vec3 rd, vec2 sp) {
  g_gd = 1E3;
  float
    B = beat()
  , iter
  , t = rayMarch(ro, rd, 0.0, iter)
  , gd = g_gd
  ;
  vec3
    ggcol = (glowCol1)/(max(gd, 1e-3))
  , skyCol = render0(ro, rd, B)
  , col = skyCol
  ;

  float tt = t/MAX_RAY_LENGTH_HI;
  tt -= 0.33;
  tt = clamp(tt, 0.0, 1.0);
  float sfo = 1.-exp(-9.*tt*tt);

  if (t < MAX_RAY_LENGTH_HI) {
    vec3
      p = ro+rd*t
    , n = normal(p)
    , r = reflect(rd, n)
    ;
    float
      fre0 = 1.0+dot(rd, n)
    , fre = fre0
    ;
    fre *= fre;

    float
      dif = max(dot(sunDir, n),0.)
    , ao = 1.0-iter/float(MAX_RAY_MARCHES_HI)
    , fo = mix(0.2, 0.5, ao)
    ;

    vec3 wp = p;
    warpWorld(wp);
    vec2 dd = dfArcs(wp);

    vec3 rcol = vec3(0.0);
    float hit = min(dd.x, dd.y);
    if (hit > .05) {
      g_gd = 1E3;
      float riter;
      float rt = rayMarch(p, r, 0.5,riter);
      float rgd = g_gd;
      vec3 rggcol = 0.5/(max(rgd*rgd, 1e-3))*(glowCol1);
      rcol = rggcol;
      rcol *= smoothstep(0.66, 0.1, tt);
      if (rt < MAX_RAY_LENGTH_HI) {
        rcol += diffColor*.2;
      } else {
        rcol += .5*render0(p, r, B);
      }
    }

    rcol += 4.0*smoothstep(0.,0.01,gd)/max(dd.x*dd.x, 0.01)*(diffColor+0.5)*glowCol1;
    col = vec3(0);
    col += dif*dif*fo*1e5*sunColor*diffColor;
    col += fre*rcol;
  }

  col = mix(col, skyCol, sfo);

  col += ggcol;

  vec3 rrd = rd*transpose(roty)*ROTX(0.027);
  float flash = dot(rrd, normalize(vec3(0.0, -0.2, -1.0)))+1.0005;
  col += (0.01*vec3(0.5, 0.25, 1.0))*smoothstep(.5, 1., B*B)/flash;


  return col;
}


vec3 effect(vec2 p, vec2 pp, vec2 q) {
  const vec3
    up = vec3(0, 1, 0)
  ;

  float
    z = 5.0*(mod(TIME, 60.)-30.0)
  ;

  vec2
    np = p + 4.0/RENDERSIZE.y;

  vec3
    ro = path(z)
  , ww = normalize(dpath(z))
  , dd = ddpath(z)
  , uu = normalize(cross(up+10.*dd, ww))
  , vv = cross(ww,uu)
  , rd = normalize(-p.x*uu + p.y*vv + 2.*ww)
  , col= render1(ro, rd, p)
  ;

  col -= 0.025*vec3(2,3,1)*(.25+length(pp));
  col *= smoothstep(1.5, 0.5, length(pp));
  col = clamp(col, 0.0, 4.0);
  vec4
    pcol=texture(syn_FinalPass,q)
  ;
  col =tanh(col);
  col= sqrt(col);
#ifndef KODELIFE
  vec4
    mcol=_loadMedia()
  ;
  col=mix(col,mcol.xyz,mcol.w*media_opacity*media_multiplier);
#endif
  col=mix(col,pcol.xyz,motion_blur);
  return col;
}

vec4 renderMain() {
  vec2 q = _uv;
  vec2 p = 2.*_uvc;
  vec3 col = effect(p, -1.+2.*q,q);

  return vec4(col, 1.0);
}

