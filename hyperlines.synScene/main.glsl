// This file is released under CC0 1.0 Universal (Public Domain Dedication).
// To the extent possible under law, Mårten Rånge has waived all copyright
// and related or neighboring rights to this work.
// See <https://creativecommons.org/publicdomain/zero/1.0/> for details.

const float
  TAU           =2.*PI
, PI_2          =.5*PI
#ifdef KODELIFE
, spectrum_boost    =1.
#endif
;

#define ROT(a)      mat2(cos(a), sin(a), -sin(a), cos(a))

vec4 spectrum(float x) {
#ifdef KODELIFE
  return vec4(1.);
#else
  return textureLod(syn_Spectrum,x,0.);
#endif
}

// License: Unknown, author: Claude Brezinski, found: https://mathr.co.uk/blog/2017-09-06_approximating_hyperbolic_tangent.html
vec3 tanh_approx(vec3 x) {
  //  Found this somewhere on the interwebs
  //  return tanh(x);
  vec3 x2 = x*x;
  return clamp(x*(27.0 + x2)/(27.0+9.0*x2), -1.0, 1.0);
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/articles/intersectors/
float ray_cylinder(vec3 ro, vec3 rd, float r) {
  // Tweaked it a bit
  float
      a=dot(rd.xy, rd.xy)
    , b=dot(ro.xy, rd.xy)
    , c=dot(ro.xy, ro.xy) - r*r
    , h=b*b - a*c
    , t
    ;

  if(h<.0) return -1.;

  h = sqrt(h);
  t = (-b+h)/a;

  return t;
}

// License: MIT, author: Pascal Gilcher, found: https://www.shadertoy.com/view/flSXRV
float atan_approx(float y, float x) {
  float
    cosatan2 = x/(abs(x)+abs(y))
  , t = PI_2-cosatan2*PI_2
  ;
  return y<0.?-t:t;
}

// License: MIT OR CC-BY-NC-4.0, author: mercury, found: https://mercury.sexy/hg_sdf/
vec2 mod_polar(inout vec2 p, float repetitions) {
  float
    angle = TAU/repetitions
  , a = atan_approx(p.y, p.x) + angle/2.
  , r = length(p)
  , c = floor(a/angle)
  ;
  a = mod(a,angle) - angle/2.;
  p = vec2(cos(a), sin(a))*r;
  if (abs(c) >= (repetitions/2.)) c = abs(c);
  return vec2(c,r);
}

// License: Unknown, author: Unknown, found: don't remember
float hash(vec2 co) {
  co+=.1234;
  return fract(sin(dot(co.xy ,vec2(12.9898,58.233)))*13758.5453);
}

// License: Unknown, author: knarkowicz, found: https://www.shadertoy.com/view/XtlSD7
vec2 crt_distort(vec2 q) {
  q = _uv*2. - 1.;
  vec2
    o = crt_effect*q.yx/vec2(6,4)
  ;
  q = q + q*o*o;
  return q*.5 + .5;
}

// License: Unknown, author: knarkowicz, found: https://www.shadertoy.com/view/XtlSD7
float vig(vec2 q) {
  float
    v = q.x*q.y*(1.0 - q.x)*(1.0 - q.y)
  ;
  v = clamp(pow(16.*v, .3), 0., 1.);
  return v;
}

float segmentz(vec3 p, float hl) {
  p.z=abs(p.z)-hl;
  float
    d0=length(p)
  , d1=length(p.xy)
  ;

  return p.z>0.?d0:d1;
}

vec3 hyperspace(vec3 RO, vec3 RD, float FO) {
  float
    H0
  , H1
  , REP
  , ci
  , d
  , n1
  , fo
  ;

  vec2
    N
  ;

  vec3
    o=mix(0e-4,3e-4, bass_thump)*(vec3(1,4,16))/(1.+1e-3-(RD.z)+RD.x*RD.x)
  , p
  , c
  ;
  vec4
    O
  ;

  mat2 
    R=ROT(rot_speed)
  ;

  for (float j=2.;j<9.;++j) {
    REP=j*j+3.;
    ci=ray_cylinder(RO,RD,j);
    p=ci*RD+RO;
    p.xy*=R;
    N=mod_polar(p.xy,REP);

    H0=hash(.123*vec2(j,N.x));
    p.z-=speed*(1.+H0*H0*H0);
    fo=exp(-2e-3*ci*ci);
    c=vec3(N.y,0,p.z);
    c-=vec3(p.xy,floor(p.z+.5));
    for(float i=-1.;i<=1.;++i) {
      n1=floor(i+p.z+.5);
      H1=hash(vec2(H0,n1));
      O=1.+sin(PI*H1+vec4(0,1,8,4)-base_color);
      d=segmentz(c-vec3(0,0,i),.35*H1*H1+.1);
      H1=spectrum_boost*spectrum(mix(.975,.025,H1)).z;
      o+=
          5e-3
        * H1
        / max(d,4e-4*ci+.02*FO)
        * fo
        * (O.w+.1)
        * O.xyz;
    }
  }

  return o;
}


vec4 renderMain() {
  vec2
    r=RENDERSIZE
  , q=crt_distort(_uv)
  , p=(2.*q-1.)*r/r.y
#ifdef KODELIFE
  , t2=.2*time*vec2(sqrt(2.),1.)
#endif
  ;

  vec3
#ifdef KODELIFE
    RO=vec3(sin(t2),-2)
  , LA=vec3(0)
  , Z =normalize(LA-RO)
  , X =normalize(cross(Z,vec3(.2*cos(t2)+vec2(0,1.),0)))
  , Y =cross(X,Z)
  , RD=normalize(2.*Z+p.y*Y-p.x*X)
#else
    RO=cam_RO
  , LA=vec3(0)
  , RD=normalize(2.*cam_Z+p.y*cam_Y-p.x*cam_X)
#endif
  , o =vec3(0)
  ;
  
  vec4
    m=_loadMedia();
  ;

  o=hyperspace(RO,RD,0.);

  o*=mix(1.,vig(q),crt_effect);
  o-=3e-2*vec3(3,2,1)*length(p+.25);
  o=max(o,0.);
  o=tanh_approx(1.125*o);
  o=sqrt(o)-.05;

  o*=mix(1.,1.5+.5*sin(_xy.y*TAU/6.), crt_effect);


  o=mix(o,m.xyz,m.w*mix(1.,dot(m.xyz,vec3(0.299, 0.587, 0.114)),media_mix_mode)*media_opacity);
  return vec4(o,1);
}
