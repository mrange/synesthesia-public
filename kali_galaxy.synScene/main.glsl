// This file is released under CC0 1.0 Universal (Public Domain Dedication).
// To the extent possible under law, Mårten Rånge has waived all copyright
// and related or neighboring rights to this work.
// See <https://creativecommons.org/publicdomain/zero/1.0/> for details.


#ifdef KODELIFE
uniform sampler2D pass0;
uniform sampler2D pass1;

const float
  global_div        = 3.
;
const vec2
  bass_limit        = vec2(0.5,1)
, base_xy           = vec2(2,1)
, base_zw           = vec2(0,0)
, base_mod          = vec2(.7,1)
, feedback_rot      = vec2(0,0)
, feedback_strength = vec2(0.1,.4)
, feedback_zoom     = vec2(0.9,0.8)
, flash_xy          = vec2(0,10)
, shape_xy          = vec2(0.5, 2)
, rotation_speed    = vec2(.1,.1)

;
const vec3
  feedback_min   = vec3(5,1,9)/9.
, feedback_max   = vec3(1)
, flash_col      = vec3(1,2,3)/3.
;
#endif

float beat() {
#ifdef KODELIFE
  return pow(1.-fract(TIME),2.);
#else
  return dot(pow(vec2(syn_BassLevel,syn_BassHits), bass_pow), bass_mix);
#endif
}

float freq(float x) {
#ifdef KODELIFE
  x=fract(x);
  return exp(-3.*x*x)*(1.-sqrt(fract(TIME)));
#else
  return smoothstep(0.4, .5, texture(syn_Spectrum,x).z);
#endif
}

mat2 rot(float a) {
  float c=cos(a),s=sin(a);
  return mat2(c,s,-s,c);
}

vec4 fpass0() {
  float
      z
    , d
    , s
    , D
    , l
    , L
    , b=smoothstep(bass_limit.x, bass_limit.y, beat())
#ifdef KODELIFE
    , T=(TIME+floor(TIME)+sqrt(fract(TIME)))*.1
#else
    , T=dot(vec2(syn_BassTime,TIME),rotation_speed)
#endif
    ;
  mat2
      M = rot(T)
    ;
  vec3
      I=normalize(vec3(_uvc, 1))
    ;
  vec4
      o=vec4(0)
    , P
    , q
    , F=vec4(flash_col*3.,0)
    , U=vec4(base_xy, base_zw)
    ;
  for(
      int i=0
    ; i < 55
    ; ++i
    ) {
    q = vec4(z*I, shape_xy.x);
    q.z -= shape_xy.y;
    q.yw *= M;
    q.wx *= M;
    q.zw *= M;
    P = q;
    s = 1.;
    d = 1e3;
    for(
        int j=0
      ; j < 7
      ; ++j
      ) {
      d = min(d, length(P) * s);
      D = dot(P, P);
      s *= D;
      P = abs(P.ywzx) / D - .1 * vec4(7, 5, 9, 5);
    }
    L=dot(q,q);
    l=sqrt(L);
    L*=L;
    P=U;
    P.w+=3. * length(q.xy);
    P = 1. + cos(base_mod.x * P + base_mod.y + 8. * l-freq(sin(l-5.)));
    o += mix(flash_xy.x, flash_xy.y,b*b)/L*F;
    o +=P.w * P / d;
    z += .6 * d + 1E-3;
  }

  return (o / (global_div*1E3));
}

vec4 fpass1() {
  float
    b=smoothstep(.5, 1., beat())
  ;
  vec2
    q0=_uv
  , pp=_uvc*mix(feedback_zoom.x,feedback_zoom.y, b)*rot(mix(feedback_rot.x,feedback_rot.y, b))
  , q1=pp*RENDERSIZE.y/RENDERSIZE+.5
  ;
  vec4
    col0=tanh(texture(pass0, q0))
  , col1=(texture(pass1, q1))
  ;
  return vec4((col0+mix(feedback_strength.x, feedback_strength.y, b*b)*mix(vec4(feedback_min,0), vec4(feedback_max,0), b)*col1).xyz,1.);
}

vec4 renderMain() {
  if (PASSINDEX==0) {
    return fpass0();
  } else if (PASSINDEX==1) {
    return fpass1();
  } else {
    return texture(pass1,_uv);
  }
}
