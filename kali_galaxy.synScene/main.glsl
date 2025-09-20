// This file is released under CC0 1.0 Universal (Public Domain Dedication).
// To the extent possible under law, Mårten Rånge has waived all copyright
// and related or neighboring rights to this work.
// See <https://creativecommons.org/publicdomain/zero/1.0/> for details.

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
  return smoothstep(0.4, .5, texture(syn_Spectrum,x).y);
#endif
}

float wave(float x) {
#ifdef KODELIFE
  x=fract(x);
  return (.5+.25*sin(x*10.+TIME));
#else
  return texture(syn_Spectrum,x).w;
#endif
}

mat2 rot(float a) {
  float c=cos(a),s=sin(a);
  return mat2(c,s,-s,c);
}

float length4(vec2 p) {
  p*=p;
  return pow(dot(p,p),.25) ;
}

vec4 renderMain() {
  float
      z
    , d
    , s
    , D
    , B=beat()
    , T=syn_BassTime/2.+TIME
    ;
  mat2
      M = rot(.1*T)
    ;
  vec3
      I=normalize(vec3(_uvc, 1))
    ;
  vec4
      o=vec4(0)
    , P
    , q
    ;
    B*=B;
  for(
      int i=0
    ; i < 55
    ; ++i
    ) {
    q = vec4(z*I, .5);
    q.z -= 2.;
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
        d = min(d, length(P.xyz) * s);
        D = dot(P, P);
        s *= D;
        P = abs(P.ywzx) / D - .1 * vec4(7, 5, 9, 5);
      }

    P = 1. + cos(.7 * vec4(2, 1, 0, 3. * length(q.xy)) + 1. + 8. * length(q)-freq(length(q)));
    o += B*0./dot(q,q)*vec4(1,2,3,0);
    o +=P.w * P / d;
    z += .6 * d + 1E-3;
  }

  return tanh(o / 3E3);
}
