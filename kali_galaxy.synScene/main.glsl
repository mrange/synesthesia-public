#ifdef KODELIFE
uniform sampler2D pass0;
uniform sampler2D pass1;
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

vec4 fpass0() {
  float
      z
    , d
    , s
    , D
    , l
    , L
    , b=smoothstep(.5, 1., beat())
#ifdef KODELIFE
    , T=TIME+floor(TIME)+sqrt(fract(TIME))
#else
    , T=syn_BassTime+TIME
#endif
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
        d = min(d, length(P) * s);
        D = dot(P, P);
        s *= D;
        P = abs(P.ywzx) / D - .1 * vec4(7, 5, 9, 5);
      }
    L=dot(q,q);
    l=sqrt(L);
    L*=L;
    P = 1. + cos(.7 * vec4(2, 1, 0, 3. * length(q.xy)) + 1. + 8. * l-freq(sin(l-5.)));
    o += 0.*b/L*vec4(1,2,3,0);
    o +=P.w * P / d;
    z += .6 * d + 1E-3;
  }

  return (o / 3E3);
}


vec4 fpass1() {
  float
    b=smoothstep(.5, 1., beat())
  ;
  vec2
    q=_uv
  , p=_uvc*mix(.9,.8, b)
  , qq=p*RENDERSIZE.y/RENDERSIZE+.5
  ;
  vec4
    col0=tanh(texture(pass0, q))
  , col1=(texture(pass1, qq))
  ;
  return vec4((col0+mix(0.00,0.06, b*b)*mix(vec4(5,1,9,0), vec4(9), b)*col1).xyz,1.);
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
