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

vec4 renderMain() {
  float
    , i
    , z=0.
    , d
    , k
    , B=beat()
    , t=dot(speed_control,vec2(TIME, syn_BassTime))
    , c=cos(t*.1)
    , s=sin(t*.1)
  ;

  vec2
    C=_xy
  , r=RENDERSIZE
  ;
  vec3
    I=normalize(vec3(C-.5*r,r.y))
  , o=vec3(0)
  ;
  vec4
    p
  , P
  ;

  mat2
    R=mat2(c,s,-s,c)
  ;

  for(i=0.;i<79.;++i) {
    p=vec4(z*I,.2),

    p.z-=3.;

    p.xw*=R;
    p.yw*=R;
    p.zw*=R;

    p*=k=9./dot(p,p);

    P=p-=.5*t;

    p=abs(p-round(p));

    d=abs(
     min(
       min(
           sqrt(min(min(dot(p.xz,p.xz),dot(p.yz,p.yz)), dot(p.xy,p.xy)))
         , length(p)-.2
         )
       , min(p.w,min(p.x,min(p.z,p.y)))+.05)
       )/k;

    p=1.+sin(P.z+log2(k)+vec4(1,3,6,0));

    o+=vec3(1,2,3)*exp(.7*k+6.*B-6.);
    o+=p.w*p.xyz/max(d,1e-3);
    z+=.8*d+1e-3;
  }

  o=tanh(o/1e4)/.9;
  o=mix(o,vec3(0),isnan(o));
  return vec4(o,1);
}
