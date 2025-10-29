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

mat2 rot(float a) {
  float
      c=cos(a)
    , s=sin(a)
    ;
  return mat2(c,s,-s,c);
}

vec4 renderMain() {
  float
    i
  , z=0.
  , d
  , k
  , B=beat()
  , t=dot(speed_control,vec2(TIME, syn_BassTime))
  , T=.1*t
  ;

  vec2
    C=_xy
  , r=RENDERSIZE
  , U=2.*_uvc
  , S=sin(U*sine_distortion.x)
  ;
  vec3
    I=normalize(vec3(C-.5*r,(1.+syn_BassHits*sine_distortion.y*dot(U,U)*S.x*S.y)*r.y))
  , o=vec3(0)
  ;
  vec4
    p
  , P
  , M
  ;

  mat2
    R0=rot(T)
  , R1=rot(rotation_speed.x*T)
  , R2=rot(rotation_speed.y*T)
  ;

  for(i=0.;i<79.;++i) {
    p=vec4(z*I,wcut),

    p.z-=3.;

    p.xw*=R0;
    p.yw*=R1;
    p.zw*=R2;

    p*=k=zoom_factor/dot(p,p);

    P=p-=translation_speed*t;

    p=abs(p-round(p));

    d=abs(
     min(
       min(
           sqrt(min(min(dot(p.xz,p.xz),dot(p.yz,p.yz)), dot(p.xy,p.xy)))
         , length(p)-.2
         )
       , min(p.w,min(p.x,min(p.z,p.y)))+.05)
       )/k;

    p=1.+sin(P.z+log2(k)+vec4(effect_base0,effect_base1));

    o+=visual_mix.x*exp(.7*k+6.*B-6.)*3.*flash_color;
    o+=visual_mix.y*p.w/max(d,1e-3)*p.xyz;
    z+=.8*d+1e-3;
  }

  o=tanh(o/1e4)/0.9;
  o=mix(o,vec3(0),isnan(o));
  M=_loadMedia();  
  o=mix(o,M.xyz,M.w*media_transparency*media_multiplier);
  return vec4(o,1);
}
