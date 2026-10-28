#define ROT(p)p=Q*dot(Q,p)+cross(Q,p)

// License: Unknown, author: Claude Brezinski, found: https://mathr.co.uk/blog/2017-09-06_approximating_hyperbolic_tangent.html
float tanh_approx(float x) {
  //  Found this somewhere on the interwebs
  //  return tanh(x);
  float x2=x*x;
  return clamp(x*(27.+x2)/(27.+9.*x2),-1., 1.);
}

// License: Unknown, author: Claude Brezinski, found: https://mathr.co.uk/blog/2017-09-06_approximating_hyperbolic_tangent.html
vec4 tanh_approx(vec4 x) {
  //  Found this somewhere on the interwebs
  //  return tanh(x);
  vec4 x2=x*x;
  return clamp(x*(27.+x2)/(27.+9.*x2),-1., 1.);
}

vec4 fpass0() {
  float d, e, s, z=0.,B=smoothstep(low_level, high_level, syn_BassHits),C=color_chroma_onbeat*B;
  vec3
    D
  , Q=u_rotor
  , H=u_hits
  ;
  vec4
    o, p, P
  , O=vec4(0)
  , R=normalize(vec4(_uvc,1.+syn_BassLevel*screen_bounce,0))
  , S=vec4(0,0,-1,mix(woff, woff_onbeat, B))
  ,PP=texelFetch(pass0, ivec2(_xy), 0)
  ;
  ROT(S.yzw);
  ROT(S.zwx);
  ROT(S.wxy);
  ROT(R.yzw);
  ROT(R.zwx);
  ROT(R.wxy);
  for(int i=0; i<77&&z<6.; ++i) {
    p=z*R+S;
    P=p;
    s=1.;
    D=vec3(1e3);
    for(int j=0;j<5;++j) {
      D=min(D,vec3(dot(p.xw,p.xw),dot(p.yw,p.yw),dot(p.zw,p.zw)));
      p-=2.*round(p*.5);
      p*=e=2./dot(p,p);
      s*=e;
    }
    p=abs(p/s);
    d=min(length(p.xy)+3e-4,min(p.z,p.w)+3e-3);
    p=e*e+3.+3.*cos(13.*dot(P.xy,P.zw)+vec4(color_chroma,C,color_tune*length(P)+4.));
    D=H/max(abs(sqrt(D)-glow_width), 1e-4);
    o=p.w/(7.*e*e*e*d+d)*p;
    o+=D.x*vec4(3,1,2,0);
    o+=D.y*vec4(1,2,3,0);
    o+=D.z*vec4(2,1,3,0);
    o*=exp(-z*z*.5);
    O+=o;
    z+=.6*d+1e-3;
    }
  ;
  O/=2e4;
  O-=length(_uv-.5)/5.*vec4(2,3,1,0);

  O=mix(O+over_glow*PP,PP,tanh_approx(dot(vec3(0.299, 0.587, 0.114),PP.xyz))*after_glow);
  return vec4(O.xyz,1);
}

vec4 fpass1() {
  vec4
    p0=texelFetch(pass0, ivec2(_xy), 0)
  , m =_loadMedia()
  , O =p0
  ;
  O=max(O,0.);
  O=tanh_approx(O);
  O/=.9;
  O=mix(O,m,m.w*media_opacity*mix(lum_opacity,1.,dot(vec3(0.299, 0.587, 0.114),m.xyz)));
  return O;
}

vec4 renderMain() {
  switch(PASSINDEX) {
  case 0: return fpass0();
  default: return fpass1();
  }
}
