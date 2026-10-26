#define R(p)p=Q*dot(Q,p)+cross(Q,p)

// License: Unknown, author: Claude Brezinski, found: https://mathr.co.uk/blog/2017-09-06_approximating_hyperbolic_tangent.html
vec4 tanh_approx(vec4 x) {
  //  Found this somewhere on the interwebs
  //  return tanh(x);
  vec4 x2 = x*x;
  return clamp(x*(27.0 + x2)/(27.0+9.0*x2), -1.0, 1.0);
}

vec4 renderMain() {
  float d,e,s,z;
  vec4 p,O=vec4(0),P,o;
  vec3 
    Q=normalize(sin(speed/vec3(11,13,17)))
  , R=normalize(vec3(_uvc,1.+syn_BassLevel/9.))
  , D
  , H=13.*smoothstep(vec3(low_level),vec3(high_level),vec3(syn_BassHits,syn_MidHits,syn_MidHighHits))
  ;
  for(int i=0; i<77; ++i) {
    p=vec4(z*R-vec3(0,0,1),.6);
    R(p.yzw);
    R(p.zwx);
    R(p.wxy);
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
    p=e*e+3.+3.*cos(13.*dot(P.xy,P.zw)+vec4(2,7,0,4));
    D=H/abs(sqrt(D)-glow_width);
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
  O=tanh_approx(O);
  O/=.9;
  
  p=_loadMedia();
  O=mix(O,p,p.w*media_opacity*mix(lum_opacity,1.,dot(vec3(0.299, 0.587, 0.114),p.xyz)));
  return vec4(O.xyz,1);
}