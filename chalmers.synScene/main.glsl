mat2 rot(float a) {
  float c=cos(a),s=sin(a);
  return mat2(c,s,-s,c);
}

float decode(vec4 c) {
  const vec4
    decoder=vec4(1,1./255.,1./65025.,0)
  ;
  const float
    halfoff=+.5/(255.)+.5/(65025.)
  ;

  return dot(c,decoder)-halfoff;
}

float dlinTex(sampler2D tex, vec2 p) {
  vec2
    sz=vec2(textureSize(tex, 0))
  , tx=clamp(p*sz-.5,vec2(0), sz-2.)
  , ntx=floor(tx)
  , ftx=fract(tx)
  ;
  ivec2
    itx=ivec2(ntx)
  ;

  vec4
    c00 = texelFetch(tex, itx+ivec2(0,0), 0)
  , c01 = texelFetch(tex, itx+ivec2(0,1), 0)
  , c10 = texelFetch(tex, itx+ivec2(1,0), 0)
  , c11 = texelFetch(tex, itx+ivec2(1,1), 0)
  ;

  float
    d00 = decode(c00)
  , d01 = decode(c01)
  , d10 = decode(c10)
  , d11 = decode(c11)
  , d   = mix(
      mix(d00, d01, ftx.y)
    , mix(d10, d11, ftx.y)
    , ftx.x
    )
  ;

  return -1.+2.*d;
}

float dtexture(sampler2D tex, vec2 p) {
  const float D=4.;
  p*=.5;
  p+=.5;
#ifdef KODELIFE
  p.y=1.-p.y;
#endif
  return D*dlinTex(tex, p);
}

float dchalmers(vec2 p) {
  p.y*=4.;
  return dtexture(t_chalmers,p);
}

float dchalmers(vec3 p, float h, float r) {
  h-=r;
  float
    d=dchalmers(p.xy)+r
  ;
  vec2
    w=vec2(d, abs(p.z)-h);
  return min(max(w.x,w.y),.0)+length(max(w,.0))-r;
}

float df(vec4 p) {
  vec4
    p0=p
  ;
  vec3
    D=sqrt(vec3(dot(p.xz,p.xz),dot(p.yz,p.yz), dot(p.xy,p.xy)))
  ;
//  p0.yz*=rot(2.*p0.x);
  p0.yz=abs(p0.yz);
  p0.yz-=.125;
  p=abs(p);
  return min(
      min(min(min(D.x,D.y),D.z)-.02, min(p.w,min(p.x,min(p.z,p.y)))+.02)
    , dchalmers(p0.xyz,.05,.02)
  );
}

// License: Unknown, author: Claude Brezinski, found: https://mathr.co.uk/blog/2017-09-06_approximating_hyperbolic_tangent.html
vec3 tanh_approx(vec3 x) {
  //  Found this somewhere on the interwebs
  //  return tanh(x);
  vec3 x2 = x*x;
  return clamp(x*(27.0 + x2)/(27.0+9.0*x2), -1.0, 1.0);
}


vec4 renderMain() {
  vec4
    O
  , p
  , P
  , RepZ=vec4(stretch,.8,.8,1e3)
  , IRepZ=1./RepZ
  , color_off=vec4(6,color_twist,8,6)
  ;

  vec2
    r2=RENDERSIZE
  , p2=2.*_uvc
  ;

  float
    i
  , z=0.
  , d
  , k
  , F=1.-syn_BassHits
  , aa=sqrt(2.)/r2.y
  ;

  mat2
    rot0=mat2(R0)
  , rot1=mat2(R1)
  , rot2=mat2(R2)
  ;

  vec3
    rd=normalize(vec3(p2,2))
  , ro=vec3(0,0,-distance)
  , o=vec3(0)
  ;

  for(i=0.;i<79.&&z<30.;++i) {
    p=vec4(z*rd+ro,w_offset),
    p.xw*=rot0;
    p.yw*=rot1;
    p.zw*=rot2;
    k=9./dot(p,p);
    p*=k;
    p-=.5*effect_time;
    P=p;
    p-=floor(p*IRepZ+.5)*RepZ;
    d=abs(df(p))/k;
    p=1.2+sin(F+P.z+log2(k)+color_off);
    o+=(flash_strength*exp(.7*k-4.*F)*vec4(4,8,12,0)+p.w*p/max(d,1e-3)).xyz;
    z+=.8*d+1e-3;
  }

  o/=4e4;

  if(logo_zoom>.1&&logo_transparency>.05) {
    p2/=logo_zoom;
    d=dchalmers(p2-logo_pos);
    d=min(d,abs(d-.005)-.00125);
    d*=logo_zoom;
    
    i=texture(syn_LevelTrail,.01+level_trail*d).y;
    o+=vec3(1,3,2)*logo_transparency*i*smoothstep(.18*logo_zoom,0.,d)*smoothstep(aa,-aa,-d);
  }
  o=tanh_approx(o);
  return vec4(o.xyz,1);
}
