
float dot2(vec2 p) {
  return dot(p,p);
}

vec4 jumpFloodPass(sampler2D tex, vec2 p, vec2 xy, int stp, ivec2 ixy) {
  ivec2
    pos
  ;
  vec2
    nseed0
  , nseed1
  , seed0
  , seed1
  ;
  float
    dist0
  , dist1
  , ndist0
  , ndist1
  ;
  int
    x
  , y
  ;
  
  vec4 
    t = texelFetch(tex, ixy, 0)
  ;
  seed0 = t.xy;
  seed1 = t.zw;
  dist0 = seed0==vec2(0) ? 1e6 : dot2(xy-seed0);
  dist1 = seed1==vec2(0) ? 1e6 : dot2(xy-seed1);
  
  for(y=-1;y<=1;++y)
  for(x=-1;x<=1;++x) {
    pos=ixy+stp*ivec2(x,y);
    if(x==0&&y==0) continue;
    /*
    if(any(lessThan(pos, ivec2(0))) || any(greaterThanEqual(pos, ivec2(resolution)))) 
      continue;
    */

    t=texelFetch(tex,pos,0);

    nseed0=t.xy;
    nseed1=t.zw;

    ndist0=dot2(xy-nseed0);
    ndist1=dot2(xy-nseed1);

    if(nseed0!=vec2(0)&&dist0>ndist0) {
      dist0=ndist0;
      seed0=nseed0;
    }
    if(nseed1!=vec2(0)&&dist1>ndist1) {
      dist1=ndist1;
      seed1=nseed1;
    }
  }
  
  return vec4(seed0,seed1);
}


vec4 renderMain() {
  vec2
    p =_uvc
  ;

  float
    aa=sqrt(.5)/RENDERSIZE.y
  , d=length(p)-.5
  ;

  return vec4(vec3(1)*smoothstep(aa,-aa,d),1);
}
