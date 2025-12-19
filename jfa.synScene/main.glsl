#ifdef KODELIFEX
const float 
  lum_cutoff = .5
;
#endif

float dot2(vec2 p) {
  return dot(p,p);
}

float texelLum(sampler2D tex, ivec2 ixy) {
  return dot(vec3(.299, .587, .114),texelFetch(tex,ixy,0).xyz);
}

vec4 image(sampler2D tex, vec2 xy, ivec2 ixy) {
   ivec2
    sz=textureSize(tex, 0)
  ;
  ixy += sz/2-ivec2(.5*RENDERSIZE);
  return texelFetch(tex,ixy,0);
}

vec4 init(vec2 xy, ivec2 ixy) {
  vec2
    seed0=vec2(0)
  , seed1=vec2(ixy)
  ;
  
  vec4
    t=image(syn_Media, xy, ixy)
  ;

  float
    l=dot(vec3(.299, .587, .114),t.xyz*t.w)
  ;
  
  if(l>lum_cutoff) {
    return vec4(seed1, seed0);  
  } else {
    return vec4(seed0, seed1);  
  }
  
}

vec4 jfa(sampler2D tex, int stp, vec2 xy, ivec2 ixy) {
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

vec4 dist(sampler2D tex, vec2 xy, ivec2 ixy) {
  vec4
    t=texelFetch(tex, ixy, 0)
  , i=image(syn_Media, xy, ixy)
  ;
  vec2 
    seed0=t.xy
  , seed1=t.zw
  ;

  float 
    dist0=distance(xy,seed0)
  , dist1=distance(xy,seed1)
  , dist =(dist0>dist1?dist0:-dist1)*2./RENDERSIZE.y
  , aa   = sqrt(2.)/RENDERSIZE. y
  ;

  vec3 
    col=vec3(smoothstep(aa,-aa,dist-.0125))
  ;
  
  if (dist < 0.0) {
    col = i.xyz*i.w;
  }
  
//  col=vec3(sin(dist*80.));
  return vec4(col, 1);
}

vec4 renderMain() {
  ivec2 
    ixy=ivec2(_xy)
  ;
  switch(PASSINDEX) {
  case 0:
    return init(_xy, ixy);
  case 1:
    return jfa(pass0, 16, _xy, ixy);
  case 2:
    return jfa(pass1, 8, _xy, ixy);
  case 3:
    return jfa(pass2, 4, _xy, ixy);
  case 4:
    return jfa(pass3, 2, _xy, ixy);
  case 5:
    return jfa(pass4, 1, _xy, ixy);
  default:
    return dist(pass5, _xy, ixy);
  }
}
