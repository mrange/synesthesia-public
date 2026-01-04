float beat() {
  return dot(pow(vec2(syn_BassLevel,syn_BassHits), bass_pow), bass_mix);
}


float dot2(vec2 p) {
  return dot(p,p);
}

float texelLum(sampler2D tex, ivec2 ixy) {
  return dot(vec3(.299, .587, .114),texelFetch(tex,ixy,0).xyz);
}

void rot(inout vec2 p, float a) {
  float
    c=cos(a)
  , s=sin(a)
  ;
  p=vec2(c*p.x+s*p.y,c*p.y-s*p.x);
}
vec4 image(sampler2D tex, vec2 xy, ivec2 ixy) {
  vec2
    p=(2.*xy-RENDERSIZE)/RENDERSIZE.y
  , sz=vec2(textureSize(tex,0))
  ;
  float
    b=beat()
  , l=length(p-dist_offset)
  ;
  
  p-=media_offset;
  
  rot(p, mix(dist_rot.x,dist_rot.y,b*l));
  p/=media_zoom*mix(dist_zoom.x,dist_zoom.y,b*l);
  p.x*=sz.y/sz.x;
  p+=.5;
  
  return textureLod(tex,clamp(p,0,1),0);
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

// License: Unknown, author: Claude Brezinski, found: https://mathr.co.uk/blog/2017-09-06_approximating_hyperbolic_tangent.html
vec3 tanh_approx(vec3 x) {
  //  Found this somewhere on the interwebs
  //  return tanh(x);
  vec3 x2 = x*x;
  return clamp(x*(27.0 + x2)/(27.0+9.0*x2), -1.0, 1.0);
}

vec4 dist(sampler2D tex, vec2 xy, ivec2 ixy) {
  vec4
    t=texelFetch(tex, ixy, 0)
  , i=image(syn_Media, xy, ixy)
  ;
  vec2 
    seed0=t.xy
  , seed1=t.zw
  , p=(2.*xy-RENDERSIZE)/RENDERSIZE.y
  ;

  float 
    dist0=distance(xy,seed0)
  , dist1=distance(xy,seed1)
  , dist =(dist0>dist1?dist0:-dist1)*2./RENDERSIZE.y
  , aa   = sqrt(2.)/RENDERSIZE. y
  , b    = beat()
  , f    = smoothstep(flash_mod.x,flash_mod.y,b)
  , freq = mix(beat_mod.x,beat_mod.y,b)
  ;
  
  //dist=min(dist,abs(dist-.3));
  vec3 
    col
  ;
  
  
  col=
     smoothstep(-aa,aa,(line_mod+cos(freq*dist))/freq*.5)
    /(1.+fade_mod*dist*dist)
    *(1.+sin(vec3(color_base,0)-TIME+dist*3.+(p.y)))
    ;

  if (dist < 0.) {
    col= mix(col,i.xyz,i.w);
    col*=2.*col;
  }

  //dist-=mix(.1,.2,b);
  col+=f*glow_mod*1e-3/max(dist*dist,1e-4);
  col=tanh_approx(col);
  col=max(col,0.);
  col=sqrt(col);

  
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
    return jfa(pass0, 512, _xy, ixy);
  case 2:
    return jfa(pass1, 256, _xy, ixy);
  case 3:
    return jfa(pass2, 128, _xy, ixy);
  case 4:
    return jfa(pass3, 64, _xy, ixy);
  case 5:
    return jfa(pass4, 32, _xy, ixy);
  case 6:
    return jfa(pass5, 16, _xy, ixy);
  case 7:
    return jfa(pass6, 8, _xy, ixy);
  case 8:
    return jfa(pass7, 4, _xy, ixy);
  case 9:
    return jfa(pass8, 2, _xy, ixy);
  case 10:
    return jfa(pass9, 1, _xy, ixy);
  default:
    return dist(pass10, _xy, ixy);
  }
}
