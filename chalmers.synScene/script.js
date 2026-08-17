function mix(a,b,c) {
  return (1-c)*a+c*b;
}
function setRotUniform(nm,a) {
  const c=Math.cos(a);
  const s=Math.sin(a);
  setUniform(nm,c,s,-s,c);
}

var R0=0;
var R1=0;
var R2=0;
var t =0;
function update(dt) {
  const tt=speed*mix(2.*syn_BassLevel*dt,dt, speed_mix);
  t+=tt;
  R0+=rot_speed_0*tt;
  R1+=rot_speed_1*tt;
  R2+=rot_speed_2*tt;
  setRotUniform("effect_time",t);
  setRotUniform("R0",R0);
  setRotUniform("R1",R1);
  setRotUniform("R2",R2);
    
}
