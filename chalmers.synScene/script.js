function mix(a,b,c) {
  return (1-c)*a+c*b;
}
function setRotUniform(nm,a) {
  const c=Math.cos(a);
  const s=Math.sin(a);
  setUniform(nm,c,s,-s,c);
}
function update(dt) {
  const t=mix(.5*syn_BassTime, speed,speed_mix);
  setRotUniform("effect_time",t);
  setRotUniform("R0",rot_speed_0*t);
  setRotUniform("R1",rot_speed_1*t);
  setRotUniform("R2",rot_speed_2*t);
    
}
