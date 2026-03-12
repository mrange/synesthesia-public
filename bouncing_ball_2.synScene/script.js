// State
var pos = [0, 2, 0];      // start near top
var vel = [1, 0, 0.5]; // small horizontal nudge so it doesn't just go straight down

const INNER_R = 1.0;
const OUTER_R = 4.0;
const BOUNCE_R = OUTER_R - INNER_R; // 3.0 — the center of the inner sphere can reach here
const GRAVITY = 9.8;
const DAMPING = 1.0;

function update(dt) {
  // Gravity
  vel[1] -= GRAVITY * dt;

  // Integrate position
  pos[0] += vel[0] * dt;
  pos[1] += vel[1] * dt;
  pos[2] += vel[2] * dt;

  // Collision with outer sphere
  dist = Math.sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
  if (dist > BOUNCE_R) {
    // Normal pointing inward (from wall toward center)
    const nx = pos[0] / dist;
    const ny = pos[1] / dist;
    const nz = pos[2] / dist;

    // Push position back to surface
    pos[0] = nx * BOUNCE_R;
    pos[1] = ny * BOUNCE_R;
    pos[2] = nz * BOUNCE_R;

    // Reflect velocity: v' = v - 2(v·n)n
    const vDotN = vel[0]*nx + vel[1]*ny + vel[2]*nz;
    vel[0] = (vel[0] - 2 * vDotN * nx) * DAMPING;
    vel[1] = (vel[1] - 2 * vDotN * ny) * DAMPING;
    vel[2] = (vel[2] - 2 * vDotN * nz) * DAMPING;
  }

  setUniform('u_sphere_pos', pos[0], pos[1], pos[2]);
  setUniform('u_angle', 0, 0);
}