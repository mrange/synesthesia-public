open System
open System.Numerics

let inline fract1 (v : float32) = v - floor v
let inline fract3 (v : Vector3) = Vector3 (fract1 v.X, fract1 v.Y, fract1 v.Z)

// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
let hsv2rgb (h : float32) (s : float32) (v : float32) : Vector3 =
  let zero  = Vector3.Zero
  let one   = Vector3.One
  let three = 3.F*one
  let xyz   = Vector3(3.F, 2.F, 1.F)/3.F
  v*Vector3.Lerp(
      one
    , Vector3.Clamp(
        Vector3.Abs(6.F*(fract3 (Vector3(h)+xyz))-three)-one
      , zero
      , one
      )
    , s
    )
let inline HSV2RGB (h, s, v) =
  let rgb = hsv2rgb (float32 h) (float32 s) (float32 v)
  printfn "HSV2RGB(%.2f,%.2f,%.2f): vec3(%.2f, %.2f, %.2f)" h s v rgb.X rgb.Y rgb.Z

HSV2RGB (0.85, 0.90, 0.5)
HSV2RGB (0.06, 0.90, 0.5)
HSV2RGB (0.03, 0.90, 0.5)
HSV2RGB (0.60, 0.90, 0.5)
HSV2RGB (0.64, 0.98, 0.5)
HSV2RGB (0.0 , 0.98, 0.5)
HSV2RGB (0.1 , 0.0 , 0.5)
