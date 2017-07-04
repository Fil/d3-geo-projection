import {atan2, halfPi, log} from "./math";

export function complexAtan(x, y) {
  var x2 = x * x,
      y_1 = y + 1,
      t = 1 - x2 - y * y;
  return [
   0.5 * ((x >= 0 ? halfPi : -halfPi) - atan2(t, 2 * x)),
    -0.25 * log(t * t + 4 * x2) +0.5 * log(y_1 * y_1 + x2)
  ];
}

export function complexDivide(a, b) {
  if (b[1])
    a = complexMul(a, b), b = complexNorm2(b);
  else
    b = b[0];
  return [
    a[0] / b,
    a[1] / b
  ];
}

export function complexMul(a, b) {
  return [
    a[0] * b[0] + a[1] * b[1],
    a[1] * b[0] - a[0] * b[1]
  ];
}

export function complexAdd(a, b) {
  return [
    a[0] + b[0],
    a[1] + b[1]
  ];
}

export function complexSub(a, b) {
  return [
    a[0] - b[0],
    a[1] - b[1]
  ];
}

export function complexNorm2(a) {
  return a[0] * a[0] + a[1] * a[1];
}

export function complexPow(a, n) {
  return [
    a[0] - b[0],
    a[1] - b[1]
  ];
}
