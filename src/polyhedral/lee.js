import {geoProjection as projection, geoStream} from "../d3-geo";
import polyhedral from "./index";
import {scan} from "../d3-array";
import {asin, degrees, epsilon, sqrt} from "../math";
import {lagrangeRaw} from "../lagrange";
import {complexAdd, complexMul, complexNorm2, complexPow} from "../complex";


var leeRaw = function(lambda, phi) {
  // return d3.geoGnomonicRaw(...arguments);
  var w = Complex([-1/2, Math.sqrt(3)/2]),
      k = Complex(0),
      h = Complex(0),
      z = Complex(d3.geoStereographicRaw(lambda, phi)).mul(Math.sqrt(2));

  // rotate to have s ~= 1
  var rot = w.clone().pow(d3.scan([0,1,2].map(
    i => -(z.clone().mul(w.clone().pow(i))).re
  )));

  var n = z.abs();

  if (n > 0.3) {


    // if |z| > 0.5, use the approx based on y = (1-z)
    // McIlroy formula 6 p6 and table for G page 16
    var y = rot.clone().mul(z).mul(-1).add(1);

    // w1 = gamma(1/3) * gamma(1/2) / 3 / gamma(5/6);
    // https://bl.ocks.org/Fil/1aeff1cfda7188e9fbf037d8e466c95c
    var w1 = 1.4021821053254548;

    var G0 = [
      1.15470053837925,
      0.192450089729875,
      0.0481125224324687,
      0.010309826235529,
      3.34114739114366e-4,
      -1.50351632601465e-3,
      -1.23044177962310e-3,
      -6.75190201960282e-4,
      -2.84084537293856e-4,
      -8.21205120500051e-5,
      -1.59257630018706e-6,
      1.91691805888369e-5,
      1.73095888028726e-5,
      1.03865580818367e-5,
      4.70614523937179e-6,
      1.4413500104181e-6,
      1.92757960170179e-8,
      -3.82869799649063e-7,
      -3.57526015225576e-7,
      -2.2175964844211e-7
    ];

    var G = Complex(0);
    for (var i = G0.length; i--;) {
      G = Complex(G0[i]).add(G.mul(y));
    }

    k = Complex(w1).add(y.sqrt().mul(-1).mul(G)).mul(rot).mul(rot)

  }

  if (n < 0.5) {

    // if |z| < 0.3
    // https://www.wolframalpha.com/input/?i=series+of+((1-z%5E3))+%5E+(-1%2F2)+at+z%3D0 (and ask for "more terms")
    // 1 + z^3/2 + (3 z^6)/8 + (5 z^9)/16 + (35 z^12)/128 + (63 z^15)/256 + (231 z^18)/1024 + O(z^21)
    // https://www.wolframalpha.com/input/?i=integral+of+1+%2B+z%5E3%2F2+%2B+(3+z%5E6)%2F8+%2B+(5+z%5E9)%2F16+%2B+(35+z%5E12)%2F128+%2B+(63+z%5E15)%2F256+%2B+(231+z%5E18)%2F1024
    // (231 z^19)/19456 + (63 z^16)/4096 + (35 z^13)/1664 + z^10/32 + (3 z^7)/56 + z^4/8 + z + constant
    var H0 = [
      1, 1/8, 3/56, 1/32, 35/1664, 63/4096, 231/19456
    ]
    var z3 = z.clone().pow(3);
    for (var i = H0.length; i--;) {
      h = Complex(H0[i]).add(h.mul(z3));
    }
    h = h.mul(z);
  }


  if (n < 0.3) return h.toVector();
  if (n > 0.5) return k.toVector();

  // in between 0.3 and 0.5, interpolate
  var t = (n - 0.3) / (0.5 - 0.3);
  return k.mul(t).add(h.mul(1 - t)).toVector();

}

// w1 = gamma(1/n) * gamma(1 - 2/n) / n / gamma(1 - 1/n)
// https://bl.ocks.org/Fil/852557838117687bbd985e4b38ff77d4
var w = [-1/2, sqrt(3)/2],
    w1 = [1.7666387502854533, 0],
    m = 0.3 * 0.3;

// Approximate \int _0 ^sm(z)  dt / (1 - t^3)^(2/3)
// sm maps a triangle to a disc, sm^-1 does the opposite
function sm_1(z) {
    var k = [0, 0];

    // rotate to have s ~= 1
    var rot = complexPow(w, scan([0, 1, 2].map(function(i) {
      return -complexMul(z, complexPow(w, [i, 0]))[0];
    })));
      
    var y = complexMul(rot, z);
    y = [1 - y[0], - y[1]];
    
    // McIlroy formula 5 p6 and table for F3 page 16
    var F0 = [
      1.44224957030741,
      0.240374928384568,
      0.0686785509670194,
      0.0178055502507087,
      0.00228276285265497,
      -1.48379585422573e-3,
      -1.64287728109203e-3,
      -1.02583417082273e-3,
      -4.83607537673571e-4,
      -1.67030822094781e-4,
      -2.45024395166263e-5,
      2.14092375450951e-5,
      2.55897270486771e-5,
      1.73086854400834e-5,
      8.72756299984649e-6,
      3.18304486798473e-6,
      4.79323894565283e-7
      -4.58968389565456e-7,
      -5.62970586787826e-7,
      -3.92135372833465e-7
    ];
    
    var F = [0, 0];
    for (var i = F0.length; i--;) F = complexAdd([F0[i],0], complexMul(F, y));

    k = complexMul(
      complexAdd(w1,
        complexMul([-F[0], -F[1]], complexPow(y, (1-2/3)))
      ),
      complexMul(rot, rot)
    );

    // when we are close to [0,0] we switch to another approximation:
    // https://www.wolframalpha.com/input/?i=(-2%2F3+choose+k)++*+(-1)%5Ek++%2F+(k%2B1)+with+k%3D0,1,2,3,4
    // the difference is _very_ tiny but necessary
    // if we want projection(0,0) === [0,0]
    var n = complexNorm2(z);
    if (n < m) {
      var H0 = [
        1, 1/3, 5/27, 10/81, 22/243 //…
       ];
      var z3 = complexPow(z, [3,0]);
      var h = [0,0];
      for (i = H0.length; i--;) h = complexAdd([H0[i],0], complexMul(h, z3));
      h = complexMul(h, z);
      k = complexAdd(complexMul(k, [n / m, 0]), complexMul(h, [1 - n / m, 0]));
    }
    
    return k;
  }


var lagrange1_2 = lagrangeRaw(0.5);
export function coxRaw(lambda, phi) {
  var s = lagrange1_2(lambda, phi);
  var t = sm_1([s[1] / 2, s[0] / 2]);
  return [t[1], t[0]]
}

// the Sphere should go *exactly* to the vertices of the triangles
// because they are singular points
function sphere() {
  var c = 2 * asin(1 / sqrt(5)) * degrees;
  return {
    type: "Polygon",
    coordinates: [
      [ [ 0,90 ],  [ -180, -c + epsilon ], [ 0, -90 ], [ 180, -c + epsilon ], [ 0,90 ] ]
    ]
  };
}

export default function() {
  var p = projection(coxRaw);

  var stream_ = p.stream;
  p.stream = function(stream) {
    var rotate = p.rotate(),
        rotateStream = stream_(stream),
        sphereStream = (p.rotate([0, 0]), stream_(stream));
    p.rotate(rotate);
    rotateStream.sphere = function() { geoStream(sphere(), sphereStream); };
    return rotateStream;
  };

  return p
    .scale(188.682)
    .center([0, 0])
    .translate([ 480, 500 * 2 /3 ]);
}

