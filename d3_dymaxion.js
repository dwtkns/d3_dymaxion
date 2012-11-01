(function() {
  var ε = 1e-6,
      π = Math.PI,
      sqrtπ = Math.sqrt(π);

  var robinsonConstants = [
    [1.0000, 0.0000],
    [0.9986, 0.0620],
    [0.9954, 0.1240],
    [0.9900, 0.1860],
    [0.9822, 0.2480],
    [0.9730, 0.3100],
    [0.9600, 0.3720],
    [0.9427, 0.4340],
    [0.9216, 0.4958],
    [0.8962, 0.5571],
    [0.8679, 0.6176],
    [0.8350, 0.6769],
    [0.7986, 0.7346],
    [0.7597, 0.7903],
    [0.7186, 0.8435],
    [0.6732, 0.8936],
    [0.6213, 0.9394],
    [0.5722, 0.9761],
    [0.5322, 1.0000]
  ];

  function sinci(x) {
    return x ? x / Math.sin(x) : 1;
  }

  function sgn(x) {
    return x > 0 ? 1 : x < 0 ? -1 : 0;
  }

  function asqrt(x) {
    return x > 0 ? Math.sqrt(x) : 0;
  }

  // Calculate F(φ+iψ|m).
  // See Abramowitz and Stegun, 17.4.11.
  function ellipticFi(φ, ψ, m) {
    var r = Math.abs(φ),
        i = Math.abs(ψ),
        sinhψ = .5 * ((sinhψ = Math.exp(i)) - 1 / sinhψ);
    if (r) {
      var cscφ = 1 / Math.sin(r),
          cotφ2 = (cotφ2 = Math.cos(r) * cscφ) * cotφ2,
          b = -(cotφ2 + m * (sinhψ * sinhψ * cscφ * cscφ + 1) - 1),
          cotλ2 = .5 * (-b + Math.sqrt(b * b - 4 * (m - 1) * cotφ2));
      return [
        ellipticF(Math.atan(1 / Math.sqrt(cotλ2)), m) * sgn(φ),
        ellipticF(Math.atan(Math.sqrt(Math.max(0, cotλ2 / cotφ2 - 1) / m)), 1 - m) * sgn(ψ)
      ];
    }
    return [
      0,
      ellipticF(Math.atan(sinhψ), 1 - m) * sgn(ψ)
    ];
  }

  // Calculate F(φ|m) where m = k² = sin²α.
  // See Abramowitz and Stegun, 17.6.7.
  function ellipticF(φ, m) {
    var a = 1,
        b = Math.sqrt(1 - m),
        c = Math.sqrt(m);
    for (var i = 0; Math.abs(c) > ε; i++) {
      if (φ % π) {
        var dφ = Math.atan(b * Math.tan(φ) / a);
        if (dφ < 0) dφ += π;
        φ += dφ + ~~(φ / π) * π;
      } else φ += φ;
      c = (a + b) / 2;
      b = Math.sqrt(a * b);
      c = ((a = c) - b) / 2;
    }
    return φ / (Math.pow(2, i) * a);
  }

//---  


function dymaxion(λ, φ) {
    // var α = Math.sqrt(4 - 3 * Math.sin(Math.abs(φ)));
    // return [
    //   2 / Math.sqrt(6 * π) * λ * α,
    //   sgn(φ) * Math.sqrt(2 * π / 3) * (2 - α)
    // ];
    var cartesian = convert_s_t_p(λ,φ);
    return [
    cartesian.x
    ,
    cartesian.y
    ]
  }


//---


  var projection = d3.geo.projection,
      projectionMutator = d3.geo.projectionMutator;

  d3.geo.dymaxion = function() { return projection(dymaxion); };

})();
