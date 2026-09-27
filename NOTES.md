# Research log

## Negative / known results (not in formulas.txt)
- Spanning trees of lattice regions: no new "round" families; square/triangle relation is known (Ciucu–Yan–Zhang, Knuth).
- Linear extensions of triangular-lattice bands: equal shifted Young diagrams (Thrall's formula), so known.
- Faxén/Bohlin pipe wall-correction coefficient, computed to 60 digits:
  k = 2.10444309167471731461161075194240325943115689225704842587604
  No integer relation found with 1, pi, pi^2, 1/pi, zeta(3), Catalan, log 2, etc. (coeffs < 1e5).
- Stimson–Jeffery drag factor for two touching equal spheres = (1/3) ∫_0^∞ [1 - 2(sinh²x - x²)/(sinh 2x + 2x)] dx
  = 0.645141425512754913652184520491451831947959684397585772191862 ; no closed form found.
- Flip-flop Grover trapping, no closed form found (PSLQ basis 1, sqrt3/pi, 1/pi, 1/pi^2):
  ruby 3.4.6.4: 0.16458341395062199972, 0.17343348160423351562;
  snub hexagonal 3^4.6: 0.23328925868537919944, 0.23188836982107697332, 0.22198045285555102298,
  0.22157633698807354436, 0.22802366367916230407.
- King's lattice flip-flop trapping: 0.32326601977799810277 (start along an axis); resistances
  not in simple bases (possibly elliptic).
