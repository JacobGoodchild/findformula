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

## Formula 19 (+1 channel, moving-shift Grover walk on the triangular lattice): structure found, no closed form
- Denominator 3 + 2(c1+c2+c3) + (c1c2+c2c3+c3c1). Its z2-quartic factors as z1 (z2^2 + s+ z2 + z1)(z2^2 + s- z2 + z1),
  s+- = (z1+1)[(u+4) +- sqrt((u+4)^2 - 20)]/2, u = 2cos k1.
- Rationalising with u + 4 = sqrt5 (v + 1/v) and w = sqrt5 v, each branch has discriminant proportional to
  (w - 1)(w^2 - w + 4): an elliptic curve, modulus m = 3/8 (or 5/8).
- The integration path runs from the branch point w = 1 to w = 1 +- 2i (at u = -2), and (1 + 2i, 2i) is NOT a torsion
  point (elliptic-log coordinates 0.0909379640545..., 0.3181240718909...). So the amplitudes involve INCOMPLETE elliptic
  integrals; PSLQ with complete K(3/8), E(3/8) (several algebraic scalings) finds nothing, as expected.

## Genus-0 diagonal decorations of the square lattice (2x2 period) — classification (lattices/plaquette_search.py)
Each plaquette of a 2x2 block gets nothing '.', '/', '\' or 'X'. Up to symmetry, the decorations whose Laplacian
spectral curve is genus 0 after removing the node (hence elementary resistances, Clausen-type tree entropies) are:
  ....(square)  ...X  /.\.  ..XX  /..\(snub square, F26/27)  /..X  .XX.(checkerboard)  /XX.  .XXX  ////(triangular)
  //\\  /\\/(Union Jack)  /X\X  /XX/  /XX\  XXXX(king, F28)
E ln det per 4-site cell (45 digits computed; 30 shown), closed forms where found:
  ....  4.66497446449310048221415130349  = 16G/pi
  ...X  5.63164375130332902671634143782   (odd factor z^2 - 10z + 1, rho = (sqrt2+sqrt3)^2)
  /.\.  5.58454671667953564832521801423   (z^2 - 82z + 1)
  ..XX  6.44942187292124365217700140508   (z^2 - 7z + 1, rho = phi^4)
  /..\  5.64342258297733934667852121599  = 40G/(3pi) + (4/3) ln(2+sqrt3)
  /..X  6.06983523543591051694249502142   (z^2 - 18z + 1, rho = (2+sqrt5)^2)
  .XX.  6.49137031560622209761046838050  = 8G/pi + 6 ln2   (line graph of Z^2)
  /XX.  6.83334058516703239962365131035
  .XXX  7.17185152094717409212058941072
  ////  6.46131894438901028187273021448  = 4 z_tri = 20 Cl2(pi/3)/pi  (my first basis wrongly had sqrt3 Cl2(pi/3))
  //\\  6.46131894438901028187273021448  (same as triangular; z^2 - 194z + 1)
  /\\/  6.29347420303065743297326391596  (= Union Jack = SS(w=2), F27)
  /X\X  7.13002116816469326074473897393
  /XX/  7.13812264332348844151004533034
  /XX\  7.17430657185560998476862293732
  XXXX  7.77495744082183413702358550582  = 80G/(3pi)
The unidentified ones need their own angle (as Phi in F27 / arccos(2w) in F28); a weighted-family fit per lattice
(as done for SS and king) should close each of them.
