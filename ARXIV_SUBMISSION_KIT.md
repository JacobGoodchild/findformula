# arXiv submission kit

Everything needed to submit the papers carefully. Recommended order: submit **Paper 2 first**
(3D sandpiles). It is the most concrete and has the most independent checks. Submit Paper 1 only
after Paper 2 has been accepted and announced.

---

## 1. Before you submit (to protect your account)

- [ ] **Read the paper in full at least once**, using the "plain-English guide" in section 6. You don't
      need to re-derive everything, but you should be able to explain what each main result says and how
      it was checked.
- [ ] **Check the PDF preview arXiv generates** after upload. Look for broken equations or missing
      references; there should be none.
- [ ] **Keep the AI disclosure.** It is in the abstract, the introduction and the AI-use statement.
      Hiding AI use is what gets people in trouble. Disclosing it is the responsible route.
- [ ] **Submit one paper, not both at once.** A single solid paper looks far better to moderators than
      two long AI-assisted papers from a new account on the same day.
- [ ] **Don't resubmit a rejected paper unchanged.** If a paper is rejected or put on hold, read the
      reason, fix it or wait, and appeal politely at most once.

## 2. Form fields: Paper 2 (submit first)

**File to upload:** `paper2_sandpiles_and_lattice_entropy.tex` (a single file with its bibliography
included; arXiv compiles it with REVTeX).

**Title:**
Exact height-one probabilities of three-dimensional Abelian sandpiles, spanning-tree entropies of the snub-square and Cairo lattices, and elementary resistor-network reductions

**Authors:** Jacob Goodchild

**Primary category:** cond-mat.stat-mech   **Cross-lists:** math-ph, math.PR

**Comments:**
7 pages, 3 tables. AI use disclosed: formulas were found and the text drafted with an AI system under the author's direction; see the AI-use statement. Code: https://github.com/JacobGoodchild/findformula

**Licence:** CC BY 4.0

**Abstract** (paste as-is; arXiv's limit is 1920 characters and this one fits):

```
We give exact closed forms for the height-one probability $P_1$ of the Abelian sandpile, spanning-tree entropies and two-point resistances on a range of lattices, using the fact that all three are controlled by the same Kirchhoff transfer currents. In three dimensions, $P_1$ on the simple-cubic, body-centred cubic, face-centred cubic and pyrochlore lattices is a polynomial in the corresponding Watson integral $W$ and $1/(\pi^2 W)$; on the simple cubic lattice $P_1=\frac{4}{3}\alpha^3(2-3\alpha)^2\approx 0.0545825$ with $\alpha=(7W+54/(\pi^2 W))/36$. For bipartite lattices of girth at least six with symmetric second neighbours we prove $P_1=(z-2)^{z-1}/[z(z-1)^{z-1}]$, giving $2/27$ for diamond. We also give exact $P_1$ for ten two-dimensional lattices, the spanning-tree entropy $\frac{10}{3\pi}G+\frac{1}{3}\ln(2+\sqrt{3})$ of the snub-square (Shastry-Sutherland) lattice and of its dual Cairo lattice, a Clausen-function entropy and elementary resistances for the Shastry-Sutherland network with any dimer conductance, an elementary formula for the axis-bond resistance of any square lattice with diagonal bonds, and a one-Green-function reduction for next-nearest-neighbour networks on bipartite lattices. Closed forms were identified by integer-relation (PSLQ) searches at 25-60 digits and confirmed by exact finite-volume computations, independent numerical integration and Monte Carlo simulation. This research, including the discovery of the formulas, the computations and the writing, was carried out with extensive use of an AI system, as described in the paper.
```

## 3. Form fields: Paper 1 (later)

**File:** `paper1_quantum_walks_and_graph_laws.tex`

**Title:**
Exact localization of Grover quantum walks on crystal lattices: transfer currents, a universal bound, and uniform-spanning-tree degree laws

**Primary category:** quant-ph   **Cross-lists:** math-ph, cond-mat.stat-mech

**Comments:** 9 pages, 2 tables. AI use disclosed: see the AI-use statement. Code: https://github.com/JacobGoodchild/findformula

**Abstract:**

```
Grover quantum walks on periodic graphs localize: a finite fraction of the wave function stays at the starting vertex because the walk operator has flat eigenvalues. We compute the long-time return probability $\bar p$ exactly for many lattices. For the flip-flop walk we show that the $\pm1$ flat-band projectors are $I-Y$ and $I-Y^Q$, where $Y$ and $Y^Q$ are the Kirchhoff transfer-current matrices of the Laplacian and signless Laplacian. This reduces trapping to lattice Green functions and yields closed forms in terms of Watson integrals (simple cubic, BCC, FCC), complete elliptic integrals (triangular, kagome, checkerboard, star, line graphs) and elementary functions of the bond anisotropy. We prove the lower bound $\bar p\ge (d-2)^2/[2d(d-1)]$ for edge-transitive lattices of vertex degree $d$, attained on 2-arc-transitive bipartite lattices (honeycomb, diamond, Laves graph). We solve lazy (Szegedy) walks on the square, triangular, honeycomb and diamond lattices as functions of the laziness. For the moving-shift walk on $\mathbb{Z}^d$ the trapping is governed by $A_d=\int_0^\infty \mathrm{erfc}(\sqrt{s})^d\,ds$, with elementary closed forms for $d\le 5$. The same transfer-current matrix gives the uniform-spanning-tree degree distribution, which is exactly $1+\mathrm{Binomial}(z-1,1/(z-1))$ on 2-arc-transitive lattices. Closed forms were identified by PSLQ at 25-80 digits and checked by Brillouin-zone diagonalization and simulation. This research, including the discovery of the formulas, the computations and the writing, was carried out with extensive use of an AI system, as described in the paper.
```

## 4. Getting an endorsement without knowing a professor

1. Start the submission on arXiv and choose the category (cond-mat.stat-mech). arXiv will then tell you
   that you need an endorsement and give you a **personal endorsement code**.
2. Find people who are allowed to endorse. Open the arXiv page of any recent paper in cond-mat.stat-mech
   on sandpiles, spanning trees or lattice Green functions. At the bottom of the abstract page there is a
   link, **"Which authors of this paper are endorsers?"**, which lists authors who can endorse you.
   Good places to look are the papers cited in Paper 2, for example arXiv:2609.10352 (kagome sandpile,
   posted this month, so its authors are active) and arXiv:1708.03493.
3. Send the email below to **one person at a time**, and wait about a week before trying the next.
   Don't mass-email; that looks like spam.
4. If nobody replies after three or four tries, use the fallback: post to Zenodo (no endorsement
   needed) and share the link when you contact experts.

### Endorsement request email (edit the bracketed parts)

> **Subject:** Question about exact sandpile height probabilities in 3D (and an arXiv endorsement request)
>
> Dear Dr [Name],
>
> I read your paper [title/arXiv number] and have a short note that I believe is closely related. It
> gives closed forms for the height-one probability of the Abelian sandpile on the simple-cubic, BCC,
> FCC and pyrochlore lattices, in terms of Watson's integrals (e.g. $P_1 = \tfrac43\alpha^3(2-3\alpha)^2 \approx 0.0545825$
> on the simple cubic lattice). Each value was checked against exact finite-box computations and
> independent numerical integration.
>
> I should be upfront: I am not an academic, and the results were found and the paper drafted with
> extensive use of an AI system. This is disclosed in full in the manuscript. That is exactly why I
> would value an expert's view before anything is posted.
>
> The PDF is attached, and all code is public at https://github.com/JacobGoodchild/findformula.
>
> Would you be willing to take a quick look? If you think the results are correct and worth sharing, I
> would be grateful for an arXiv endorsement for cond-mat.stat-mech (my code is [CODE]). If they are
> already known, or you see a problem, I would equally appreciate hearing that.
>
> Thank you for your time.
>
> Kind regards,
> Jacob Goodchild
> jacob.l.goodchild@gmail.com

## 5. If something goes wrong

- **"On hold":** normal for new authors. It can take a few days. Don't resubmit.
- **Rejected or reclassified:** read the moderator's note. Reclassification to another category is fine.
  For a rejection, wait and consider Zenodo instead of appealing repeatedly.
- **Someone finds an error or shows a result is already known:** post a replacement version (v2) that
  corrects it and credits them. That is normal and respected practice, not a failure.

## 6. Plain-English guide to Paper 2 (so you can explain it)

- **What a sandpile is:** grains are dropped on a grid. When a site holds too many, it "topples" and passes
  one grain to each neighbour. After a long time the pile settles into a critical state. $P_1$ is the chance that
  a given site holds the smallest possible number of grains.
- **Main idea:** Majumdar and Dhar (1991) showed that $P_1$ is a small determinant built from electrical
  resistances between nearby points of an infinite grid of 1-ohm resistors. The new work is finding those
  resistances *exactly* in 3D, where they turn out to be simple combinations of a famous constant
  (Watson's integral, 1939). The determinant then collapses to a short formula.
- **How it was checked:** (1) the same method reproduces the known 2D answers (square, honeycomb, kagome);
  (2) exact calculations on finite boxes get closer and closer to the formula as the box grows;
  (3) simulations of actual sandpiles agree within their error bars; (4) the resistances were recomputed
  from scratch by numerical integration and match to about 15 digits (`verification/indep.py`).
- **What might not be new:** the simple-cubic case could be implicit in Majumdar and Dhar's own
  d-dimensional analysis, and the paper says so. BCC, FCC, pyrochlore and diamond ($2/27$) were not found
  anywhere in our searches.
- **The spanning-tree entropy:** counts how many ways there are to connect all points of a lattice with no
  loops. The snub-square answer mixes two classic constants, Catalan's constant and $\ln(2+\sqrt3)$. It was
  found numerically to 60 digits and then matched by an integer-relation search.
