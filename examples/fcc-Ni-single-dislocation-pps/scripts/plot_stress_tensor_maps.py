#!/usr/bin/env python3
"""Full per-atom stress tensor of the relaxed single-dislocation slabs: uncorrected box (as built by
dcreator + lego-cut) vs corrected period vector (Lx', xy').

Frame: x = glide direction (periodic), y = line direction (periodic), z = glide-plane normal (free).
Per-atom stress = LAMMPS stress/atom divided by the perfect-crystal atomic volume a0^3/4 (exact far
from the core and the surfaces, approximate in the core). Positive = tensile. GPa.

Outputs (OUTDIR):
  stress-tensor-maps_<name>.png        6 components, y-averaged x-z maps, rows = uncorrected / corrected
  stress-tensor-farfield-profiles.png  sigma_ij(z) averaged over |x - x_line| > Lx/4, all characters
  stress-tensor-farfield-summary.dat   top / bottom / mean far-field values per component
usage: plot_stress_tensor_maps.py CELLDIR RELAXDIR OUTDIR NAME [NAME ...] [--a0 3.52]
"""
import argparse, os, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

COMP = ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']
# diverging: blue (compressive) - neutral gray - red (tensile); symmetric limits always
CMAP = LinearSegmentedColormap.from_list('div', ['#184f95', '#6da7ec', '#f0efec', '#ec7a6e', '#a8201a'])
C_UNCORR, C_CORR = '#eb6834', '#2a78d6'          # categorical slots 2 and 1
INK, INK2 = '#0b0b0b', '#52514e'


def read_vars(f):
    v = {}
    for l in open(f):
        m = re.match(r'variable\s+(\S+)\s+equal\s+(\S+)', l)
        if m:
            v[m.group(1)] = float(m.group(2))
    return v


def read_dump(f, rotated):
    """Return per-atom data in the (glide x, line y, normal z) frame.
    rotated=True: the corrected cells were relaxed in the LAMMPS frame X = line, Y = -glide (see
    apply_period_vector.py); positions and stresses are mapped back here."""
    with open(f) as fh:
        L = fh.read().split('\n')
    i = L.index(next(l for l in L if l.startswith('ITEM: BOX BOUNDS')))
    tilt = 'xy' in L[i]
    bx = [list(map(float, L[i + k].split())) for k in (1, 2, 3)]
    xy = bx[0][2] if tilt else 0.0
    xlo = bx[0][0] - min(0.0, xy) if tilt else bx[0][0]
    xhi = bx[0][1] - max(0.0, xy) if tilt else bx[0][1]
    ylo, yhi = bx[1][:2]
    cols = L[i + 4].split()[2:]
    a = np.loadtxt(L[i + 5:])
    raw = {c: a[:, k] for k, c in enumerate(cols)}
    S = [raw[f'c_s[{k}]'] for k in range(1, 7)]          # xx yy zz xy xz yz in the dump frame
    d = {'z': raw['z'], 'c_cna': raw['c_cna']}
    if rotated:
        Lg = yhi - ylo
        d['xs'] = np.mod(Lg - (raw['y'] - ylo), Lg)
        d['_s'] = [S[1], S[0], S[2], -S[3], -S[5], S[4]]
    else:
        Lg = xhi - xlo
        d['xs'] = np.mod(raw['x'] - xlo, Lg)
        d['_s'] = S
    return d, Lg, xy


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('celldir'); ap.add_argument('relaxdir'); ap.add_argument('outdir')
    ap.add_argument('names', nargs='+'); ap.add_argument('--a0', type=float, default=3.52)
    ap.add_argument('--bin', type=float, default=2.5)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    omega = a.a0 ** 3 / 4.0
    plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                         'xtick.color': INK2, 'ytick.color': INK2, 'axes.linewidth': 0.6})
    summary = []
    prof = {}
    for name in a.names:
        v = read_vars(os.path.join(a.celldir, f'box-correction_{name}.lmp'))
        data = {}
        for var in ('uncorrected', 'corrected'):
            d, Lx, xy = read_dump(os.path.join(a.relaxdir, f'relaxed_{name}_{var}_stress-tensor.dump'),
                                  rotated=(var == 'corrected'))
            for k, c in enumerate(COMP):
                d['s' + c] = d['_s'][k] / omega / 1e4
            data[var] = (d, Lx, xy)
        zg, xl_built = v['z_glide'], v['x_line']
        zlo_i, zhi_i = v['z_fix_lo'] + 4, v['z_fix_hi'] - 4      # interior, away from the 2D layers

        # ---- maps
        fig, axs = plt.subplots(2, 6, figsize=(17, 6.4), constrained_layout=True)
        maps = {}
        for var in ('uncorrected', 'corrected'):
            d, Lx, xy = data[var]
            xl = xl_built * Lx / v['Lx_built'] if var == 'corrected' else xl_built
            # z bins = one (111) plane each (bins straddling planes alias into stripes); x bins ~a.bin
            d111 = a.a0 / np.sqrt(3.0)
            ez = np.arange(d['z'].min() - 0.5 * d111, d['z'].max() + d111, d111)
            ex = np.linspace(0, Lx, int(round(Lx / a.bin)) + 1)
            n, _, _ = np.histogram2d(d['xs'], d['z'], [ex, ez])
            hcp, _, _ = np.histogram2d(d['xs'], d['z'], [ex, ez], weights=(d['c_cna'] == 2).astype(float))
            for c in COMP:
                s, _, _ = np.histogram2d(d['xs'], d['z'], [ex, ez], weights=d['s' + c])
                with np.errstate(invalid='ignore', divide='ignore'):
                    maps[(var, c)] = (s / n, hcp / n, ex, ez, xl, Lx)
        for j, c in enumerate(COMP):
            # symmetric limit from the interior far from the core, shared by both rows
            vals = []
            for var in ('uncorrected', 'corrected'):
                m, _, ex, ez, xl, Lx = maps[(var, c)]
                X, Z = np.meshgrid(0.5 * (ex[1:] + ex[:-1]), 0.5 * (ez[1:] + ez[:-1]), indexing='ij')
                dx = np.abs(X - xl); dx = np.minimum(dx, Lx - dx)
                sel = (np.hypot(dx, Z - zg) > 25) & (Z > zlo_i) & (Z < zhi_i) & np.isfinite(m)
                vals.append(np.abs(m[sel]))
            lim = max(np.percentile(np.concatenate(vals), 99), 0.05)
            for i, var in enumerate(('uncorrected', 'corrected')):
                m, h, ex, ez, xl, Lx = maps[(var, c)]
                ax = axs[i, j]
                im = ax.pcolormesh(ex, ez, m.T, cmap=CMAP, vmin=-lim, vmax=lim, shading='flat', rasterized=True)
                ax.contour(0.5 * (ex[1:] + ex[:-1]), 0.5 * (ez[1:] + ez[:-1]), np.nan_to_num(h).T, [0.3],
                           colors=INK, linewidths=0.6)
                for zz in (v['z_fix_lo'], v['z_fix_hi']):
                    ax.axhline(zz, color=INK2, lw=0.5, ls=':')
                ax.set_aspect('equal'); ax.set_xlim(0, Lx)
                ax.set_title(f'$\\sigma_{{{c}}}$ — {var}', fontsize=9, color=INK)
                if j == 0: ax.set_ylabel('z (Å)')
                if i == 1: ax.set_xlabel('x (Å)')
                if i == 1:
                    cb = fig.colorbar(im, ax=axs[:, j], orientation='horizontal', shrink=0.85, pad=0.02)
                    cb.set_label(f'$\\sigma_{{{c}}}$ (GPa), clipped at ±{lim:.2f}', fontsize=8)
        fig.suptitle(f'{name}: y-averaged per-atom stress; black contour = stacking fault (hcp atoms); '
                     f'dotted = 2D-dynamic layers. Top row: box as built; bottom row: Lx\' = {v["Lx_corr"]:.4f} Å, '
                     f'δ = {v["xy_corr"]:+.4f} Å (period vector (Lx\', δ, 0))', fontsize=9.5, color=INK)
        fig.savefig(os.path.join(a.outdir, f'stress-tensor-maps_{name}.png'), dpi=150)
        plt.close(fig)

        # ---- far field
        for var in ('uncorrected', 'corrected'):
            d, Lx, xy = data[var]
            xl = xl_built * Lx / v['Lx_built'] if var == 'corrected' else xl_built
            dx = np.abs(d['xs'] - xl); dx = np.minimum(dx, Lx - dx)
            far = dx > Lx / 4
            top = far & (d['z'] > zg + 15) & (d['z'] < zhi_i)
            bot = far & (d['z'] < zg - 15) & (d['z'] > zlo_i)
            row = [name, var]
            for c in COMP:
                t, b = d['s' + c][top].mean(), d['s' + c][bot].mean()
                row += [t, b, 0.5 * (t + b)]
            summary.append(row)
            ez = np.arange(d['z'].min() - 0.1, d['z'].max() + 2.04, 2.0322)
            zc = 0.5 * (ez[1:] + ez[:-1]); n, _ = np.histogram(d['z'][far], ez)
            prof[(name, var)] = (zc, {c: np.histogram(d['z'][far], ez, weights=d['s' + c][far])[0] / np.maximum(n, 1)
                                      for c in COMP}, zg, v)

    # ---- profile figure: rows = characters, cols = components
    fig, axs = plt.subplots(len(a.names), 6, figsize=(16, 3.0 * len(a.names)), sharey='row', constrained_layout=True)
    axs = np.atleast_2d(axs)
    for i, name in enumerate(a.names):
        for j, c in enumerate(COMP):
            ax = axs[i, j]
            ax.axvline(0, color=INK2, lw=0.5)
            for var, col in (('uncorrected', C_UNCORR), ('corrected', C_CORR)):
                zc, p, zg, v = prof[(name, var)]
                inside = (zc > v['z_fix_lo']) & (zc < v['z_fix_hi'])
                ax.plot(p[c][inside], zc[inside], color=col, lw=2, label=var)
            ax.axhline(zg, color=INK2, lw=0.6, ls='--')
            ax.grid(color='#e4e3df', lw=0.5)
            if i == 0: ax.set_title(f'$\\sigma_{{{c}}}$ (GPa)', color=INK)
            if j == 0: ax.set_ylabel(f'{name}\nz (Å)', fontsize=8.5)
    axs[0, 0].legend(frameon=False, fontsize=8, loc='lower left')
    fig.suptitle('Far-field stress profiles, averaged over |x − x_line| > Lx/4 (dashed: glide plane). '
                 'Corrected box: top and bottom halves equal and opposite, zero mean.', color=INK)
    fig.savefig(os.path.join(a.outdir, 'stress-tensor-farfield-profiles.png'), dpi=150)
    plt.close(fig)

    with open(os.path.join(a.outdir, 'stress-tensor-farfield-summary.dat'), 'w') as fh:
        fh.write('# far-field mean per-atom stress (GPa), |x - x_line| > Lx/4, 15 A off the glide plane, '
                 'outside the 2D layers + 4 A\n# name variant ' +
                 ' '.join(f's{c}_top s{c}_bot s{c}_mean' for c in COMP) + '\n')
        for r in summary:
            fh.write(f'{r[0]} {r[1]} ' + ' '.join(f'{x:+.4f}' for x in r[2:]) + '\n')
    print(open(os.path.join(a.outdir, 'stress-tensor-farfield-summary.dat')).read())


if __name__ == '__main__':
    main()
