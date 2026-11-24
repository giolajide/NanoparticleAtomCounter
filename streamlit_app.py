"""NanoparticleAtomCounter web interface. Uses the existing CLI unchanged."""
from __future__ import annotations

import ast
import hashlib
import io
import json
import logging
import math
import subprocess
import sys
import tempfile
from pathlib import Path
from urllib.request import Request, urlopen

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st
from NanoparticleAtomCounter.count_utilities import DEFAULT_FACETS
from NanoparticleAtomCounter.data import (
    atomic_numbers, reference_states, molar_volumes, interplanar_dist,
)

ROOT = Path(__file__).resolve().parent
DOCS = ROOT / 'docs'
REPO = 'https://github.com/giolajide/NanoparticleAtomCounter'
PAPER = 'https://doi.org/10.26434/chemrxiv.15001082/v1'
PYPI = 'https://pypi.org/project/NanoparticleAtomCounter/'
CITATION = ('Olajide, G.; Szilvási, T. NanoparticleAtomCounter: A Python/Web-based '
            'package to convert nanoparticle geometry from Transmission Electron '
            'Microscopy into atom counts. ChemRxiv (preprint). '
            'https://doi.org/10.26434/chemrxiv.15001082/v1')
HEADERS = ['r (A)', 'R (A)', 'Theta', 'Element', 'Interface Facet', 'Surface Facet']
COUNTS = ['Perimeter', 'Interface', 'Surface', 'Total']
SITES = ['Bulk', 'Surface', 'Interface', 'Perimeter']
COLORS = ['#505762', '#347bb8', '#d49a17', '#c83b3e']
LABELS = ['Bulk', 'Surface', 'Interfacial', 'Perimeter']
MAX_ROWS = 1_000_000
EXAMPLE = dict(element='Pt', radius_type='Footprint radius (r)', unit='Å',
               radius=25.0, theta=109.0, interface='Default (unknown)',
               surface='Default (unknown)', mode='volume')
LOG = logging.getLogger(__name__)


def supported_elements():
    """Use the package data, not a hard-coded list of metals."""
    result = []
    for el in sorted(set(molar_volumes) & set(interplanar_dist)):
        try:
            if float(molar_volumes[el]) > 0 and any(float(d)>0 for d in interplanar_dist[el].values()):
                result.append(el)
        except (KeyError, TypeError, ValueError, IndexError):
            continue
    return result


ELEMENTS = supported_elements()


def default_facet(element):
    try:
        facet = DEFAULT_FACETS[reference_states[atomic_numbers[element]]['symmetry']]
        return facet if facet in interplanar_dist[element] else None
    except (KeyError, TypeError, IndexError):
        return None


def facet_text(value):
    if value is None or (isinstance(value, float) and math.isnan(value)) or not str(value).strip():
        return ''
    try:
        v = ast.literal_eval(str(value).strip())
        if not isinstance(v, (list, tuple)) or len(v) != 3 or any(type(x) is not int for x in v) or not any(v):
            raise ValueError
        return str(tuple(v))
    except (ValueError, SyntaxError, TypeError):
        raise ValueError('Facets must be three integer Miller indices, e.g. (1, 1, 1), or blank.') from None


def demo_data():
    return pd.DataFrame([(r, np.nan, t, 'Pt', '', '') for r,t in
                         [(15,80),(18,95),(20,109),(22,105),(25,109),(27,120),
                          (30,95),(32,110),(35,130),(40,100),(45,115),(50,125)]], columns=HEADERS)


def sample_bytes():
    p = DOCS / 'sample_input.csv'
    return p.read_bytes() if p.is_file() else demo_data().to_csv(index=False).encode()


def read_input(raw, suffix):
    if suffix == '.csv':
        return pd.read_csv(io.BytesIO(raw))
    if suffix in ('.xlsx', '.xls'):
        return pd.read_excel(io.BytesIO(raw), sheet_name=0)
    raise ValueError('Choose a CSV, XLS or XLSX file.')


def normalize_input(frame):
    """Validate and produce the exact six-column CLI schema, resolving r per row."""
    df = frame.copy().reset_index(drop=True)
    df.columns = [str(c).strip() for c in df.columns]
    if df.columns.duplicated().any():
        raise ValueError('Duplicate column headers are not allowed.')
    if set(df.columns) - set(HEADERS):
        raise ValueError('Unexpected columns: ' + ', '.join(sorted(set(df.columns) - set(HEADERS))))
    if not {'Theta','Element'} <= set(df.columns) or not {'r (A)','R (A)'} & set(df.columns):
        raise ValueError('Required columns: Theta, Element and at least one of r (A) or R (A).')
    if df.empty or len(df) > MAX_ROWS:
        raise ValueError(f'Provide between 1 and {MAX_ROWS:,} data rows.')
    for col in HEADERS:
        if col not in df:
            df[col] = '' if 'Facet' in col else np.nan
    df = df[HEADERS]
    def reject(mask, message):
        mask = np.asarray(mask, dtype=bool)
        if mask.any():
            rows = ', '.join(map(str, np.flatnonzero(mask)[:8]+1))
            raise ValueError(f'{message} Data row(s): {rows}.')
    numeric = {}
    for c in ('Theta','r (A)','R (A)'):
        orig = df[c]
        missing = orig.isna() | orig.astype(str).str.strip().eq('')
        numeric[c] = pd.to_numeric(orig, errors='coerce').astype(float)
        if c == 'Theta':
            reject(~np.isfinite(numeric[c]) | (numeric[c] <= 0) | (numeric[c] >= 180), 'Contact angle must be strictly between 0° and 180°.')
        elif c == 'r (A)':
            reject(~missing & ~np.isfinite(numeric[c]), 'Footprint radius is not a finite number.')
    r = numeric['r (A)']; R = numeric['R (A)']; theta = numeric['Theta']
    use_R = r.isna() | r.eq(0)
    reject(~use_R & (r <= 0), 'Footprint radius must be positive.')
    reject(use_R & (~np.isfinite(R) | (R <= 0)), 'A positive curvature radius is required where footprint radius is absent or zero.')
    resolved = r.where(~use_R, R * np.sin(np.deg2rad(theta)))
    reject(~np.isfinite(resolved) | (resolved <= 0), 'Resolved footprint radius must be finite and positive.')
    df['r (A)'] = resolved
    # A supplied R is ignored if r is present, as in the original CLI.
    df['R (A)'] = np.nan
    df['Theta'] = theta
    df['Element'] = df['Element'].fillna('').astype(str).str.strip().str.capitalize()
    reject(~df['Element'].isin(ELEMENTS), 'Unsupported or missing element.')
    for c in ('Interface Facet','Surface Facet'):
        # Parse each distinct string once; avoids per-particle Python parsing overhead.
        raw = df[c].fillna('').astype(str).str.strip()
        mapping = {}
        for value in raw.unique():
            try:
                mapping[value] = facet_text(value)
            except ValueError as err:
                reject(raw.eq(value), f'{c}: {err}')
        df[c] = raw.map(mapping)
        for el, facet in df[['Element',c]].drop_duplicates().itertuples(index=False, name=None):
            if not facet and default_facet(el) is None:
                reject(df['Element'].eq(el) & df[c].eq(facet), f'{c} is required for {el}; no default is available.')
            if facet and facet not in interplanar_dist[el]:
                reject(df['Element'].eq(el) & df[c].eq(facet), f'{c} {facet} is not tabulated for {el}.')
    return df


def run_counter(df, mode):
    """Run the existing CLI in an isolated temporary directory."""
    if mode not in ('volume','area'):
        raise ValueError('Unknown counting mode.')
    with tempfile.TemporaryDirectory(prefix='nac_') as tmp:
        in_path, out_path = Path(tmp)/'input.csv', Path(tmp)/'output.csv'
        df.to_csv(in_path, index=False)
        proc = subprocess.run([sys.executable,'-m','NanoparticleAtomCounter.cli.atom_count',
                               '--input',str(in_path),'--output',str(out_path),'--mode',mode],
                              capture_output=True,text=True,check=False,timeout=600)
        if proc.returncode:
            # No raw input, local filesystem paths or tracebacks shown to other services.
            detail = (proc.stderr.strip().splitlines() or ['No error detail returned.'])[-1]
            raise RuntimeError('Calculation engine failed: ' + detail[:500])
        result = pd.read_csv(out_path)
    if not set(COUNTS) <= set(result) or len(result) != len(df):
        raise RuntimeError('Unexpected calculation output schema or row count.')
    return result[COUNTS].apply(pd.to_numeric, errors='raise')


def enrich_result(raw, original, mode):
    """Retain raw counts. Bulk is a residual; never clip a nonphysical result."""
    out = raw.copy()
    out['Bulk'] = out['Total'] - out[['Surface','Interface','Perimeter']].sum(axis=1, min_count=3)
    valid = np.isfinite(out[COUNTS+['Bulk']]).all(axis=1) & out[COUNTS+['Bulk']].ge(0).all(axis=1) & out['Total'].gt(0)
    out['Physical partition'] = valid
    for site in SITES:
        out[site+' fraction (%)'] = (100*out[site]/out['Total'].where(valid)).where(valid)
    meta = original[HEADERS].copy().reset_index(drop=True)
    meta.insert(0,'Particle',np.arange(1,len(meta)+1))
    meta['Counting mode'] = mode
    return pd.concat([meta,out.reset_index(drop=True)],axis=1)


def analytics_config():
    try:
        cfg = st.secrets['analytics']
        if cfg.get('enabled',False) and all(cfg.get(k) for k in ('url','key','token')):
            return dict(cfg)
    except (KeyError, FileNotFoundError):
        pass
    except Exception:
        LOG.warning('Analytics configuration unavailable.')
    return None


def is_test_mode():
    """Skip aggregate analytics for sessions opened with ?test=1.

    This is not authentication and grants no privileged access.
    """
    return st.query_params.get('test') == '1'


def track(event, workflow, source, particles=0):
    """Only deliberate interactions. No file contents, names, IPs or user IDs."""
    if is_test_mode():
        return
    cfg = analytics_config()
    if not cfg:
        return
    payload = dict(p_token=cfg['token'], p_event=event, p_workflow=workflow,
                   p_source=source,p_particles=int(particles))
    try:
        # Legacy anon JWT needs Authorization; newer publishable keys use apikey only.
        headers = {'apikey':cfg['key'],'Content-Type':'application/json'}
        if cfg['key'].startswith('eyJ'):
            headers['Authorization'] = 'Bearer '+cfg['key']
        req = Request(cfg['url'].rstrip('/')+'/rest/v1/rpc/nac_record_event',
                      data=json.dumps(payload).encode(),method='POST',headers=headers)
        with urlopen(req, timeout=3) as response:
            response.read()
    except Exception:
        LOG.warning('Usage analytics write failed; calculation unaffected.')


def geometry_svg(radius_A, theta, radius_type):
    """True spherical-cap cross-section, not an atomistic reconstruction."""
    t = math.radians(theta)
    R = radius_A / math.sin(t) if radius_type.startswith('Footprint') else radius_A
    r = R*math.sin(t); h = R*(1-math.cos(t)); zc = -R*math.cos(t)
    # Normalize before plotting to avoid extremes; scale uniformly in x/y.
    wide = R if theta >= 90 else r
    scale = min(180/wide, 210/h)
    cx, base = 250, 258
    px = lambda x: cx+x*scale
    py = lambda z: base-z*scale
    xz = [(R*math.sin(a),zc+R*math.cos(a)) for a in np.linspace(-t,t,151)]
    points = ' '.join(f'{px(x):.3f},{py(z):.3f}' for x,z in xz)
    # Direction of tangent into particle at right edge: (-cos theta, +sin theta).
    ax, ay = px(r), base
    length = 39
    tx,ty = ax-length*math.cos(t), ay-length*math.sin(t)
    ar = 24
    arc = ' '.join(f'{ax-ar*math.cos(a):.3f},{ay-ar*math.sin(a):.3f}' for a in np.linspace(0,t,40))
    center_y = py(zc)
    radius_line = (f'<circle cx="{cx}" cy="{center_y}" r="3" fill="#ae4646"/>'
                   f'<line x1="{cx}" y1="{center_y}" x2="{ax}" y2="{base}" stroke="#ae4646" stroke-width="2" stroke-dasharray="5 3"/>') if 30 <= center_y <= 282 else ''
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 500 330">
    <rect x="1" y="1" width="498" height="328" rx="12" fill="#f1f5f9"/>
    <polygon points="{points}" fill="#c7def3" stroke="#28679d" stroke-width="2.5"/>
    <rect x="35" y="{base}" width="430" height="17" fill="#9ca3af" rx="2"/>
    <line x1="{px(0)}" y1="{base}" x2="{px(0)}" y2="{py(h)}" stroke="#526579" stroke-dasharray="4 4"/>
    <line x1="{px(0)}" y1="{base-2}" x2="{ax}" y2="{base-2}" stroke="#11845e" stroke-width="3"/>
    {radius_line}
    <line x1="{ax}" y1="{ay}" x2="{tx}" y2="{ty}" stroke="#a85c10" stroke-width="2"/>
    <polyline points="{arc}" fill="none" stroke="#a85c10" stroke-width="2"/>
    <text x="250" y="24" text-anchor="middle" font-family="sans-serif" font-size="13" fill="#24384d">Contact angle = {theta:g}°</text>
    <text x="250" y="297" text-anchor="middle" font-family="sans-serif" font-size="13" fill="#11845e">Footprint radius r = {r:.3g} Å</text>
    <text x="250" y="315" text-anchor="middle" font-family="sans-serif" font-size="13" fill="#ae4646">Curvature radius R = {R:.3g} Å · height = {h:.3g} Å</text>
    </svg>'''



def atom_types_svg():
    """Symmetric illustrative cross-section, not a reconstruction or atom count.

    All contact-layer atoms are gold (interior) or red (perimeter). Upper
    layers are classified geometrically in the schematic only. The support
    is deliberately light gray and never counted in nanoparticle populations.
    """
    palette = dict(zip(LABELS, COLORS))
    palette['Support'] = '#b9bdc4'
    radius = 12
    dx = 24
    dz = dx * math.sqrt(3) / 2
    cx, base = 180, 182
    row_sizes = (12, 13, 12, 13, 12, 11, 10, 7, 4)
    # Doubled horizontal lattice coordinate permits exact neighbor checks.
    lattice = {(2*i-(n-1), j) for j,n in enumerate(row_sizes) for i in range(n)}
    neighbors = ((-2,0),(2,0),(-1,-1),(1,-1),(-1,1),(1,1))
    keys = dict(zip(palette, ('bulk','surface','interface','perimeter','support')))
    defs = ''.join(
        f'<radialGradient id="{keys[name]}" cx="30%" cy="25%" r="80%">'
        f'<stop offset="0%" stop-color="#ffffff" stop-opacity="0.60"/>'
        f'<stop offset="48%" stop-color="{color}"/>'
        f'<stop offset="100%" stop-color="{color}"/></radialGradient>'
        for name,color in palette.items()
    )
    # Use a solid palette-colored base under the translucent highlight.
    def atom(x,y,site):
        return (f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{radius}" '
                f'fill="{palette[site]}"/>'
                f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{radius}" '
                f'fill="url(#{keys[site]})" stroke="#000000" stroke-opacity="0.18" stroke-width="0.6"/>')
    support = ''.join(atom(cx+(i-(n-1)/2)*dx,base+(j+1)*dz,'Support')
                      for j,n in ((1,15),(0,14)) for i in range(n))
    atoms = []
    for j,n in enumerate(row_sizes):
        for i in range(n):
            u = 2*i-(n-1)
            if j == 0:
                site = 'Perimeter' if i in (0,n-1) else 'Interfacial'
            else:
                site = 'Bulk' if all((u+a,j+b) in lattice for a,b in neighbors) else 'Surface'
            atoms.append(atom(cx+u*dx/2,base-j*dz,site))
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 360 240" '
            f'role="img" aria-label="Cross-section of a supported nanoparticle showing bulk, surface, interfacial and perimeter atoms">'
            f'<defs>{defs}</defs>{support}{"".join(atoms)}</svg>')


def color_dot(color):
    return (f'<span aria-hidden="true" style="display:inline-block;width:11px;height:11px;'
            f'border-radius:50%;background:{color};border:1px solid #8888;'
            f'margin-right:7px;vertical-align:middle;"></span>')


def shared_site_legend():
    items = ''.join(f'<span style="white-space:nowrap;">{color_dot(color)}{label}</span>'
                    for label,color in zip(LABELS,COLORS))
    st.markdown('<div style="display:flex;justify-content:center;flex-wrap:wrap;'
                'gap:10px 22px;font-size:0.9rem;margin:4px 0 14px;">'+items+'</div>',
                unsafe_allow_html=True)


def population_table(composition):
    """Compact, theme-compatible table; all interpolated labels are internal."""
    palette = dict(zip(LABELS,COLORS))
    rows = ''.join(
        f'<tr><td>{color_dot(palette[row["Site"]])}{row["Site"]}</td>'
        f'<td style="text-align:right;">{row["Atoms"]:,.0f}</td>'
        f'<td style="text-align:right;">{row["Fraction (%)"]:.2f}</td></tr>'
        for row in composition.to_dict('records')
    )
    st.markdown('<div style="overflow-x:auto;"><table style="width:100%;font-size:0.88rem;'
                'font-variant-numeric:tabular-nums;"><thead><tr><th>Site</th>'
                '<th style="text-align:right;">Atoms</th><th style="text-align:right;">%</th>'
                '</tr></thead><tbody>'+rows+'</tbody></table></div>',unsafe_allow_html=True)


def population_donut(composition):
    """Responsive arc radii prevent clipping; one shared external legend."""
    data = composition.assign(Order=range(len(composition)))
    arcs = alt.Chart(data).mark_arc(
        innerRadius=alt.ExprRef(expr='min(width, height) * 0.25'),
        outerRadius=alt.ExprRef(expr='min(width, height) * 0.45'),
    ).encode(
        theta=alt.Theta('Atoms:Q'),
        color=alt.Color('Site:N',scale=alt.Scale(domain=LABELS,range=COLORS),legend=None),
        order=alt.Order('Order:Q',sort='ascending'),
        tooltip=['Site',alt.Tooltip('Atoms:Q',format=',.0f'),
                 alt.Tooltip('Fraction (%):Q',format='.2f')],
    )
    center = pd.DataFrame({'Total':[metric_value(composition['Atoms'].sum())],
                           'Label':['Total atoms']})
    base = alt.Chart(center)
    total = base.mark_text(fontWeight='bold',fontSize=alt.ExprRef(
        expr='min(21, min(width, height)*0.75/max(1,length(datum.Total)))'),dy=-8).encode(
        text='Total:N')
    label = base.mark_text(fontSize=12,dy=12).encode(text='Label:N')
    return (arcs + total + label).properties(height=260).configure_view(stroke=None)

def signature(df, mode):
    return hashlib.sha256(df.to_csv(index=False).encode()+mode.encode()).hexdigest()


def mode_picker(key):
    return st.radio('Counting mode', ['volume','area'], horizontal=True,key=key,
                    help='Volume mode estimates interfacial/perimeter populations from volumes; area mode uses footprint areas. Both estimate total from volume and exposed surface from area.')


def execute(frame, mode, workflow, source, sig):
    try:
        with st.spinner('Calculating atom populations…'):
            df = normalize_input(frame)
            raw = run_counter(df, mode)
            # Display the user's original radius definitions, not the normalized CLI file.
            original = frame.copy().reset_index(drop=True)
            original.columns = [str(c).strip() for c in original.columns]
            original = original.reindex(columns=HEADERS)
            result = enrich_result(raw,original,mode)
        st.session_state[workflow+'_result'] = dict(data=result, signature=sig, source=source, mode=mode)
        track('calculation_completed',workflow,source,len(result))
    except (ValueError, RuntimeError, OSError, subprocess.TimeoutExpired) as err:
        st.session_state.pop(workflow+'_result',None)
        st.error(str(err))
    except Exception:
        st.session_state.pop(workflow+'_result',None)
        st.error('Calculation could not be completed. Check the input and try again.')
        LOG.warning('Calculation failed with an unexpected error.')


def metric_value(value):
    return f'{value:,.0f}' if np.isfinite(value) else 'Unavailable'


def render_results(workflow, sig):
    record = st.session_state.get(workflow+'_result')
    if not record:
        return
    if record['signature'] != sig:
        st.info('Inputs changed. Run the calculation again to display updated results.')
        return
    out = record['data']
    valid = out['Physical partition']
    if (~valid).any():
        st.warning(f'{int((~valid).sum()):,} row(s) have negative, non-finite or zero-total estimates. Raw counts are retained; these rows are excluded from fraction charts. Check small particles, contact angles and facets.')
    st.subheader('Calculated results')
    st.caption(f"{len(out):,} nanoparticle(s) · {record['mode']} mode · {record['source']} input")
    if workflow == 'single':
        cells = st.columns(5)
        for col, name in zip(cells,['Total','Bulk','Surface','Interface','Perimeter']):
            col.metric('Interfacial' if name=='Interface' else name,metric_value(out.iloc[0][name]))
    else:
        cells = st.columns(3)
        cells[0].metric('Particles analyzed',f'{len(out):,}')
        cells[1].metric('Total atoms (valid partitions)',metric_value(out.loc[valid,'Total'].sum()) if valid.any() else np.nan)
        denom = out.loc[valid,'Total'].sum()
        cells[2].metric('Surface fraction (valid atoms)',f"{100*out.loc[valid,'Surface'].sum()/denom:.2f}%" if denom>0 else 'Unavailable')
    st.caption('Bulk = Total − Surface − Interface − Perimeter (derived, not a native CLI output). Interface and surface exclude perimeter in the package model.')
    if valid.any():
        aggregate = out.loc[valid,SITES].sum()
        composition = pd.DataFrame({'Site':LABELS,'Atoms':aggregate.to_numpy()})
        composition['Fraction (%)'] = 100*composition['Atoms']/composition['Atoms'].sum()
        illustration, left, right = st.columns([1, 1, 1.2])
        with illustration:
            st.markdown('**Atom-type definitions**')
            st.image(atom_types_svg(),use_container_width=True)
            st.caption('Schematic cross-section · light gray = support (not counted).')
        with left:
            st.markdown('**Atom population fractions**')
            st.altair_chart(population_donut(composition),use_container_width=True)
        with right:
            st.markdown('**Population breakdown**')
            population_table(composition)
            if workflow=='batch':
                st.caption('Atom-weighted fractions: sum of each site population divided by sum of total atoms over valid rows.')
        shared_site_legend()
    if workflow=='batch':
        st.markdown('**Population distributions**')
        fields = ['Total','Bulk','Surface','Interface','Perimeter'] + [s+' fraction (%)' for s in SITES]
        field = st.selectbox('Quantity to plot',fields,key='distribution_field')
        values = out.loc[valid,field].to_numpy(dtype=float)
        if len(values):
            hist, edges = np.histogram(values,bins=min(40,max(5,int(np.sqrt(len(values))))))
            bins = pd.DataFrame({'From':edges[:-1],'To':edges[1:],'Particles':hist})
            st.altair_chart(alt.Chart(bins).mark_bar(color='#347bb8').encode(
                x=alt.X('From:Q',bin='binned',title=field),x2='To:Q',y='Particles:Q',tooltip=['From','To','Particles']).properties(height=230),use_container_width=True)
        else:
            st.caption('No valid partitions to plot.')
    st.download_button('Download results (.csv)',out.to_csv(index=False).encode(),
                       file_name=f'atom_count_{workflow}_output.csv',mime='text/csv',
                       key=workflow+'_download',on_click=track,
                       args=('results_downloaded',workflow,record['source'],len(out)))
    st.download_button('Download original four-column CLI output (.csv)',out[COUNTS].to_csv(index=False).encode(),
                       file_name='atom_count_output.csv',mime='text/csv',key=workflow+'_legacy_download',
                       on_click=track,args=('results_downloaded',workflow,record['source'],len(out)))
    with st.expander('View results table',expanded=workflow=='batch'):
        st.dataframe(out.head(1000),hide_index=True,use_container_width=True)
        if len(out)>1000:
            st.caption('First 1,000 rows shown; CSV contains every row.')


def reset_single():
    for k,v in EXAMPLE.items():
        st.session_state['single_'+k] = v
    st.session_state['single_example_requested'] = True
    track('example_selected','single','example')


def facet_select(label,el,key):
    opts=(['Default (unknown)'] if default_facet(el) else [])+sorted(interplanar_dist[el])
    if st.session_state.get(key) not in opts:
        st.session_state[key]=opts[0]
    return st.selectbox(label,opts,key=key)


def single_page():
    st.caption('Enter a geometry directly, or run a built-in platinum example in one click.')
    for k,v in EXAMPLE.items():
        st.session_state.setdefault('single_'+k,v)
    st.button('Try example: Pt, 25 Å, 109°',on_click=reset_single,key='single_example_button')
    left,right=st.columns([1,1])
    with left:
        st.subheader('Nanoparticle parameters')
        el=st.selectbox('Element',ELEMENTS,key='single_element')
        rt=st.radio('Radius definition',['Footprint radius (r)','Curvature radius (R)'],key='single_radius_type')
        a,b=st.columns([2,1])
        radius=a.number_input('Radius',min_value=0.001,format='%.3f',step=1.0,key='single_radius')
        unit=b.selectbox('Unit',['Å','nm'],key='single_unit',help='Changing units changes how the displayed number is interpreted.')
        theta=st.number_input('Contact angle (°)',min_value=0.1,max_value=179.9,step=1.0,format='%.2f',key='single_theta')
        mode=mode_picker('single_mode')
        with st.expander('Advanced settings: optional facets'):
            interface=facet_select('Interface facet',el,'single_interface')
            surface=facet_select('Surface facet',el,'single_surface')
            st.caption('Unknown facets use the package defaults where available; otherwise select a tabulated facet.')
        calculate=st.button('Calculate atom populations',type='primary',key='single_calculate',use_container_width=True)
    radius_A=radius*(10 if unit=='nm' else 1)
    footprint=radius_A if rt.startswith('Footprint') else radius_A*math.sin(math.radians(theta))
    with right:
        st.subheader('Nanoparticle geometry')
        st.image(geometry_svg(radius_A,theta,rt),use_container_width=True)
        st.caption('Spherical-cap cross-section, drawn to scale within the panel. Not an atomistic reconstruction.')
        if footprint<5:
            st.warning('Footprint radius below 5 Å: spherical-cap approximation may be unreliable.')
    frame=pd.DataFrame([[radius_A if rt.startswith('Footprint') else np.nan,
                        radius_A if rt.startswith('Curvature') else np.nan,theta,el,
                        '' if interface.startswith('Default') else interface,
                        '' if surface.startswith('Default') else surface]],columns=HEADERS)
    sig=signature(frame,mode)
    is_example=(el=='Pt' and rt==EXAMPLE['radius_type'] and math.isclose(radius_A,25) and math.isclose(theta,109)
                and interface.startswith('Default') and surface.startswith('Default'))
    if calculate or st.session_state.pop('single_example_requested',False):
        execute(frame,mode,'single','example' if is_example else 'custom',sig)
    render_results('single',sig)


def batch_page():
    st.caption('Analyze a collection of supported nanoparticles from CSV or Excel. Radius values are in Å; Theta is in degrees.')
    a,b=st.columns(2)
    with a:
        st.download_button('Download sample input (.csv)',sample_bytes(),file_name='sample_input.csv',mime='text/csv',
                           on_click=track,args=('sample_downloaded','batch','example'),key='batch_sample_download')
    with b:
        st.download_button('Download blank template (.csv)',(','.join(HEADERS)+'\n').encode(),file_name='input_template.csv',mime='text/csv',
                           on_click=track,args=('template_downloaded','batch','example'),key='batch_template_download')
    source_choice=st.radio('Data source',['Upload my dataset','Built-in example'],horizontal=True,key='batch_source_choice')
    mode=mode_picker('batch_mode')
    uploaded=None
    if source_choice=='Upload my dataset':
        uploaded=st.file_uploader('Choose CSV, XLS or XLSX',type=['csv','xls','xlsx'],key='batch_upload')
        if uploaded is None:
            st.info('Upload a dataset or switch to the built-in example to get started.')
            return
        file_sig=hashlib.sha256(uploaded.getbuffer()).hexdigest()
        suffix=Path(uploaded.name).suffix.lower()
        cache=st.session_state.get('batch_parsed')
        if not cache or cache[0]!=(file_sig,suffix):
            try:
                frame=read_input(uploaded.getvalue(),suffix)
                # Byte-identical distributed samples are not counted as custom uploads.
                demo=(file_sig==hashlib.sha256(sample_bytes()).hexdigest() or file_sig==hashlib.sha256(demo_data().to_csv(index=False).encode()).hexdigest())
                source='example' if demo else 'custom'
                st.session_state['batch_parsed']=((file_sig,suffix),frame,source)
                track('dataset_uploaded','batch',source)
            except Exception as err:
                st.error(f'Unable to read the file: {str(err)[:400]}')
                return
        _,frame,source=st.session_state['batch_parsed']
        sig=file_sig+suffix+mode
    else:
        frame=demo_data();source='example';sig='built_in_'+mode
    st.markdown('**Input preview**')
    st.dataframe(frame.head(10),hide_index=True,use_container_width=True)
    st.caption(f'{len(frame):,} rows · First 10 displayed. The first worksheet is used for Excel files.')
    if st.button('Calculate all nanoparticles' if source=='custom' else 'Run example dataset',type='primary',key='batch_calculate'):
        execute(frame,mode,'batch',source,sig)
    render_results('batch',sig)


def resources():
    st.divider()
    st.subheader('Methods and resources')
    with st.expander('How to obtain geometry from TEM images'):
        st.markdown('Measure the particle–support contact width and divide by two to obtain the footprint radius **r**. Alternatively fit a circle to the particle outline to estimate curvature radius **R**. Measure the contact angle through the particle. The model uses r = R sin θ. Use Å or nm for single input; batch radii must be in Å.')
        a,b=st.columns(2)
        for col,name,caption in [(a,'Acute_1.png','θ < 90°'),(b,'Obtuse_1.png','θ > 90°')]:
            if (DOCS/name).is_file():
                col.image(str(DOCS/name),caption=caption,use_container_width=True)
    with st.expander('How atom populations are estimated'):
        st.markdown('The package approximates the nanoparticle as a spherical cap. The total is estimated from volume and bulk atomic density. The interface and perimeter are estimated from volumes or footprint areas depending on mode; the exposed non-interface surface is estimated from area and facet-dependent density. Bulk is obtained by subtraction. See the preprint for equations and validation.')
        st.image(atom_types_svg(),caption='Schematic cross-section · light gray = support (not counted).',width=360)
        shared_site_legend()
    with st.expander('Batch input format and facet conventions'):
        st.code(','.join(HEADERS),language='text')
        st.markdown('Provide Theta, Element and either radius for every row. If both radii are supplied, a positive r takes precedence. Zero or blank r uses R. Optional facet columns can be omitted or left blank. Use a tabulated three-index facet such as `(1, 1, 1)`; quote fields containing commas in CSV. Do not use extra columns. The web app resolves mixed radius definitions before calling the CLI.')
    with st.expander('Model assumptions, accuracy and limitations'):
        st.markdown('Estimates assume a monometallic spherical-cap geometry and bulk-derived material constants. Faceting, reconstruction, strain, non-spherical shapes and very small particles may invalidate these assumptions. Estimated populations are not direct atomic observations. The model does not identify catalytically active sites or calculate TOFs directly. Consult the preprint for validation and use independent measurements when needed.')
    st.markdown('**Using this tool in your research?**')
    st.code(CITATION,language='text')
    st.markdown(f'[Read the preprint]({PAPER}) · [Source and issues]({REPO}) · [PyPI]({PYPI})')
    st.caption('Developed by Gbolagade Olajide and Tibor Szilvási · The University of Alabama · MIT license')


def main():
    st.set_page_config(page_title='NanoparticleAtomCounter',page_icon='🧮',layout='wide',initial_sidebar_state='collapsed')
    st.title('NanoparticleAtomCounter')
    st.markdown('**From TEM-derived nanoparticle geometry to catalytic site populations.**')
    st.write('Estimate total, bulk, surface, interfacial and perimeter atom populations in supported monometallic nanoparticles. No explicit atomistic model or installation required.')
    st.markdown(f'[Documentation]({REPO}#readme) · [GitHub]({REPO}) · [Preprint]({PAPER}) · [Python package]({PYPI})')
    if is_test_mode():
        st.info('Test mode active — usage statistics will not be recorded.')
    # Preserve parameter widgets when switching workflows (Streamlit otherwise removes
    # state for widgets that were not rendered on the current run).
    for key in ['single_'+k for k in EXAMPLE]+['batch_source_choice','batch_mode','distribution_field']:
        if key in st.session_state:
            st.session_state[key] = st.session_state[key]
    # Radio instead of st.tabs: only the selected workflow executes on reruns.
    workflow=st.radio('Workflow',['Single nanoparticle','Batch analysis'],horizontal=True,label_visibility='collapsed',key='workflow')
    if workflow=='Single nanoparticle':
        single_page()
    else:
        batch_page()
    resources()


if __name__=='__main__':
    main()
