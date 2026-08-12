"""
Biomass Pyrolysis Digital Twin
==============================
Author : Hesam Pero
Stack  : Python · Streamlit · Plotly
Experiment residence time: 15 minutes
"""

import math
import warnings
warnings.filterwarnings('ignore')

import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ──────────────────────────────────────────────────────────────────────────────
#  PAGE CONFIG
# ──────────────────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Pyrolysis Digital Twin",
    page_icon="🌿",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ──────────────────────────────────────────────────────────────────────────────
#  CUSTOM CSS
# ──────────────────────────────────────────────────────────────────────────────
st.markdown("""
<style>
    .stApp { background-color: #ffffff; color: #1a1a1a; }
    [data-testid="stSidebar"] {
        background-color: #f0f7f0;
        border-right: 1px solid #4a7c59;
        color: #1a3a1a;
    }
    [data-testid="stSidebar"] * { color: #1a3a1a !important; }
    [data-testid="stSidebar"] h2 { color: #2d6a4f !important; }
    .banner {
        background: linear-gradient(135deg, #e8f5e9, #c8e6c9);
        padding: 18px 26px; border-radius: 12px;
        border-left: 5px solid #4a7c59; margin-bottom: 20px;
    }
    .banner h1 { color: #1b5e20; font-family: monospace; letter-spacing: 2px;
                 font-size: 1.6rem; margin: 0; }
    .banner p  { color: #2e7d32; font-size: 0.82rem; margin: 5px 0 0;
                 font-family: monospace; }
    .metric-card {
        background: #f1f8f1; border: 1px solid #4a7c59;
        border-radius: 10px; padding: 12px 16px; text-align: center;
    }
    .metric-val { font-size: 1.5rem; font-weight: bold; font-family: monospace; }
    .metric-lbl { font-size: 0.72rem; color: #ffffff; margin-top: 3px; }
    .err-box {
        background: #fff0f0; border: 1px solid #ff4d4d;
        border-radius: 8px; padding: 12px 16px;
        color: #cc0000; font-weight: bold; margin-bottom: 12px;
        font-family: monospace;
    }
    .sec-hdr {
        font-family: monospace; font-size: 0.85rem;
        color: #2d6a4f; letter-spacing: 1px;
        border-bottom: 1px solid #4a7c59; padding-bottom: 4px;
        margin: 18px 0 10px;
    }
    .info-box {
        background: #f1f8f1; border: 1px solid #4a7c59;
        border-radius: 8px; padding: 10px 16px;
        font-family: monospace; font-size: 13px;
        color: #1b5e20; margin-bottom: 12px;
    }
    #MainMenu, footer { visibility: hidden; }
    header { visibility: visible !important; background-color: transparent !important; }
    [data-testid="collapsedControl"] { display: block !important; visibility: visible !important; }
</style>
""", unsafe_allow_html=True)
# ──────────────────────────────────────────────────────────────────────────────
#  CONSTANTS
# ──────────────────────────────────────────────────────────────────────────────
MOLAR_MASS = {'H2': 2.016, 'CH4': 16.043, 'CO2': 44.009,
              'CO': 28.01, 'C2H4': 28.054, 'C2H6': 30.07, 'H2O': 18.015}

LHV_KWH = {'H2': 33.3, 'CH4': 13.9, 'CO2': 0.0,
            'CO': 2.81, 'C2H4': 13.096, 'C2H6': 13.191}

GAS_SPECIES = ['H2', 'CH4', 'CO2', 'CO', 'C2H4', 'C2H6']

GAS_COLORS = ['#e8f4f8', '#ffd166', '#ef476f', '#06d6a0', '#118ab2', '#073b4c']

BIOMASS_COLORS = {
    'Cellulose':     '#4CAF50',
    'Hemicellulose': '#FFC107',
    'Lignin':        '#9C27B0',
    'Ash':           '#9E9E9E',
    'Others':        '#FF5722',
}

CLR = dict(
    biomass='#4a7c59', heart='#2d6a4f', biochar='#1a1a2e',
    tar='#c9722b',     aqueous='#2176ae', gas='#c9a800',
    loss='#9e2a2b',    elec='#6a0572',
    bg='#ffffff',      text='#1a1a1a',   grid='#e0e0e0',
)

LAYOUT_BASE = dict(
    paper_bgcolor='#ffffff',
    plot_bgcolor='#ffffff',
    font=dict(color='#000000', family='monospace', size=13),
)

# ──────────────────────────────────────────────────────────────────────────────
#  cp VALUES — fixed from Excel Constants sheet Row 19
#  Temperature: 535.65 K (mean between room temp and HTT=500°C)
# ──────────────────────────────────────────────────────────────────────────────
cp_H2   = 14.442830645423767
cp_CH4  =  3.0222403938728393
cp_CO   =  1.07178079693991
cp_CO2  =  1.0372651776155373
cp_C2H4 =  2.333051279469571
cp_C2H6 =  2.7237546431758823
cp_H2O  =  1.9766043522946775

# ──────────────────────────────────────────────────────────────────────────────
#  MATHEMATICAL ENGINE  (JMP regressions — validated against Excel)
# ──────────────────────────────────────────────────────────────────────────────
def run_model(cel_pct, hemi_pct, lig_pct, ash_pct, HTT):
    oth_pct = 100 - cel_pct - hemi_pct - lig_pct - ash_pct
    cel  = cel_pct  / 100
    hemi = hemi_pct / 100
    lig  = lig_pct  / 100
    ash  = ash_pct  / 100
    oth  = oth_pct  / 100

    cn  = (cel  - 0.3)   / 0.565
    hn  = (hemi - 0.05)  / 0.565
    ln  = (lig  - 0.05)  / 0.565
    an  = (ash  - 0.005) / 0.565
    on  = (oth  - 0.03)  / 0.565
    T   = (HTT  - 650)   / 200

    # ── MASS BALANCE ──────────────────────────────────────────────────────
    I3 = (100.07696893*cn + 28.950376027*hn + 31.433368789*ln
          + 173.24727769*an - 23.99459131*on
          - 6.410677039*cn*T - 4.338525041*hn*T - 7.781643193*ln*T
          - 128.2443251*an*cn + 1.7040594563*T*T
          - 124.6929765*hn*cn - 112.9481343*ln*cn
          - 76.35203162*cn*cn - 534.8065211*an*an)

    I10 = math.exp(
          -29.9999602*ln - 97.12667684*an + 12.593132802*on
          - 3.906706984*cn*T - 3.443436681*hn*T - 3.895145492*ln*T
          - 5.244721856*an*T
          + 96.974248524*an*cn + 97.126882131*an*hn + 143.23023662*an*ln
          + 0.8306168801*T*T + 2.8712773853*T
          + 25.761508721*ln*cn + 34.2544335543*ln*hn
          + 2.6715461702*cn*cn + 31.907813392*ln*ln + 117.79141504*an*an)

    I14 = (-80.90048365*cn + 15.745437532*hn + 12.963984562*ln
           + 53.160706607*on
           - 13.66319331*cn*T - 11.99377561*hn*T - 11.67588784*ln*T
           + 130.8847331*an*cn
           - 5.59357011*T*T
           + 122.1574317*hn*cn + 120.0343408*ln*cn
           + 97.881342875*cn*cn)

    I18 = 100 - I14 - I10 - I3

    # ── BIOCHAR ───────────────────────────────────────────────────────────
    I4 = (-162.5057886*cn - 294.0093995*hn - 35.37091753*an
          + 252.32393076*on + 7.3133995143*ln*T
          + 566.80541264*an*hn + 688.4628714*hn*cn
          + 361.45673268*ln*cn + 466.45177908*ln*hn
          + 247.81796519*cn*cn + 349.06058638*hn*hn + 81.9783256*ln*ln)

    I5 = (-15.84947478*an + 4.2020698881*on
          - 1.542833187*cn*T - 1.635174427*hn*T - 1.812258897*ln*T
          + 14.607931576*an*cn + 31.143748976*an*hn + 22.723814867*an*ln
          + 0.2599514383*T*T + 0.4272613591*T*T*T
          + 3.7469305279*hn*cn + 3.0889333841*ln*cn + 2.3207868942*ln*hn
          + 2.4088741972*cn*cn + 1.7206687364*ln*ln)

    I6 = (45.996502831*cn + 42.612413868*ln + 301.73430335*an
          + 16.636098788*on
          - 8.547139184*cn*T - 7.50272576*hn*T - 7.740743509*ln*T
          - 239.4503722*an*cn + 276.63070088*an*hn
          + 1.6713741407*T*T - 31.51342926*ln*cn
          - 1345.109709*an*an)

    I7 = (90.635198223*on
          - 626.3108033*an*hn - 544.9369016*an*ln
          + 238.29830027*hn*cn + 140.10772128*ln*cn + 194.72664989*ln*hn
          - 9.921383266*T
          + 113.68808255*ln*ln + 1354.6559674*an*an)

    I8 = (85.108788292*hn + 6.6022914528*ln + 382.68366539*an
          - 5.167216364*on + 14.999921222*an*T
          - 743.2042477*an*hn - 312.046245*an*ln
          - 106.8141973*hn*cn - 71.61900629*ln*hn)

    I9 = (0.4507287165*on
          - 0.258762865*cn*T - 0.187412695*hn*T
          - 0.261243935*ln*T - 0.341510863*an*T
          + 1.3387205768*an*cn + 0.6674713435*an*ln
          + 0.0477768542*T*T + 0.0629265598*T*T*T
          + 0.2417285248*ln*cn + 0.5133932375*ln*hn
          + 0.416321673*cn*cn + 0.7372986154*hn*hn + 0.2784252064*ln*ln)

    # ── TAR ───────────────────────────────────────────────────────────────
    I11 = (-470.7491011*hn + 234.80940495*on
           - 10.89833872*ln*T + 1030.1633129*an*hn
           - 48.66917338*T*T - 37.31435606*T*T*T
           + 488.7759288*hn*cn + 141.50057855*ln*cn + 612.65929356*ln*hn
           + 24.102481225*T
           + 93.118307643*cn*cn + 615.34505423*hn*hn + 85.868931308*ln*ln)

    I12 = math.exp(
          -2.948160732*ln - 47.55296133*an + 7.6205515382*on
          - 6.786652577*cn*T - 5.666710868*hn*T - 6.673456567*ln*T
          - 6.932339492*an*T
          + 51.422379606*an*cn + 78.039572841*an*hn + 69.779775891*an*ln
          - 0.576884247*T*T + 4.6251718751*T
          + 3.2553621984*cn*cn + 3.2519857416*ln*ln)

    # ── AQUEOUS ───────────────────────────────────────────────────────────
    I15 = (144.07763825*hn + 59.227141304*ln - 324.0226646*an
           + 33.837213312*on
           + 75.614077094*cn*T + 90.120888858*hn*T + 82.65854275*ln*T
           + 346.19037974*an*cn - 567.2045341*an*hn
           - 23.55709667*T*T - 29.34208733*T*T*T
           + 92.572476195*ln*cn - 120.7692073*ln*hn
           - 39.43890503*T + 2077.8484171*an*an)

    I16 = math.exp(
          -2.750576975*hn - 13.61759629*an + 4.0643985699*on
          - 2.122374314*cn*T - 0.918186366*hn*T - 1.468039658*ln*T
          + 42.884855363*an*hn + 16.790284277*an*ln
          + 5.2715671295*hn*cn + 3.0738167897*cn*cn)

    I17 = math.exp(
          -5.486667622*an + 1.1862299906*on
          - 2.95852876*cn*T - 1.53044809*hn*T - 2.119624719*ln*T
          + 16.030974119*an*hn - 0.597514131*T*T
          + 1.6225791635*cn*cn)

    I19 = 100 - I17 - I12 - I6

    # ── GAS SPECIES ───────────────────────────────────────────────────────
    I33 = (50.855159479*cn + 21.646502907*ln + 275.9693742*an
           + 4.3027055929*on
           + 19.961414236*cn*T + 10.623075841*hn*T + 19.891669938*ln*T
           + 34.500211966*an*T
           - 317.3419658*an*cn - 179.6776216*an*ln - 30.45879796*ln*hn
           - 42.55445372*cn*cn - 581.5419592*an*an)

    I34 = (-124.6443434*cn - 197.8981344*hn - 188.1134822*ln
           - 249.8742062*an + 155.04051965*on
           + 13.311471091*an*T
           + 349.33334106*an*cn + 701.52748458*an*hn + 532.81532125*an*ln
           - 3.172169954*T*T
           + 311.05930624*hn*cn + 310.89568005*ln*cn + 398.20533775*ln*hn
           + 155.08078189*cn*cn + 206.93726426*hn*hn + 201.01125112*ln*ln)

    I35 = (77.687731093*cn + 85.803262875*hn + 40.169831958*ln
           + 644.10469991*an - 44.22969501*on
           + 47.687946473*cn*T + 53.313145469*hn*T + 43.404151825*ln*T
           - 597.9665234*an*cn - 1307.605344*an*hn - 775.1536572*an*ln
           - 11.44612943*T*T*T - 45.1290387*T
           - 93.63825504*cn*cn)

    I36 = (-237.9225838*hn - 282.5214296*an + 120.0412427*on
           - 47.10820214*cn*T - 56.07720167*hn*T - 48.13515755*ln*T
           + 280.56731448*an*cn + 842.35747754*an*hn + 391.66019226*an*ln
           + 2.4572238598*T*T + 11.069453743*T*T*T
           + 171.78587293*hn*cn + 304.54636978*ln*hn
           + 29.028806444*T
           + 88.23001841*cn*cn + 339.71714902*hn*hn + 28.513836366*ln*ln)

    I37 = math.exp(
          -21.25453061*cn - 10.59734423*hn - 15.06147195*an
          + 13.798540129*on
          - 0.924312937*cn*T - 0.619626563*ln*T
          + 46.27633626*an*cn - 0.815055262*T*T
          + 32.939700947*hn*cn + 22.576791381*ln*cn + 13.828381556*ln*hn
          + 0.5795252254*T
          + 21.908118863*cn*cn + 17.379936999*hn*hn + 0.5434783467*ln*ln
          + 53.811983324*an*an)

    I38 = math.exp(
          -2.949302411*hn - 3.093539556*ln + 1.0221979581*on
          - 0.725272009*cn*T
          + 39.403582609*an*hn + 20.404976472*an*ln
          - 0.777644968*T*T - 0.792633552*T*T*T
          + 1.8422012948*cn*cn + 1.8780691224*ln*ln
          - 68.654788896*an*an)

    # I33–I38 are direct JMP vol% outputs — use as-is, no normalisation
    vol_pct = [I33, I34, I35, I36, I37, I38]

    # wt% from vol% × molar mass, then normalise
    MM = [MOLAR_MASS[s] for s in GAS_SPECIES]
    wt_raw  = [v * m for v, m in zip(vol_pct, MM)]
    wt_sum  = sum(wt_raw)
    wt_pct  = [w / wt_sum * 100 for w in wt_raw]
    gas_yields = [I18 * w / 100 for w in wt_pct]

    # ── ENERGY BALANCE ────────────────────────────────────────────────────
    O4 = (20569.574376*ln - 40588.17312*an + 30168.980335*on
          + 56600.930392*an*cn + 37822.238406*hn*cn
          + 17652.766485*ln*cn + 24695.159497*ln*hn
          + 13196.638859*cn*cn + 24915.661943*hn*hn + 202199.21498*an*an)

    O5 = (-91724.17531*an + 56812.209158*on
          + 3215.5189824*cn*T + 2709.2335479*ln*T
          + 245300.29378*an*hn + 126788.53215*an*ln
          - 1383.349521*T*T
          + 83696.657992*hn*cn + 68409.814186*ln*cn + 56124.660921*ln*hn
          - 1732.256896*T
          + 34249.084397*cn*cn + 31300.891765*ln*ln)

    O6  = O5 * I3 * 0.01 * 1000 * 0.001 * 0.001
    O7  = (-166507.3282*hn + 98821.570032*on
           + 351831.42247*an*hn
           - 17733.40785*T*T - 16525.00372*T*T*T
           + 135888.97323*hn*cn + 57398.20583*ln*cn + 217379.71751*ln*hn
           + 9076.653401*T
           + 52017.182345*cn*cn + 267227.42437*hn*hn + 36088.8843*ln*ln)

    O8   = O7 * 0.001 * 0.001 * 1000 * I10 * 0.01
    O10  = I14 * 0.01 * 2257 * 1000 * 0.001 * 0.001
    O11  = I14 * 0.01 * 1000 * cp_H2O * (HTT - 25) * 0.001 * 0.001
    O12  = sum(wt_pct[i] * LHV_KWH[GAS_SPECIES[i]] for i in range(6)) / 100 * 3600
    O13  = O12 * 0.001 * 0.001 * 1000 * I18 * 0.01
    O14  = (wt_pct[0]*0.01*cp_H2
          + wt_pct[1]*0.01*cp_CH4
          + wt_pct[2]*0.01*cp_CO2
          + wt_pct[3]*0.01*cp_CO
          + wt_pct[4]*0.01*cp_C2H4
          + wt_pct[5]*0.01*cp_C2H6)
    O15  = I18 * 0.01 * 1000 * O14 * (HTT - 25) * 0.001 * 0.001

    O16  = (-24.8335076*an + 2.4880741103*on
            - 4.067270118*hn*T
            + 98.753352012*an*hn + 54.458690797*an*ln
            + 2.2074234994*T*T + 6.9986304371*ln*cn
            + 4.8447561648*T
            + 13.265134491*cn*cn + 4.6531497438*ln*ln)

    O17  = (O4 * 1000 * 0.001 * 0.001 + O16) - O6 - O8 - O13 - O15 - O10 - O11

    return {
        'cel': cel_pct, 'hemi': hemi_pct, 'lig': lig_pct,
        'ash': ash_pct, 'oth': oth_pct, 'HTT': HTT,
        'biochar_yield': I3, 'tar_yield': I10,
        'aqueous_yield': I14, 'gas_yield': I18,
        'bc_C': I4, 'bc_H': I5, 'bc_yC': I6, 'bc_yN': I7,
        'bc_ash': I8, 'bc_HC': I9,
        'tar_C': I11, 'tar_yC': I12,
        'aq_yN': I15, 'aq_C': I16, 'aq_yC': I17,
        'gas_yC': I19,
        'gas_vol_pct': dict(zip(GAS_SPECIES, vol_pct)),
        'gas_wt_pct':  dict(zip(GAS_SPECIES, wt_pct)),
        'gas_yields':  dict(zip(GAS_SPECIES, gas_yields)),
        'bm_LHV': O4, 'bc_LHV': O5,
        'bm_energy': O4 * 1000 * 0.001 * 0.001,
        'bc_energy': O6, 'tar_HHV': O7, 'tar_energy': O8,
        'aq_cond': O10, 'aq_sens': O11,
        'gas_LHV': O12, 'gas_chem': O13, 'gas_sens': O15,
        'elec_input': O16, 'losses': O17,
    }

# ──────────────────────────────────────────────────────────────────────────────
#  CHART BUILDERS
# ──────────────────────────────────────────────────────────────────────────────
def chart_biomass_pie(r):
    labels = list(BIOMASS_COLORS.keys())
    values = [r['cel'], r['hemi'], r['lig'], r['ash'], r['oth']]
    colors = list(BIOMASS_COLORS.values())
    fig = go.Figure(go.Pie(
        labels=labels, values=values,
        marker=dict(colors=colors, line=dict(color='#0f1117', width=2)),
        hole=0.42,
        textinfo='label+percent',
        textposition='outside',
        textfont=dict(size=12),
        hovertemplate='<b>%{label}</b><br>Share: %{value:.1f}%<extra></extra>',
        pull=[0.03]*5, rotation=90,
    ))
    fig.update_layout(
        **LAYOUT_BASE,
        title=dict(text='🌾  Biomass Composition',
                   font=dict(size=14, color='#a8d8a8'), x=0.5, xanchor='center'),
        height=340, margin=dict(l=10, r=10, t=50, b=10),
        showlegend=True,
        legend=dict(bgcolor='rgba(0,0,0,0)', font=dict(size=11),
                    orientation='v', x=1.0, y=0.5),
    )
    return fig
def chart_mass_sankey(r):
    labels = [
        'Biomass\n1.00 kg', 'Reactor',
        f'Biochar\n{r["biochar_yield"]/100:.2f} kg ({r["biochar_yield"]:.2f}%)',
        f'Tar\n{r["tar_yield"]/100:.2f} kg ({r["tar_yield"]:.2f}%)',
        f'Aqueous\n{r["aqueous_yield"]/100:.2f} kg ({r["aqueous_yield"]:.2f}%)',
        f'Gas\n{r["gas_yield"]/100:.2f} kg ({r["gas_yield"]:.2f}%)',
    ]
    hover = [
        (f'Cellulose: {r["cel"]:.1f}%<br>Hemicellulose: {r["hemi"]:.1f}%<br>'
         f'Lignin: {r["lig"]:.1f}%<br>Ash: {r["ash"]:.1f}%<br>Others: {r["oth"]:.1f}%'),
        f'HTT: {r["HTT"]} °C',
        (f'C content: {r["bc_C"]:.2f}%<br>H content: {r["bc_H"]:.2f}%<br>'
         f'C yield: {r["bc_yC"]:.2f}%<br>N yield: {r["bc_yN"]:.2f}%<br>'
         f'Ash: {r["bc_ash"]:.2f}%<br>H/C molar: {r["bc_HC"]:.3f}'),
        f'C content: {r["tar_C"]:.2f}%<br>C yield: {r["tar_yC"]:.2f}%',
        f'N yield: {r["aq_yN"]:.2f}%<br>C yield: {r["aq_yC"]:.2f}%',
        (f'C yield: {r["gas_yC"]:.2f}%<br>'
         + '<br>'.join(f'{s}: {r["gas_wt_pct"][s]:.2f} wt% | {r["gas_vol_pct"][s]:.2f} vol%'
                       for s in GAS_SPECIES)),
    ]
    fig = go.Figure(go.Sankey(
        arrangement='snap',
        node=dict(
            pad=8, thickness=28,
            label=labels,
            color=[CLR['biomass'], '#c9722b', CLR['biochar'],
                   CLR['tar'], CLR['aqueous'], CLR['gas']],
            customdata=hover,
            hovertemplate='%{customdata}<extra>%{label}</extra>',
            line=dict(color='#000000', width=0.5),
        ),
        link=dict(
            source=[0, 1, 1, 1, 1],
            target=[1, 2, 3, 4, 5],
            value=[100, r['biochar_yield'], r['tar_yield'],
                   r['aqueous_yield'], r['gas_yield']],
            color=['rgba(45,106,79,0.45)', 'rgba(26,26,46,0.7)',
                   'rgba(201,114,43,0.7)', 'rgba(33,118,174,0.7)',
                   'rgba(232,197,71,0.7)'],
            hovertemplate='Flow: %{value:.2f}%<extra></extra>',
        ),
    ))
    fig.update_layout(
        paper_bgcolor='#ffffff',
        plot_bgcolor='#ffffff',
        font=dict(color='#000000', size=14, family='Arial Black'),
        title=dict(text='⚖️  Mass Balance — per 1 kg Biomass',
                   font=dict(size=14, color='#1b5e20'), x=0.5, xanchor='center'),
        height=340, margin=dict(l=10, r=10, t=50, b=10),
    )
    return fig


def chart_yield_bar(r):
    names  = ['Biochar', 'Tar', 'Aqueous', 'Gas']
    values = [r['biochar_yield'], r['tar_yield'], r['aqueous_yield'], r['gas_yield']]
    colors = [CLR['biochar'], CLR['tar'], CLR['aqueous'], CLR['gas']]
    fig = go.Figure(go.Bar(
        x=values, y=names, orientation='h',
        marker_color=colors,
        text=[f'{v:.2f}%' for v in values],
        textposition='outside',
        textfont=dict(color=CLR['text']),
        hovertemplate='%{y}: %{x:.3f}%<extra></extra>',
    ))
    fig.update_layout(
        **LAYOUT_BASE,
        xaxis=dict(title='Yield [wt%]', gridcolor=CLR['grid'], range=[0, 108]),
        yaxis=dict(gridcolor=CLR['grid']),
        title=dict(text='📊  Mass Balance',
                   font=dict(size=14, color=CLR['text']), x=0.5, xanchor='center'),
        height=250, margin=dict(l=20, r=20, t=50, b=20),
        showlegend=False,
    )
    return fig


def chart_gas_pies(r):
    fig = make_subplots(
        rows=1, cols=2,
        specs=[[{'type': 'pie'}, {'type': 'pie'}]],
        subplot_titles=['Gas Composition [wt%]', 'Gas Composition [vol%]'],
    )
    common = dict(
        marker=dict(colors=GAS_COLORS, line=dict(color='#0f1117', width=1)),
        hole=0.4, textinfo='label+percent', textposition='outside',
    )
    fig.add_trace(go.Pie(
        labels=GAS_SPECIES,
        values=[r['gas_wt_pct'][s] for s in GAS_SPECIES],
        hovertemplate='<b>%{label}</b><br>%{value:.2f} wt%<extra></extra>',
        **common,
    ), row=1, col=1)
    fig.add_trace(go.Pie(
        labels=GAS_SPECIES,
        values=[r['gas_vol_pct'][s] for s in GAS_SPECIES],
        hovertemplate='<b>%{label}</b><br>%{value:.2f} vol%<extra></extra>',
        **common,
    ), row=1, col=2)
    fig.update_layout(
        **LAYOUT_BASE,
        title=dict(text='🔬  Gas Phase Composition',
                   font=dict(size=14, color='#c9722b'), x=0.5, xanchor='center'),
        height=500, margin=dict(l=60, r=60, t=60, b=100),
        legend=dict(bgcolor='rgba(0,0,0,0)'),
    )
    fig.update_annotations(font_color=CLR['text'])
    return fig


def chart_energy_sankey(r):
    bm_e  = r['bm_energy']
    elec  = r['elec_input']
    bc_e  = r['bc_energy']
    tar_e = r['tar_energy']
    g_tot = r['gas_chem'] + r['gas_sens']
    a_tot = r['aq_cond']  + r['aq_sens']
    loss  = max(r['losses'], 0.001)

    labels = [
        f'Biomass LHV\n{bm_e:.3f} MJ',
        f'Electrical\nInput\n{elec:.3f} MJ',
        'Total\nEnergy\nInput',
        f'Biochar\n{bc_e:.3f} MJ',
        f'Tar\n{tar_e:.3f} MJ',
        f'Gas\n{g_tot:.3f} MJ',
        f'Aqueous\n{a_tot:.3f} MJ',
        f'Losses\n{loss:.3f} MJ',
    ]
    hover = [
        f'Biomass LHV: {r["bm_LHV"]:.1f} J/g<br>Chemical energy: {bm_e:.4f} MJ/kg',
        f'Electrical heating correction<br>{elec:.4f} MJ/kg',
        f'Total input: {bm_e+elec:.4f} MJ/kg',
        f'Biochar LHV: {r["bc_LHV"]:.1f} J/g<br>Energy stored: {bc_e:.4f} MJ/kg',
        f'Tar HHV: {r["tar_HHV"]:.1f} J/g<br>Energy in tar: {tar_e:.4f} MJ/kg',
        f'Chemical: {r["gas_chem"]:.4f} MJ/kg<br>Sensible: {r["gas_sens"]:.4f} MJ/kg',
        f'Condensation: {r["aq_cond"]:.4f} MJ/kg<br>Sensible: {r["aq_sens"]:.4f} MJ/kg',
        f'Non-recovered losses: {loss:.4f} MJ/kg',
    ]
    fig = go.Figure(go.Sankey(
        arrangement='snap',
        node=dict(
            pad=8, thickness=26,
            label=labels,
            color=[CLR['biomass'], CLR['elec'], CLR['heart'],
                   CLR['biochar'], CLR['tar'], CLR['gas'],
                   CLR['aqueous'], CLR['loss']],
            customdata=hover,
            hovertemplate='%{customdata}<extra>%{label}</extra>',
            line=dict(color='#000000', width=0.5),
        ),
        link=dict(
            source=[0, 1, 2, 2, 2, 2, 2],
            target=[2, 2, 3, 4, 5, 6, 7],
            value=[bm_e, elec, bc_e, tar_e, g_tot, a_tot, loss],
            color=['rgba(74,124,89,0.5)', 'rgba(106,5,114,0.5)',
                   'rgba(26,26,46,0.7)', 'rgba(201,114,43,0.7)',
                   'rgba(232,197,71,0.6)', 'rgba(33,118,174,0.6)',
                   'rgba(158,42,43,0.7)'],
            hovertemplate='Flow: %{value:.4f} MJ/kg<extra></extra>',
        ),
    ))
    fig.update_layout(
        paper_bgcolor='#ffffff',
        plot_bgcolor='#ffffff',
        font=dict(color='#000000', size=14, family='Arial Black'),
        title=dict(text='⚡  Energy Balance — MJ per kg Biomass',
                   font=dict(size=14, color='#c9722b'), x=0.5, xanchor='center'),
        height=420, margin=dict(l=20, r=20, t=60, b=20),
    )
    return fig

# ──────────────────────────────────────────────────────────────────────────────
#  SIDEBAR
# ──────────────────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("## ⚙️ Input Parameters")
    st.markdown("---")
    st.markdown("**Biomass Composition [wt%]**")
    cel  = st.slider("🌿 Cellulose",     20.0, 60.0, 54.2, 0.1)
    hemi = st.slider("🍂 Hemicellulose",  5.0, 35.0,  9.4, 0.1)
    lig  = st.slider("🪵 Lignin",         5.0, 60.0, 23.2, 0.1)
    ash  = st.slider("⚪ Ash",            0.0, 10.0,  0.4, 0.1)
    total = cel + hemi + lig + ash
    oth   = 100.0 - total

    if oth < 0:
        st.error(f"❌ Exceeds 100% by **{-oth:.1f}%**")
    elif oth > 17:
        st.error(f"❌ \"Others\" = {oth:.1f}% exceeds the model's valid maximum of 17%")
    else:
        st.success(f"✅ Others: **{oth:.1f}%**")
        st.caption(
            'The fraction "Others" describes components in the biomass '
            'other than Cellulose, Hemicellulose, Lignin, and Ash (e.g., proteins).'
        )

    st.markdown("---")
    st.markdown("**Process Temperature**")
    HTT = st.slider("🌡️ HTT [°C]", 450, 850, 500, 10)
    st.markdown("---")
    st.markdown(
        "<div style='font-size:0.75rem;font-family:monospace;'>"
        "Valid range: HTT 450–850 °C<br>"
        "Basis: 1 kg dry biomass<br>"
        "Residence time: 15 minutes"
        "</div>",
        unsafe_allow_html=True,
    )

# ──────────────────────────────────────────────────────────────────────────────
#  MAIN PANEL
# ──────────────────────────────────────────────────────────────────────────────
st.markdown("""
<div class="banner">
  <h1>🌿 BIOMASS PYROLYSIS DIGITAL TWIN</h1>
  <p>Mass &amp; Energy Balance Model · Co-current of solids and gas phase during pyrolysis · Residence time 15 minutes · Highest Treatment Temperature (HTT) 450–850 °C</p>
</div>
""", unsafe_allow_html=True)

# Guard
if oth < 0 or oth > 17:
    if oth < 0:
        msg = f"Component sum exceeds 100% by {-oth:.1f}%. Please reduce Cellulose, Hemicellulose, Lignin, or Ash."
    else:
        msg = f'"Others" = {oth:.1f}% exceeds the model\'s valid maximum of 17%. Please adjust the inputs so Others ≤ 17%.'
    st.markdown(
        f'<div class="err-box">⚠️ {msg}</div>',
        unsafe_allow_html=True,
    )
    st.stop()

r = run_model(cel, hemi, lig, ash, HTT)

# ── Quick summary ─────────────────────────────────────────────────────────────
st.markdown(
    f'<div class="info-box">📋 <b>Quick Summary</b> &nbsp;|&nbsp; '
    f'Biochar: <b>{r["biochar_yield"]:.2f}%</b> &nbsp;'
    f'Tar: <b>{r["tar_yield"]:.2f}%</b> &nbsp;'
    f'Aqueous: <b>{r["aqueous_yield"]:.2f}%</b> &nbsp;'
    f'Gas: <b>{r["gas_yield"]:.2f}%</b> &nbsp;|&nbsp; '
    f'Biochar LHV: <b>{r["bc_LHV"]:.0f} J/g</b> &nbsp;'
    f'Biomass LHV: <b>{r["bm_LHV"]:.0f} J/g</b></div>',
    unsafe_allow_html=True,
)

# ── Metric cards ──────────────────────────────────────────────────────────────
st.markdown('<div class="sec-hdr">📋 PRODUCT YIELD SUMMARY</div>', unsafe_allow_html=True)
c1, c2, c3, c4 = st.columns(4)
cards = [
    (c1, f"{r['biochar_yield']:.2f}%", "Biochar Yield",  "#3a3a5e", "#a8a8e8"),
    (c2, f"{r['tar_yield']:.2f}%",     "Tar Yield",      "#3a1a0a", "#e8a87c"),
    (c3, f"{r['aqueous_yield']:.2f}%", "Aqueous Yield",  "#0a1a3a", "#7cb8e8"),
    (c4, f"{r['gas_yield']:.2f}%",     "Gas Yield",      "#3a3a0a", "#e8d87c"),
]
for col, val, lbl, bg, fg in cards:
    col.markdown(
        f'<div class="metric-card" style="background:{bg};">'
        f'<div class="metric-val" style="color:{fg};">{val}</div>'
        f'<div class="metric-lbl">{lbl}</div></div>',
        unsafe_allow_html=True,
    )

st.markdown("<br>", unsafe_allow_html=True)

# ── Yield bar ─────────────────────────────────────────────────────────────────
st.plotly_chart(chart_yield_bar(r), use_container_width=True, key="bar")

# ── Biomass pie + Mass Sankey ─────────────────────────────────────────────────
st.markdown('<div class="sec-hdr">⚖️ MASS BALANCE</div>', unsafe_allow_html=True)
col_pie, col_sankey = st.columns([1, 1.6])
with col_pie:
    st.plotly_chart(chart_biomass_pie(r), use_container_width=True, key="bm_pie")
with col_sankey:
    st.plotly_chart(chart_mass_sankey(r), use_container_width=True, key="sankey")

# ── Gas pies ──────────────────────────────────────────────────────────────────
st.markdown('<div class="sec-hdr">🔬 GAS PHASE COMPOSITION</div>', unsafe_allow_html=True)
st.plotly_chart(chart_gas_pies(r), use_container_width=True, key="gas")

# ── Energy Sankey ─────────────────────────────────────────────────────────────
st.markdown('<div class="sec-hdr">⚡ ENERGY BALANCE</div>', unsafe_allow_html=True)
st.plotly_chart(chart_energy_sankey(r), use_container_width=True, key="energy")

# ── Energy summary table ──────────────────────────────────────────────────────
st.markdown('<div class="sec-hdr">📋 ENERGY SUMMARY TABLE</div>', unsafe_allow_html=True)
e1, e2 = st.columns(2)
with e1:
    st.markdown(f"""
    | Input | Value |
    |-------|-------|
    | Biomass LHV | `{r['bm_LHV']:.2f} J/g` |
    | Biomass Chemical Energy | `{r['bm_energy']:.4f} MJ/kg` |
    | Electrical Heating Input | `{r['elec_input']:.4f} MJ/kg` |
    | **Total Input** | **`{r['bm_energy']+r['elec_input']:.4f} MJ/kg`** |
    """)
with e2:
    st.markdown(f"""
    | Output | Value |
    |--------|-------|
    | Biochar Chemical Energy | `{r['bc_energy']:.4f} MJ/kg` |
    | Tar Chemical Energy | `{r['tar_energy']:.4f} MJ/kg` |
    | Gas (Chemical Energy + Sensible Heat) | `{r['gas_chem']+r['gas_sens']:.4f} MJ/kg` |
    | Aqueous (Condensation + Sensible Heat) | `{r['aq_cond']+r['aq_sens']:.4f} MJ/kg` |
    | Non-recovered Losses | `{r['losses']:.4f} MJ/kg` |
    """)

# ── Footer ────────────────────────────────────────────────────────────────────
st.markdown("---")
st.markdown(
    "<div style='text-align:center;font-family:monospace;font-size:0.75rem;"
    "color:#2d6a4f;padding:10px;'>"
    "🌿 Biomass Pyrolysis Digital Twin · Built by <b>Hesam Pero</b> · Python & Streamlit<br>"
    "This digital twin is based on an empirical pyrolysis model. The pyrolysis experiments have been performed within the project &quot;PyFlex - Pyrolysetechnologie als intersektorale Flexibilit&auml;tsma&szlig;nahme im Energiesystem Deutschland&quot; funded by the German Federal Ministry for Economic Affairs and Energy (BMWE) · Residence time: 15 minutes · HTT: 450–850 °C · 1 kg Biomass"
    "</div>",
    unsafe_allow_html=True,
)

