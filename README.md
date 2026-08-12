# 🌿 Biomass Pyrolysis Digital Twin

An interactive web-based dashboard modelling the complete **mass and energy balance** of a slow pyrolysis process for any biomass feedstock, on a **1 kg dry basis**.

Developed at the **Institute for Sustainable Energy Systems (INES), Hochschule Offenburg** under the supervision of **Prof. Heide Biollaz**, as part of the PyFlex project funded by the German Federal Ministry for Economic Affairs and Energy (BMWE).

---

## What it does

Given a biomass biochemical composition and process temperature, the model computes and visualises in real time:

- **Mass balance** — yields of biochar, tar, aqueous phase, and gas
- **Biochar sub-parameters** — carbon/hydrogen content, H/C molar ratio, ash content, nitrogen yield
- **Tar sub-parameters** — carbon content and carbon yield
- **Aqueous sub-parameters** — nitrogen yield, carbon content, carbon yield
- **Gas speciation** — H₂, CH₄, CO₂, CO, C₂H₄, C₂H₆ in both wt% and vol%
- **Energy balance** — chemical energy (LHV/HHV), sensible heat, condensation heat, electrical heating input, and non-recovered losses in MJ/kg

---

## Experimental Conditions

| Parameter | Value |
|---|---|
| Process type | Slow pyrolysis |
| Residence time | 15 minutes |
| HTT range | 450 – 850 °C |
| Basis | 1 kg dry biomass |
| Reactor system | PYREKA / PyFlex |

---

## Input Parameters

| Parameter | Range |
|---|---|
| Cellulose content | 20 – 60 % |
| Hemicellulose content | 5 – 35 % |
| Lignin content | 5 – 60 % |
| Ash content | 0 – 10 % |
| Others | Auto-calculated (100 − sum), valid up to 17% |
| Highest Treatment Temperature (HTT) | 450 – 850 °C |

---

## Mathematical Model

The model uses **multivariate polynomial and exponential regression equations** derived from experimental PYREKA/PyFlex slow pyrolysis data processed using JMP statistical software. All inputs are normalised before regression using the JMP model scaling parameters. Heat capacity values are taken directly from the experimental reference constants.

**Validation against experimental reference (Cel=54.2%, Hemi=9.4%, Lig=23.2%, Ash=0.4%, HTT=500°C):**

| Output | Model | Reference | Status |
|---|---|---|---|
| Biochar yield | 22.2729% | 22.2729% | ✅ |
| Tar yield | 4.2329% | 4.2329% | ✅ |
| Aqueous yield | 23.2297% | 23.2297% | ✅ |
| Gas yield | 50.2644% | 50.2644% | ✅ |
| Biomass LHV | 18777.15 J/g | 18777.15 J/g | ✅ |
| Non-recovered losses (O17) | 7.5113 MJ/kg | 7.5113 MJ/kg | ✅ |

---

## Tech Stack

| Layer | Tools |
|---|---|
| Web framework | Streamlit |
| Visualisation | Plotly (Sankey, Pie, Bar) |
| Physics engine | Python · math · NumPy |
| Regression source | JMP Statistical Software |

---

## Run Locally

```bash
pip install -r requirements.txt
streamlit run app_real.py
```

---

## Deploy on University Server

```bash
streamlit run app_real.py --server.port 8501 --server.address 0.0.0.0
```

Access via browser:
