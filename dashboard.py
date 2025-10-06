import streamlit as st
import lightkurve as lk
from lightkurve import TessTargetPixelFile
import matplotlib.pyplot as plt
from urllib.request import urlopen
from io import BytesIO
from zipfile import ZipFile
import os
import ssl
import numpy as np
import tempfile

ssl._create_default_https_context = ssl._create_unverified_context

# Constants
AU = 149.6e6
SolarConstant = 1361.0
sigma = 5.67e-8
CO2Concentration = 0.0004
CH4Concentration = 1.8e-6

# Page configuration
st.set_page_config(
    page_title="Exoplanet Light Curve Analyzer",
    page_icon="🌍",
    layout="wide"
)

st.title("🔭 Exoplanet Light Curve Analyzer")
st.markdown("*Tools for exoplanet research and temperature estimation of TESS & KEPLER light curves*")

# Initialize session state
if 'lc' not in st.session_state:
    st.session_state.lc = None
if 'obj_name' not in st.session_state:
    st.session_state.obj_name = None
if 'stellar_radius' not in st.session_state:
    st.session_state.stellar_radius = None
if 'star_temp' not in st.session_state:
    st.session_state.star_temp = None
if 'planet_radius' not in st.session_state:
    st.session_state.planet_radius = None
if 'planet_distance' not in st.session_state:
    st.session_state.planet_distance = None
if 'period' not in st.session_state:
    st.session_state.period = None
if 'flux' not in st.session_state:
    st.session_state.flux = None
if 'irradiance' not in st.session_state:
    st.session_state.irradiance = None

# Helper functions
def rolling_avg(arr):
    ret = np.cumsum(arr, dtype=float)
    ret[3:] = ret[3:] - ret[:-3]
    return ret[2:] / 3

def within(arr, roller):
    m = np.mean(roller)
    s = np.std(roller)
    c = 0
    for item in arr:
        if item < m + s and item > m - 2 * s:
            c += 1
    if c > len(arr) * 0.50 and s < 200:
        return True
    return False

def get_rgb(temperature):
    if temperature < 200:
        r, g, b = 0, 0, 255
    elif temperature < 273:
        r = 0
        g = int(255 * (temperature - 200) / 73)
        b = 255
    elif temperature < 400:
        r = 0
        g = 255
        b = int(255 * (400 - temperature) / 127)
    elif temperature < 600:
        r = int(255 * (temperature - 400) / 199)
        g = 255
        b = 0
    else:
        r = 255
        g = max(int(255 * (1000 - temperature) / 400), 0)
        b = 0
    return (r, g, b)

def rgbtohex(rgb):
    r, g, b = rgb
    return f"#{r:02x}{g:02x}{b:02x}"

def calculate_surface_temperature(sigma, flux, co2, ch4, radius, pressure, density, heat):
    temp = ((flux * (1.0 - co2 - ch4)) / (16*np.pi*sigma)) ** 0.25
    temperature = temp * (1.0 + ((co2 + ch4) * pressure * 1000) / (density * heat))
    eqtemp = (flux/(16*np.pi*sigma))**0.25
    return eqtemp, temperature

# Tabs
tab1, tab2, tab3, tab4 = st.tabs(["📊 Light Curve Analysis", "🌡️ Temperature Simulation", "❓ Help", "⚙️ Settings"])

# TAB 1: LIGHT CURVE ANALYSIS
with tab1:
    col1, col2, col3 = st.columns([1, 1, 2])
    
    with col1:
        st.subheader("Download & Load")
        
        # Download section
        st.markdown("**Download from MAST**")
        pg_id = st.text_input("Product Group ID:", key="pg_id")
        if st.button("Download and Extract"):
            if len(pg_id) > 2:
                try:
                    # Create progress bar and status text
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    
                    # Download phase
                    status_text.text("📥 Downloading from MAST...")
                    progress_bar.progress(10)
                    
                    http_response = urlopen(f'https://mast.stsci.edu/api/v0.1/Download/bundle.zip?previews=false&obsid={pg_id}')
                    
                    # Read with progress tracking
                    file_size = int(http_response.headers.get('Content-Length', 0))
                    chunk_size = 8192
                    downloaded = 0
                    chunks = []
                    
                    while True:
                        chunk = http_response.read(chunk_size)
                        if not chunk:
                            break
                        chunks.append(chunk)
                        downloaded += len(chunk)
                        
                        if file_size > 0:
                            # Update progress (10-70% for download)
                            progress = int(10 + (downloaded / file_size * 60))
                            progress_bar.progress(min(progress, 70))
                    
                    file_data = b''.join(chunks)
                    progress_bar.progress(70)
                    
                    # Extract phase
                    status_text.text("📂 Extracting files...")
                    progress_bar.progress(80)
                    
                    zipfile = ZipFile(BytesIO(file_data))
                    download_path = os.path.join(os.getcwd(), "downloads")
                    os.makedirs(download_path, exist_ok=True)
                    zipfile.extractall(path=download_path)
                    
                    progress_bar.progress(100)
                    status_text.text("✅ Download complete!")
                    st.success(f"Files extracted to: {download_path}")
                    
                except Exception as e:
                    st.error(f"Download failed: {str(e)}")
            else:
                st.error("Please enter a valid Product Group ID")
        
        st.markdown("---")
        
        # File upload
        st.markdown("**Load FITS File**")
        uploaded_file = st.file_uploader("Choose a FITS file", type=['fits'])
        is_lc_file = st.checkbox("LC File (uncheck for TPF)", value=False)
        auto_save = st.checkbox("Auto Save figures", value=False)
        
        if uploaded_file is not None:
            try:
                # Save uploaded file temporarily
                with tempfile.NamedTemporaryFile(delete=False, suffix='.fits') as tmp_file:
                    tmp_file.write(uploaded_file.getvalue())
                    tmp_path = tmp_file.name
                
                if is_lc_file:
                    st.session_state.lc = lk.read(tmp_path)
                    st.session_state.obj_name = st.session_state.lc.meta["OBJECT"]
                    st.session_state.stellar_radius = st.session_state.lc.meta.get("RADIUS", 1.0)
                    st.session_state.star_temp = st.session_state.lc.meta.get("TEFF", 5778)
                else:
                    tpf = TessTargetPixelFile(tmp_path)
                    st.session_state.obj_name = tpf.get_header()["OBJECT"]
                    st.session_state.stellar_radius = tpf.get_header().get("RADIUS", 1.0)
                    st.session_state.star_temp = tpf.get_header().get("TEFF", 5778)
                    st.session_state.lc = tpf.to_lightcurve(aperture_mask=tpf.pipeline_mask)
                
                st.success(f"Loaded: {st.session_state.obj_name}")
                
                # Clean up temp file
                os.unlink(tmp_path)
                
            except Exception as e:
                st.error(f"Error loading file: {str(e)}")
    
    with col2:
        st.subheader("Analysis Parameters")
        
        if st.session_state.obj_name:
            st.info(f"**Current Object:** {st.session_state.obj_name}")
            st.write(f"**Stellar Radius:** {st.session_state.stellar_radius:.3f} R☉")
            st.write(f"**Stellar Temp:** {st.session_state.star_temp:.0f} K")
        
        stellar_mass = st.number_input("Stellar Mass (M☉)", min_value=0.1, max_value=100.0, value=1.0, step=0.1)
        
        if st.button("🔍 Analyze Light Curve", type="primary"):
            if st.session_state.lc is None:
                st.error("Please load a file first!")
            else:
                with st.spinner("Analyzing light curve..."):
                    try:
                        # Period calculations
                        upper_bound = 7
                        periodx = np.linspace(1, upper_bound, 10000)
                        bls = st.session_state.lc.to_periodogram(method="bls", period=periodx, frequency_factor=500)
                        st.session_state.period = bls.period_at_max_power
                        transit = bls.transit_time_at_max_power
                        
                        # Folding
                        phase_folded = st.session_state.lc.fold(period=st.session_state.period, epoch_time=transit)["flux"].value
                        moving_avg = rolling_avg(phase_folded)
                        st.session_state.flux = np.mean(st.session_state.lc.flux)
                        
                        # Check if planet detected
                        if within(phase_folded, moving_avg):
                            period_years = st.session_state.period.value / 365
                            st.session_state.planet_distance = (stellar_mass * period_years ** 2) ** (1 / 3)
                            
                            # Calculate planet radius
                            mean, minx = np.mean(moving_avg), np.min(moving_avg)
                            transit_depth = abs(mean - minx)
                            st.session_state.planet_radius = np.sqrt(transit_depth / st.session_state.flux.value) * st.session_state.stellar_radius * 109.076
                            
                            # Calculate stellar irradiance
                            star_rad = st.session_state.stellar_radius * 6.9e8
                            area = 4 * np.pi * star_rad**2
                            st.session_state.irradiance = sigma * area * st.session_state.star_temp**4
                            
                            st.success("✅ Potential exoplanet detected!")
                        else:
                            st.session_state.planet_radius = None
                            st.session_state.planet_distance = None
                            st.warning("⚠️ No planet detected in this light curve")
                        
                        # Save figures if requested
                        if auto_save and st.session_state.obj_name:
                            os.makedirs(st.session_state.obj_name, exist_ok=True)
                            
                    except Exception as e:
                        st.error(f"Analysis error: {str(e)}")
        
        # Display results
        if st.session_state.period is not None:
            st.markdown("---")
            st.markdown("**Results**")
            st.write(f"**Period:** {st.session_state.period.value:.2f} days")
            st.write(f"**Flux:** {st.session_state.flux:.2f}")
            
            if st.session_state.planet_radius is not None:
                st.write(f"**Planet Radius:** {st.session_state.planet_radius:.4f} R⊕")
                st.write(f"**Semi-major Axis:** {st.session_state.planet_distance:.4f} AU")
    
    with col3:
        st.subheader("Visualizations")
        
        if st.session_state.lc is not None:
            try:
                # Light curve plot
                fig1, ax1 = plt.subplots(figsize=(8, 4))
                st.session_state.lc.plot(ax=ax1)
                ax1.set_title("Light Curve")
                st.pyplot(fig1)
                plt.close()
                
                # Phase folded plot
                if st.session_state.period is not None:
                    fig2, ax2 = plt.subplots(figsize=(8, 4))
                    st.session_state.lc.fold(period=st.session_state.period).scatter(ax=ax2)
                    ax2.set_title("Phase Folded Light Curve")
                    st.pyplot(fig2)
                    plt.close()
                    
            except Exception as e:
                st.error(f"Plotting error: {str(e)}")

# TAB 2: TEMPERATURE SIMULATION
with tab2:
    st.subheader("🌡️ Exoplanet Temperature Simulator")
    st.info("Earth's reference values: Distance=1 AU, Radius=1 R⊕, Pressure=100 kPa, Density=1.2 kg/m³, Heat=700 J/kg·K, CO₂=400 ppm, CH₄=1.8 ppm, Solar Luminosity=3.8×10²⁶ W")
    
    col1, col2, col3 = st.columns([1, 1, 1])
    
    with col1:
        st.markdown("**Planetary Parameters**")
        distance = st.number_input("Distance (AU)", min_value=0.01, value=1.0, step=0.01)
        radius = st.number_input("Radius (Earth radii)", min_value=0.1, value=1.0, step=0.1)
        pressure = st.number_input("Pressure (kPa)", min_value=0.1, value=100.0, step=1.0)
        density = st.number_input("Density (kg/m³)", min_value=0.1, value=1.2, step=0.1)
    
    with col2:
        st.markdown("**Atmospheric Parameters**")
        specific_heat = st.number_input("Specific Heat (J/kg·K)", min_value=1.0, value=700.0, step=10.0)
        co2 = st.number_input("CO₂ Concentration (ppm)", min_value=0.0, value=400.0, step=10.0)
        ch4 = st.number_input("CH₄ Concentration (ppm)", min_value=0.0, value=1.8, step=0.1)
        irr = st.number_input("Star Luminosity (W)", min_value=1e20, value=3.8e26, format="%.2e", step=1e25)
    
    col_btn1, col_btn2 = st.columns(2)
    
    with col_btn1:
        if st.button("🔄 Transfer Data from Analysis", use_container_width=True):
            if st.session_state.planet_distance is not None and st.session_state.planet_radius is not None:
                st.session_state.temp_distance = float(st.session_state.planet_distance)
                st.session_state.temp_radius = float(st.session_state.planet_radius)
                st.session_state.temp_irr = float(st.session_state.irradiance)
                st.rerun()
            else:
                st.error("No planet data available to transfer!")
    
    with col_btn2:
        calculate = st.button("🧮 Calculate Temperature", type="primary", use_container_width=True)
    
    if calculate:
        try:
            # Calculate flux
            distance_m = distance * 1.49e11
            radius_m = radius * 6730000
            F = irr / (distance_m ** 2)
            
            # Greenhouse parameters
            CO2Parameter = co2 / 1000000 * 5.35
            CH4Parameter = ch4 / 1000000 * 0.036
            
            # Calculate temperatures
            eq_temp, surface_temp = calculate_surface_temperature(
                sigma, F, CO2Parameter, CH4Parameter, radius_m, pressure, density, specific_heat
            )
            
            # Calculate ESI
            esi = (1 - abs((288 - surface_temp) / (288 + surface_temp))) * \
                  (1 - abs((1 - distance) / (1 + distance)) ** 3) * \
                  (1 - abs((1 - radius) / (1 + radius)) ** 3)
            
            with col3:
                st.markdown("**Results**")
                st.metric("Equilibrium Temperature", f"{eq_temp:.2f} K")
                st.metric("Surface Temperature", f"{surface_temp:.2f} K")
                st.metric("Earth Similarity Index", f"{esi:.3f}")
                
                # Temperature visualization
                color = rgbtohex(get_rgb(surface_temp))
                st.markdown(f"""
                <div style='width: 150px; height: 150px; border-radius: 50%; 
                background-color: {color}; margin: 20px auto; border: 2px solid white;'>
                </div>
                """, unsafe_allow_html=True)
            
            # Comparison plot
            st.markdown("---")
            st.subheader("Comparison with Earth")
            
            parameters = ['Distance', 'Radius', 'Pressure', 'Density', 'Sp. Heat', 'CO₂', 'CH₄']
            exoplanet_values = [distance, radius, pressure, density, specific_heat, co2, ch4]
            earth_values = [1.0, 1.0, 100.0, 1.2, 700.0, 400.0, 1.8]
            normalized = [e/earth for e, earth in zip(exoplanet_values, earth_values)]
            
            fig, ax = plt.subplots(figsize=(10, 5))
            x_pos = np.arange(len(parameters))
            width = 0.35
            
            ax.bar(x_pos - width/2, normalized, width, label='Exoplanet (normalized)', color='steelblue')
            ax.bar(x_pos + width/2, [1.0]*len(parameters), width, label='Earth', color='green', alpha=0.7)
            
            ax.set_ylabel('Normalized Values')
            ax.set_title('Exoplanet vs Earth Parameters')
            ax.set_xticks(x_pos)
            ax.set_xticklabels(parameters, rotation=45, ha='right')
            ax.legend()
            ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
            
            st.pyplot(fig)
            plt.close()
            
        except Exception as e:
            st.error(f"Calculation error: {str(e)}")

# TAB 3: HELP
with tab3:
    st.subheader("📚 Help & Documentation")
    
    with st.expander("🔍 Finding and Downloading TESS Light Curves"):
        st.markdown("""
        1. Go to the **MAST archive** (link below)
        2. Select the collection **'MAST Catalogs'**
        3. Select the mission **'TESS Ctl v8.01'**
        4. Click on **'Advanced Search'** and restrict your search to close stars (ca <30pc)
        5. Look through all the TIC observations by:
           - Selecting 'MAST Observations by Object Name or RA/Dec'
           - Searching for TIC + the TIC IDs from the first list
        6. Observations with a lightcurve file will have a symbol at the beginning of their row
        7. Click **'Show Details'** → **'Details'** → Copy the **Product Group ID**
        8. Paste the Product Group ID in the download box and wait for the download
        """)
    
    with st.expander("📊 Light Curve Analysis"):
        st.markdown("""
        **How to Use:**
        1. **Load File:** Upload a TESS/KEPLER observation FITS file
           - Check 'LC File' for light curve files
           - Uncheck for 'TPF' (Target Pixel File)
        
        2. **Analyze:** Click 'Analyze Light Curve' to detect planets
           - The tool estimates orbital period
           - Identifies potential exoplanets through transit detection
           - Calculates planet radius and semi-major axis
        
        3. **Stellar Mass:** Not included in LC files by default
           - Default is 1 Solar Mass if left blank
           - Find accurate values on the MAST website
        
        4. **Auto Save:** Check to save figures as PNG files
           - Saved in the working directory
           - Organized by object name
        """)
    
    with st.expander("🌡️ Temperature Simulator"):
        st.markdown("""
        **Parameters:**
        - **Distance (AU):** Semi-major axis from host star
        - **Radius (Earth radii):** Planet size compared to Earth
        - **Pressure (kPa):** Atmospheric pressure
        - **Density (kg/m³):** Atmospheric density
        - **Specific Heat (J/kg·K):** Atmospheric specific heat capacity
        - **CO₂ Concentration (ppm):** Carbon dioxide in atmosphere
        - **CH₄ Concentration (ppm):** Methane in atmosphere
        - **Star Luminosity (W):** Total energy output of host star
        
        **Results:**
        - **Equilibrium Temperature:** Blackbody temperature based on stellar flux
        - **Surface Temperature:** Accounts for atmospheric greenhouse effects
        - **ESI:** Earth Similarity Index (0-1 scale)
        
        **Transfer Data:** Automatically fills parameters from light curve analysis
        """)
    
    with st.expander("ℹ️ Credits & References"):
        st.markdown("""
        **GUI and Simulation:**
        - Marco Leonardi, 2023
        - University of Bologna
        - marcoleonarditredici@gmail.com
        
        **Lightkurve Library:**
        - Lightkurve Collaboration
        - Scientists from STScI, NASA, and ESA
        - [ADS Reference](https://ui.adsabs.harvard.edu/abs/2018ascl.soft12013L/abstract)
        
        **Streamlit Version:** Adapted from original Tkinter GUI
        """)
    
    st.markdown("---")
    st.markdown("**Useful Links:**")
    st.markdown("- [MAST Archive](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html)")
    st.markdown("- [Lightkurve Documentation](https://docs.lightkurve.org/index.html)")

# TAB 4: SETTINGS
with tab4:
    st.subheader("⚙️ Settings")
    
    st.markdown("**Theme**")
    st.info("Use the theme switcher in Streamlit's menu (☰ → Settings → Theme)")
    
    st.markdown("---")
    st.markdown("**About This Application**")
    st.write("Version: 2.0 (Streamlit)")
    st.write("Original: Tkinter GUI by Marco Leonardi")
    st.write("Purpose: Exoplanet detection and characterization from light curve data")
