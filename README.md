# Streamlit Web App for Light Curve Analysis of TESS // KEPLER Observations
<img width="1395" height="703" alt="Screenshot 2025-10-06 at 23 08 48" src="https://github.com/user-attachments/assets/b8e3b27b-c796-4bd3-8dae-2cf82218c284" />

# Exoplanet Light Curve Analysis Tool

A modern web-based Python application for analyzing TESS and Kepler light curves to detect exoplanets and estimate their physical properties, including temperature modeling.

![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![Framework](https://img.shields.io/badge/framework-Streamlit-red.svg)
![Dependencies](https://img.shields.io/badge/dependencies-lightkurve%2C%20streamlit%2C%20matplotlib-orange.svg)

## 🌟 Features

- **Modern Web Interface**: Clean, responsive Streamlit-based UI accessible from any browser
- **Light Curve Analysis**: Upload and analyze TESS/Kepler observation files
- **Automatic Transit Detection**: Uses Box Least Squares (BLS) algorithm to identify potential exoplanets
- **Progress Tracking**: Real-time progress bars for downloads and analysis
- **Period Estimation**: Calculates orbital periods and folded light curves
- **Physical Parameter Estimation**: Estimates planet radius, semi-major axis, and orbital distance
- **Temperature Modeling**: Simulates exoplanet surface temperatures based on atmospheric parameters
- **Earth Similarity Index (ESI)**: Calculates habitability metrics
- **Interactive Visualizations**: Dynamic plots with matplotlib integration
- **MAST Archive Integration**: Direct download capability for TESS observations with progress tracking
- **Tab-based Navigation**: Organized workflow with separate analysis, simulation, help, and settings tabs
- **Session State Management**: Data persists across interactions

## 📋 Requirements

### Python Dependencies
```
streamlit
lightkurve
matplotlib
numpy
```

### System Requirements
- Python 3.8 or higher
- Internet connection (for downloading observations)
- ~100MB free disk space for temporary files
- Modern web browser (Chrome, Firefox, Safari, Edge)

## 🚀 Installation

1. **Clone the repository:**
```bash
git clone https://github.com/yourusername/exoplanet-analysis-tool.git
cd exoplanet-analysis-tool
```

2. **Install required packages:**
```bash
pip install streamlit lightkurve matplotlib numpy
```

3. **Run the application:**
```bash
streamlit run lightcurve_analysis.py
```

The app will automatically open in your default web browser at `http://localhost:8501`

## 🌐 Deployment

### Deploy to Streamlit Cloud (Free)
1. Push your code to GitHub
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Connect your GitHub repository
4. Click "Deploy"

Your app will be live with a shareable URL in minutes!

### Local Network Access
Share your app on your local network:
```bash
streamlit run lightcurve_analysis.py --server.address 0.0.0.0
```

## 📖 Usage Guide

### 1. Light Curve Analysis Tab

#### Loading Data
- **Option A - Upload File**: Click "Browse files" to upload a local FITS file
  - Drag and drop supported
  - Progress indicator during upload
- **Option B - MAST Download**: 
  1. Get Product Group ID from [MAST Archive](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html)
  2. Enter ID in the download field
  3. Click "Download and Extract"
  4. Real-time progress bar shows download and extraction status

#### File Types
- **Unchecked (TPF)**: Target Pixel Files - raw pixel data
- **Checked (LC)**: Light Curve files - processed photometry

#### Analysis Process
1. Upload or download your observation file
2. Set stellar mass (defaults to 1 solar mass if left blank)
3. Click "🔍 Analyze Light Curve"
4. Review results displayed in organized columns:
   - Orbital period
   - Planet radius (in Earth radii)
   - Semi-major axis (in AU)
   - Received flux
5. View interactive plots:
   - Raw light curve
   - Phase-folded light curve

#### Results Display
Results are organized into three columns:
- **Left**: Download and file loading controls
- **Center**: Analysis parameters and computed results
- **Right**: Visualization plots

### 2. Temperature Simulation Tab

#### Input Parameters
Organized into intuitive sections:

**Planetary Parameters:**
- **Distance (AU)**: Semi-major axis from host star
- **Radius (Earth radii)**: Planet size relative to Earth

**Atmospheric Parameters:**
- **Pressure (kPa)**: Atmospheric pressure
- **Density (kg/m³)**: Atmospheric density
- **Specific Heat (J/kg·K)**: Atmospheric heat capacity
- **CO2 Concentration (ppm)**: Carbon dioxide levels
- **CH4 Concentration (ppm)**: Methane levels
- **Star Luminosity (W)**: Host star's energy output

#### Interactive Features
- **Transfer Data Button**: Automatically populate simulation with discovered planet parameters
- **Calculate Button**: Compute temperatures and ESI
- **Visual Feedback**: Color-coded planet visualization based on temperature
- **Comparison Chart**: Side-by-side bar chart comparing exoplanet to Earth

#### Results Display
- Equilibrium Temperature (blackbody temperature)
- Surface Temperature (with greenhouse effects)
- Earth Similarity Index (ESI)
- Visual planet representation with temperature-based coloring
- Normalized parameter comparison chart

### 3. Earth Reference Values
Provided in the interface for easy reference:
- Distance: 1 AU
- Radius: 1 Earth radius
- Pressure: 100 kPa
- Density: 1.2 kg/m³
- Specific Heat: 700 J/kg·K
- CO2: 400 ppm
- CH4: 1.8 ppm
- Solar Luminosity: 3.8×10²⁶ W

### 4. Help Tab
Comprehensive documentation with expandable sections:
- Finding and downloading TESS light curves
- Light curve analysis instructions
- Temperature simulator guide
- Credits and references
- Quick links to MAST and Lightkurve documentation

### 5. Settings Tab
- Theme customization (via Streamlit's native theme switcher)
- Application information
- Version details

## 🔬 Scientific Background

### Transit Detection Method
The tool uses the Box Least Squares (BLS) periodogram to identify periodic dimming events in stellar light curves, which indicate planetary transits.

### Physical Parameter Calculations
- **Planet Radius**: Derived from transit depth and stellar radius using the formula:
  ```
  R_planet = √(transit_depth) × R_star × 109.076
  ```
- **Semi-major Axis**: Calculated using Kepler's third law:
  ```
  a = (M_star × P²)^(1/3)
  ```
- **Surface Temperature**: Uses energy balance equations with greenhouse effect modeling:
  ```
  T_eq = (F / (16πσ))^0.25
  T_surface = T_eq × (1 + greenhouse_factor)
  ```

### Earth Similarity Index (ESI)
Quantifies how Earth-like an exoplanet is based on:
- Surface temperature similarity
- Orbital distance comparison
- Size comparison

Formula incorporates weighted differences in each parameter.

## 📊 Output and Data Management

### Session State
The Streamlit version maintains data across interactions:
- Loaded light curves persist between tabs
- Calculated parameters available for temperature simulation
- Results remain visible after computation

### File Management
- Uploaded files are processed in temporary directories
- Downloads saved to `downloads/` subdirectory
- No manual cleanup required for temporary files

## 🔧 Technical Details

### Key Algorithms
- **BLS Periodogram**: For transit detection
- **Phase Folding**: For signal enhancement
- **Rolling Average**: For noise reduction (3-point moving average)
- **Stefan-Boltzmann Law**: For temperature calculations
- **Greenhouse Effect Modeling**: CO₂ and CH₄ radiative forcing

### Architecture Improvements
- **Reactive UI**: Streamlit's reactive framework eliminates manual threading
- **Session State**: Built-in state management for data persistence
- **Progress Bars**: Native support for download progress tracking
- **Auto-refresh**: UI updates automatically on user interaction

### Data Sources
- **TESS**: Transiting Exoplanet Survey Satellite
- **Kepler/K2**: NASA's planet-hunting missions
- **MAST Archive**: Mikulski Archive for Space Telescopes

## 🎨 Advantages Over Tkinter Version

### User Experience
- ✅ Modern, clean web interface
- ✅ Responsive design works on any screen size
- ✅ No installation required for end users (when deployed)
- ✅ Accessible from any device with a browser
- ✅ Native support for file drag-and-drop

### Development
- ✅ Simpler codebase (~20% less code)
- ✅ No manual threading required
- ✅ Built-in progress indicators
- ✅ Easier to maintain and extend
- ✅ Free cloud deployment options

### Deployment
- ✅ One-click deployment to Streamlit Cloud
- ✅ Shareable URL for collaboration
- ✅ Automatic HTTPS
- ✅ No server configuration needed

## 🤝 Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for:
- Bug reports
- Feature requests
- Code improvements
- Documentation updates
- UI/UX enhancements

### Development Setup
```bash
git clone https://github.com/yourusername/exoplanet-analysis-tool.git
cd exoplanet-analysis-tool
pip install -r requirements.txt
streamlit run lightcurve_analysis.py
```

## 📝 Version History

### Version 2.0 (Streamlit)
- Complete rewrite using Streamlit framework
- Modern web-based interface
- Real-time progress tracking for downloads
- Improved visualization layout
- Session state management
- Cloud deployment ready

### Version 1.0 (Tkinter)
- Original desktop GUI application
- Basic light curve analysis
- Temperature simulation
- Local file management

## 🙏 Acknowledgments

- **Lightkurve Collaboration**: For the excellent lightkurve library
- **NASA Kepler/K2 Guest Observer Office**: For mission support
- **Space Telescope Science Institute (STScI)**: For data archiving
- **University of Bologna**: Project affiliation
- **Streamlit Team**: For the amazing web framework

## 📜 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 📧 Contact

**Marco Leonardi**  
University of Bologna / University of Leiden  
marcoleonarditredici@gmail.com

## 🔗 Useful Links

- [MAST Archive](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html)
- [Lightkurve Documentation](https://docs.lightkurve.org/index.html)
- [Streamlit Documentation](https://docs.streamlit.io)
- [TESS Mission](https://tess.mit.edu/)
- [Kepler/K2 Mission](https://www.nasa.gov/mission_pages/kepler/main/index.html)

---

*Originally developed as a Tkinter GUI, now reimagined as a modern web application with Streamlit.*
