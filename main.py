import tkinter as tk
import lightkurve as lk
from lightkurve import TessTargetPixelFile
from tkinter import filedialog, ttk
import tkinter.messagebox as msgbox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt
from urllib.request import urlopen
from io import BytesIO
import threading
from zipfile import ZipFile
import webbrowser
import sv_ttk
import os
import ssl
import numpy as np

ssl._create_default_https_context = ssl._create_unverified_context

# -----------------------------------------------------------------------------------------------
# Tools for exoplanet research and temperature estimating of TESS & KEPLER light curves
# This is a python Graphical User Interface (GUI) that can display a TESS/KEPLER lightcurve using the 
# lightkurve python library {NASA Kepler/K2 Guest Observer Office and the Space Telescope Science Institute (STScI)}
# It features automatic planetary transit detection, estimate of period, distance, radius of the potential exoplanet, 
# and a secondary GUI to simulate planetary temperature given various input parameters, comparing these to Earth's.

# Marco Leonardi 2023, x esame astrobiologia
# -----------------------------------------------------------------------------------------------


# Constants for Reference
AU = 149.6e6  			 	# Astronomical Unit in kilometers
SolarConstant = 1361.0   	# Solar constant in W/m^2
sigma = 5.67e-8          	# Stefen Boltzmann Costant

# These are just for reference, not used in calculations
CO2Concentration = 0.0004   # Concentration of CO2 in Earth's atmosphere (400 ppm)
CH4Concentration = 1.8e-6	# Concentration of CH4 in Earth's atmosphere (1.8 ppm)

root = tk.Tk()
root.title("Lightcurve Analysis")

# This allows for the different tabs at the top
notebook = ttk.Notebook(root)
notebook.grid(row=1, column=0, columnspan=4, rowspan=11, sticky="nsew")

# Tab 1: Light Curve Analysis
lc_tab = ttk.Frame(notebook)
notebook.add(lc_tab, text='Light Curve Analysis')

# Tab 2: Temperature Simulation
tempsim_tab = ttk.Frame(notebook)
notebook.add(tempsim_tab, text='Temperature Simulation')

# Tab 3: Help
help_tab = ttk.Frame(notebook)
notebook.add(help_tab, text='Help')

# Tab 4: Settings
settings_tab = ttk.Frame(notebook)
notebook.add(settings_tab, text="Settings")


# -----------------------------------------------------------------------------------------------
# Functions used in the finder GUI
# -----------------------------------------------------------------------------------------------
def download_and_unzip():
	product_group_id = pg_entry.get()
	if len(product_group_id) > 2:
		print("Downloading...")
		http_response = urlopen('https://mast.stsci.edu/api/v0.1/Download/bundle.zip?previews=false&obsid=' + str(product_group_id))
		print("Extracting...")
		zipfile = ZipFile(BytesIO(http_response.read()))
		zipfile.extractall(path=f"{os.getcwd()}/downloads")
		print("Download Complete")
	else:
		msgbox.showinfo(title="No File Found", message="Please input a valid product ID from a TESS/KEPLER observation")

def callback(url):
	webbrowser.open(url)

def load_file():
	file_path = filedialog.askopenfilename(filetypes=[("Fits Files", "*.fits")])
	path_label.configure(text=file_path)
	try:
		tpf = TessTargetPixelFile(path_label["text"])
		obj = tpf.get_header()["OBJECT"]
		lc = tpf.to_lightcurve(aperture_mask=tpf.pipeline_mask)
	except:
		lc = lk.read(path_label["text"])
		obj = lc.meta["OBJECT"]
		
	current_label.configure(text=obj)

def analyse_all():
    if path_label["text"] == " ":
        current_label.configure(text="Please load a file")
        msgbox.showinfo(title="No File Found", message="Please select a valid file")
        return
    if varfile.get() == 0:
        tpf = TessTargetPixelFile(path_label["text"])
        obj = tpf.get_header()["OBJECT"]
        stellar_radius = tpf.get_header()["RADIUS"]
        stellar_radius_lab.configure(text=tpf.get_header()["RADIUS"])
        star_temp_lab.configure(text=tpf.get_header()["TEFF"])
        lc = tpf.to_lightcurve(aperture_mask=tpf.pipeline_mask)
        saved_name.configure(text=obj)
    else:
        lc = lk.read(path_label["text"])
        obj = lc.meta["OBJECT"]
        stellar_radius = lc.meta["RADIUS"]
        stellar_radius_lab.configure(text=lc.meta["RADIUS"])
        star_temp_lab.configure(text=lc.meta["TEFF"])

    # Period calculations
    upper_bound = 7
    periodx = np.linspace(1, upper_bound, 10000)
    bls = lc.to_periodogram(method="bls", period=periodx, frequency_factor=500)
    period = bls.period_at_max_power
    transit = bls.transit_time_at_max_power
    period_label.configure(text=f"Estimated Period: {period:.2f}")

    # Folding
    phase_folded = lc.fold(period=period, epoch_time=transit)["flux"].value
    moving_avg = rolling_avg(phase_folded)

    fig = Figure(figsize=(5, 4), dpi=100)
    fig2 = Figure(figsize=(5, 4), dpi=100)
    ax = fig.add_subplot(111)
    ax2 = fig2.add_subplot(111)

    # Plot all 3 figures
    lc.plot(ax=ax)
    lc.fold(period=period, epoch_time=transit).scatter(ax=ax2)

    # light curve canvas
    canvas = FigureCanvasTkAgg(fig, master=right_frame)
    canvas.get_tk_widget().delete("all")
    canvas.draw()
    canvas.get_tk_widget().grid(row=1, column=2)
    lc_title.configure(text="Light curve")
    
    # ph folding canvas
    canvas2 = FigureCanvasTkAgg(fig2, master=right_frame)
    canvas2.get_tk_widget().delete("all")
    canvas2.draw()
    canvas2.get_tk_widget().grid(row=3, column=2)
    phase_title.configure(text="Phase Folding")

    if saver.get() == 1:
        if not os.path.exists(obj):
            os.makedirs(obj)
        lightp = os.path.join(obj, f"{obj}_lightcurve.png")
        foldedp = os.path.join(obj, f"{obj}_folded.png")
        fig.savefig(lightp)
        fig2.savefig(foldedp)
        saved.configure(text=f"Saved as .png in {str(os.getcwd()) + '/' + str(obj)}")

    if within(phase_folded, moving_avg):
        if len(stellar_mass_entry.get()) > 0:
            star_mass = float(stellar_mass_entry.get())  # * 1.89 * 10**30
        else:
            star_mass = 1  # 1.89 * 10**30
        period_years = period.value / 365
        planet_distance = (star_mass * period_years ** 2) ** (1 / 3)
        flux = np.mean(lc.flux)
        flux_label.configure(text=f"Received Flux: {flux:.2f}")

        # Calculate planet radius using the transit depth and stellar radius
        mean, minx = np.mean(moving_avg), np.min(moving_avg)
        transit_depth2 = abs(mean - minx)
        planet_radius = np.sqrt(
            transit_depth2 / flux.value) * stellar_radius * 109.076  # 109.076 is the conversion sun-rad to earth-rad

        # Update the GUI with the calculated planet information
        planet_radius_lab.configure(text=f"Estimated Radius: {planet_radius:.4f} Earth radii")
        planet_distance_lab.configure(text=f"Estimated Semi-major axis: {planet_distance:.4f} AU")
        saved_radius.configure(text=planet_radius)
        saved_distance.configure(text=planet_distance)

    else:
        period_label.configure(text="No planet was found orbiting this star")
        planet_distance_lab.configure(text="")
        planet_radius_lab.configure(text="")
        flux_label.configure(text=f"Received Flux: {np.mean(lc.flux):.2f}")

def analyse_all_threaded():
    t1 = threading.Thread(target=analyse_all)
    t1.start()

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

# -----------------------------------------------------------------------------------------------
# Functions used in the Temperature Simulator
# -----------------------------------------------------------------------------------------------

def calculate_surface_temperature(sigma, flux, co2, ch4, radius, pressure, density, heat):
    temp = ((flux * (1.0 - co2 - ch4)) / (16*np.pi*sigma)) ** 0.25
    temperature = temp * (1.0 + ((co2 + ch4) * pressure * 1000) / (density * heat))
    eqtemp = (flux/(16*np.pi*sigma))**0.25
    return eqtemp, temperature

def get_rgb(temperature):
    # You can customize the color mapping based on your preference
    if temperature < 200:
        r = 0
        g = 0
        b = 255
    elif temperature < 273:
        # Blue to cyan gradient
        r = 0
        g = int(255 * (temperature - 200) / 73)
        b = 255
    elif temperature < 400:
        # Cyan to green gradient
        r = 0
        g = 255
        b = int(255 * (400 - temperature) / 127)
    elif temperature < 600:
        # Green to yellow gradient
        r = int(255 * (temperature - 400) / 199)
        g = 255
        b = 0
    else:
        # Yellow to red gradient
        r = 255
        g = max(int(255 * (1000 - temperature) / 400), 0)
        b = 0

    return (r, g, b)

def rgbtohex(rgb):
    r, g, b, = rgb
    return f"#{r:02x}{g:02x}{b:02x}"

def temp_sim():    
    def calculate_temperature():
        try:
            distance = float(distance_entry.get()) * 1.49e11
            radius = float(radius_entry.get()) * 6730000
            pressure = float(pressure_entry.get())
            density = float(density_entry.get())
            heat = float(specificHeat_entry.get())
            co2 = float(co2_entry.get())
            ch4 = float(ch4_entry.get())
            irr = float(irr_entry.get())

            F = irr / (distance ** 2)
            
            CO2Parameter = co2 / 1000000 * 5.35  # Scaling factor for CO2 greenhouse effect (Earth-like climate) - https://www.ipcc.ch/report/ar5/wg1/
            CH4Parameter = ch4 / 1000000 * 0.036  # Scaling factor for CH4 greenhouse effect (Earth-like climate) - https://www.ipcc.ch/report/ar5/wg1/
            eqt, temperature = calculate_surface_temperature(sigma, F, CO2Parameter, CH4Parameter, radius, pressure, density,
                                                        heat)

            # Esimtate ESI
            esi = (1 - abs((288 - temperature) / (288 + temperature))) * (
                        1 - abs((1 - distance/1.49e11) / (1 + distance/1.49e11)) ** 3) * (1 - abs((1 - radius/6730000) / (1 + radius/6730000)) ** 3)
            planet_ESI.configure(text=f"ESI: {esi:.2f}")

            # Update the result label
            result_eq.configure(text="Equilibrium Temperature {:.2f} °K".format(eqt))
            result_label.configure(text="Est. Surface Temperature: {:.2f} °K".format(temperature))
            saved_temperature.configure(text=temperature)
            # Update the circle size and color
            circle_size = min(int(radius / 673000), 95)
            circle_color = rgbtohex(get_rgb(temperature))
            canvas4.coords(circle, 125 - circle_size, 125 - circle_size, 125 + circle_size, 125 + circle_size)
            canvas4.itemconfigure(circle, width=circle_size, fill=circle_color)

            plot_comparison()

        except ValueError:
            msgbox.showinfo(title="Bad Parameters", message="Invalid input parameters")

    def plot_comparison():
        parameters = ['D', 'R', 'P', 'Rho', 'Q', 'CO2', 'CH4']
        exoplanet_values = [float(distance_entry.get()), float(radius_entry.get()), float(pressure_entry.get()),
                            float(density_entry.get()),
                            float(specificHeat_entry.get()), float(co2_entry.get()), float(ch4_entry.get())]
        earth_values = [1.0, 1.0, 100.0, 1.2, 700.0, 400.0, 1.8]
        ones = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        ex_values = [x/e for x,e in zip(exoplanet_values, earth_values)]

        fig, ax3 = plt.subplots(figsize=(3, 4))
        bar_width = 0.3  # Adjust the bar width

        # Set the positions of the bars
        x_pos = np.arange(len(parameters))

        ax3.bar(x_pos, ex_values, width=bar_width, label='Exoplanet')
        ax3.bar(x_pos + bar_width, ones, width=bar_width, label='Earth')

        ax3.set_ylabel('Values')
        ax3.set_title('Exoplanet vs Earth')
        ax3.set_xticks(x_pos + bar_width / 2)
        ax3.set_xticklabels(parameters, rotation=45, ha='right')  # Rotate x-axis labels by 45 degrees
        ax3.legend()

        # Clear previous content in the bar_canvas
        bar_canvas.delete('all')

        # Create a FigureCanvasTkAgg instance and display it in the bar_canvas
        canvas3 = FigureCanvasTkAgg(fig, master=bar_canvas)
        canvas3.draw()
        canvas3.get_tk_widget().grid(row=0, column=2)

        # Add a toolbar
        toolbar = NavigationToolbar2Tk(canvas3, bar_canvas)
        toolbar.update()
        canvas3.get_tk_widget().grid(row=0, column=2)

    def compare_to_earth():

        if period_label["text"] == "No planet was found orbiting this star" or period_label['text'] == "Estimated Period: ":
            msgbox.showinfo(title="No planet", message="No planet found to compare")
            return
        else:
            star_rad = float(stellar_radius_lab['text']) * 6.9e8
            star_temp = float(star_temp_lab['text']) 

            # Calculate the surface temperature
            area = 4* np.pi * star_rad**2
            
            irr = (sigma * area * star_temp**4) # / ((saved_distance['text'] * 1.49e11)**2)
            # Clear entries
            distance_entry.delete(0, tk.END)
            radius_entry.delete(0, tk.END)
            pressure_entry.delete(0, tk.END)
            density_entry.delete(0, tk.END)
            specificHeat_entry.delete(0, tk.END)
            co2_entry.delete(0, tk.END)
            ch4_entry.delete(0, tk.END)
            irr_entry.delete(0, tk.END)

            # Input found planet
            distance_entry.insert(0, float(saved_distance["text"]))
            radius_entry.insert(0, float(saved_radius["text"]))
            pressure_entry.insert(0, 100)
            density_entry.insert(0, 1.2)
            specificHeat_entry.insert(0, 700)
            co2_entry.insert(0, 400)
            ch4_entry.insert(0, 1.8)
            irr_entry.insert(0, irr)

            calculate_temperature()
    

    # HEADER
    temp_frame = ttk.Frame(tempsim_tab)
    temp_frame.grid(row=0, column=3)
    temp_tab = ttk.Frame(tempsim_tab)
    temp_tab.grid(row=0, column=0, sticky="n")
    sim_header = ttk.Label(temp_tab, text="The 'help' tab contains Earth's values for reference.")
    sim_header.grid(row=0, column=0, sticky="n")
    
    
    # Result label
    result_label = ttk.Label(temp_tab, text=" ", font=("Consolas", 14))
    result_label.grid(row=11, column=0, columnspan=2)
    result_eq = ttk.Label(temp_tab, text=" ", font=("Consolas", 14))
    result_eq.grid(row=12, column=0, columnspan=2)

    # Planet ESI
    planet_ESI = ttk.Label(temp_tab, text="ESI: ", font=("Consolas", 14))
    planet_ESI.grid(row=13, column=0, columnspan=2)

    # Distance input
    distance_label = ttk.Label(temp_tab, text="Distance (AU):")
    distance_label.grid(row=1, column=0)
    distance_entry = ttk.Entry(temp_tab)
    distance_entry.grid(row=1, column=1)

    # Radius input
    radius_label = ttk.Label(temp_tab, text="Radius (Earth radii):")
    radius_label.grid(row=2, column=0)
    radius_entry = ttk.Entry(temp_tab)
    radius_entry.grid(row=2, column=1)

    # Pressure input
    pressure_label = ttk.Label(temp_tab, text="Pressure (kPa)")
    pressure_label.grid(row=3, column=0)
    pressure_entry = ttk.Entry(temp_tab)
    pressure_entry.grid(row=3, column=1)

    # Density input
    density_label = ttk.Label(temp_tab, text="Density (kg/m^-3)")
    density_label.grid(row=4, column=0)
    density_entry = ttk.Entry(temp_tab)
    density_entry.grid(row=4, column=1)

    # Heat input
    heat_label = ttk.Label(temp_tab, text="Specific Heat (J/Kg*K)")
    heat_label.grid(row=5, column=0)
    specificHeat_entry = ttk.Entry(temp_tab)
    specificHeat_entry.grid(row=5, column=1)

    # CO2 input
    co2_label = ttk.Label(temp_tab, text="Co2 Concentration (ppm \\ mg/L)")
    co2_label.grid(row=6, column=0)
    co2_entry = ttk.Entry(temp_tab)
    co2_entry.grid(row=6, column=1)

    # CH4 input
    ch4_label = ttk.Label(temp_tab, text="Methane Concentration (ppm \\ mg/L)")
    ch4_label.grid(row=7, column=0)
    ch4_entry = ttk.Entry(temp_tab)
    ch4_entry.grid(row=7, column=1)

    # Irradiance input
    irr_label = ttk.Label(temp_tab, text="Star Luminosity (W)")
    irr_label.grid(row=8, column=0)
    irr_entry = ttk.Entry(temp_tab)
    irr_entry.grid(row=8, column=1)

    # Calculate button
    calculate_button = ttk.Button(temp_tab, text="Calculate", command=calculate_temperature)
    calculate_button.grid(row=9, column=0, columnspan=2)

    # Other button
    earth_comp = ttk.Button(temp_tab, text="Transfer Data", command=compare_to_earth)
    earth_comp.grid(row=10, column=0, columnspan=2)

    # Circle widget
    canvas4 = tk.Canvas(temp_frame, width=250, height=250)
    circle = canvas4.create_oval(0, 0, 250, 250, outline="")
    canvas4.grid(row=1, column=2, columnspan=1)

    # Bar graph canvas
    bar_canvas = tk.Canvas(temp_frame, width=250, height=250)
    bar_canvas.grid(row=0, column=2, columnspan=1)



# -----------------------------------------------------------------------------------------------
# Help tab:
# -----------------------------------------------------------------------------------------------
def display_help_window():
    def display_general_help():
        general_help_text = """Finding and Downloading TESS Light Curve (LC) files:
        1. Go to the MAST archive, link in the bottom right
        2. Select the collection 'Mast Catalogs'
        3. Select the mission 'TESS Ctl v8.01'
        4. Click on 'Advanced Search' and restrict your search to close stars (ca <30pc)
        5. Look through all the TIC observations by:
         - Selecting 'MAST Observations by Object Name or RA/Dec'
         - Searching for TIC+ the TIC IDs of the first list
        6. Observations with a lightcurve file will have a symbol at the beginning of their row
        7. Click 'Show Details' -> 'Details' -> 'Copy the Product Group ID'
        8. Paste the Product Group ID in the download box and wait for the download."""

        msgbox.showinfo(title="Download Help", message=general_help_text)

    def display_lc_help():
        lc_help_text = """Light Curve Analysis Help

    1. Load File:
       - You can load a TESS/KEPLER observation FITS file or a lightkurve file.
       - Click 'Load File' to choose the file from your local directory.
       - Check 'LC File' for light curve files, uncheck for 'TPF' pixel files. 

    2. Analyze Light Curve:
       - After loading a file, click 'Analyse Light Curve' to estimate the period and identify potential exoplanets.
       - The light curve, periodogram, folded light curve, and other information will be displayed
       - Stellar Mass is not included in these LC files, so you may have to tweak it manually. The values are available on the MAST website.
       - Default stellar mass (if left blank) is 1 Solar Mass).

    3. Auto Save figures:
       - If you check this box, the figures will be automatically saved as PNG files in the working directory.
       - This will not override the saved files for different observations, but it will for the same ones."""

        msgbox.showinfo(title="Light Curve Analysis Help", message=lc_help_text)

    def display_temp_sim_help():
        temp_sim_help_text = """Exoplanet Temperature Simulator Help

    1. Distance (AU):
       - Enter the semi-major axis distance of the exoplanet from its host star in Astronomical Units (AU).

    2. Radius (Earth radii):
       - Enter the radius of the exoplanet in Earth radii.

    3. Pressure (kPa):
       - Enter the atmospheric pressure of the exoplanet in kilopascals (kPa).

    4. Density (kg/m^-3):
       - Enter the atmospheric density of the exoplanet in kilograms per cubic meter (kg/m^-3).

    5. Specific Heat (J/Kg*K):
       - Enter the specific heat of the exoplanet's atmosphere in Joules per kilogram per Kelvin (J/Kg*K).

    6. CO2 Concentration (ppm):
       - Enter the concentration of carbon dioxide (CO2) in the exoplanet's atmosphere in parts per million (ppm).

    7. Methane Concentration (ppm):
       - Enter the concentration of methane (CH4) in the exoplanet's atmosphere in parts per million (ppm).

    8. Star Luminosity (Watts):
       - Enter the luminosity of the host star in Watts.

    9. Calculate:
       - Click 'Calculate' to estimate the surface temperature of the exoplanet.

    10. Transfer Data:
       - Click 'Transfer Data' to automatically fill the parameters with the estimated values from the light curve analysis, 
       and fill the remaining atmospheric parameters with Earth's.
    
    11. Results: 
	   - Equilibrium Temperature is the black body temperature of a planet a distance D from the host star; it depends only on 
	   the flux of the star and the planet's distance. 
	   - The esimated surface temperature takes into account your atmospheric parameters inputs. """

        msgbox.showinfo(title="Temperature Simulator Help", message=temp_sim_help_text)

    def cred():
        creds = """ GUI and Sim:
		Marco Leonardi 2023 
		University of Bologna
		marcoleonarditredici@gmail.com
		------------------------------
		The lightkurve library was made by the Lightkurve
		Collaboration group that includes scientists from the STScI, NASA and ESA
		https://ui.adsabs.harvard.edu/abs/2018ascl.soft12013L/abstract"""
        msgbox.showinfo(title="Credits", message=creds)

    # General Help Button
    general_help_button = ttk.Button(help_tab, text="Finding and Downloading LCs", command=display_general_help)
    general_help_button.grid(row=0, column=1, padx = 10, pady = 10)

    # Light Curve Analysis Help Button
    lc_help_button = ttk.Button(help_tab, text="Light Curve Analysis Help", command=display_lc_help)
    lc_help_button.grid(row=1, column=1, padx = 10, pady = 10)

    # Temperature Simulator Help Button
    temp_sim_help_button = ttk.Button(help_tab, text="Temperature Simulator Help", command=display_temp_sim_help)
    temp_sim_help_button.grid(row=2, column=1, padx = 10, pady = 10)

    # Credits
    credits_button = ttk.Button(help_tab, text="Credits", command=cred)
    credits_button.grid(row=3, column=1, padx = 10, pady = 10)


# -----------------------------------------------------------------------------------------------
# Settings / Themes:
# -----------------------------------------------------------------------------------------------

def settings():	

	def on_closing():
		result = msgbox.askyesno(title="Quit", message="Are you sure you want to exit?")
		if result:
			root.destroy()

        
	mode_switch = ttk.Button(settings_tab, text="Change Theme", command=sv_ttk.toggle_theme)
	mode_switch.grid(row=0, column=0, padx=10, pady=10)
	root.protocol("WM_DELETE_WINDOW", on_closing)


# -----------------------------------------------------------------------------------------------
# Threading: (not doing anything atm, i removed it; not sure if it helps)
# -----------------------------------------------------------------------------------------------

def displayTabContents():
    t1 = threading.Thread(target=temp_sim)
    t2 = threading.Thread(target=display_help_window)
    t3 = threading.Thread(target=settings)

    t1.start()
    t2.start()
    t3.start()


# -----------------------------------------------------------------------------------------------
# GUI Section: buttons, labels and text boxes. Very messy
# -----------------------------------------------------------------------------------------------

left_frame = tk.Frame(lc_tab)
left_frame.grid(row=0, column=0)

# Download name entry
pg_label = ttk.Label(left_frame, text="Product Group ID:")
pg_label.grid(row=0, column=0)
pg_entry = ttk.Entry(left_frame)
pg_entry.grid(row=1, column=0)
# Path
path_label = ttk.Label(lc_tab, text=" ")
# Download button
download_button = ttk.Button(left_frame, text="Download and Extract", command=download_and_unzip)
download_button.grid(row=2, column=0)
# Load
load_button = ttk.Button(left_frame, text="Load File", command=load_file)
load_button.grid(row=3, column=0, pady=20)
# Current workspace
current_label = ttk.Label(left_frame, text="No file selected.")
current_label.grid(row=4, column=0)
# File type
varfile = tk.IntVar()
checkbox = ttk.Checkbutton(left_frame, text="LC File", variable=varfile, onvalue=1, offvalue=0)
checkbox.grid(row=5, column=0)
# Save figures
saver = tk.IntVar()
save_button = ttk.Checkbutton(left_frame, text="Auto Save figures", variable=saver, onvalue=1, offvalue=0)
save_button.grid(row=6, column=0)
saved = tk.Label(left_frame, text="")
saved.grid(row=7, column=0)
# Analyse
lc_button = ttk.Button(left_frame, text="Analyse Light Curve", command=analyse_all)
lc_button.grid(row=8, column=0, pady=5) 

# Result graph canvas
right_frame = tk.Frame(lc_tab)
right_frame.grid(row=0, column=2, rowspan=10, padx=50)

canvas = tk.Canvas(right_frame, bd=0, highlightthickness=0, width=250, height=250)
canvas.grid(row=1, column=2, rowspan=4, padx=10, pady=10)  
canvas2 = tk.Canvas(right_frame, bd=0, highlightthickness=0, width=250, height=250)
canvas2.grid(row=3, column=2, rowspan=4, padx=10, pady=10)
lc_title = ttk.Label(right_frame, text="")
lc_title.grid(row=0, column=2, sticky="n")
phase_title = ttk.Label(right_frame, text="")
phase_title.grid(row=2, column=2, sticky="n")

#Additional
central_frame = ttk.Frame(lc_tab)
central_frame.grid(row=0, column=1, padx=50, sticky="n")

stellar_mass = ttk.Label(central_frame, text="Stellar Mass (in M_s)")
stellar_mass.grid(row=0, column=1)
stellar_mass_entry = ttk.Entry(central_frame)
stellar_mass_entry.grid(row=1, column=1)
# Result stuff
period_label = ttk.Label(central_frame, text="Estimated Period: ")
period_label.grid(row=2, column=1, sticky="w")
planet_distance_lab = ttk.Label(central_frame, text="Estimated Semi-major axis: ")
planet_distance_lab.grid(row=3, column=1, sticky="w")
planet_radius_lab = ttk.Label(central_frame, text="Estimated Radius: ")
planet_radius_lab.grid(row=4, column=1, sticky="w")
flux_label = ttk.Label(central_frame, text="Received Flux: ")
flux_label.grid(row=5, column=1, sticky="w")

# Savers
saved_name = ttk.Label(lc_tab, text="")
saved_radius = ttk.Label(lc_tab, text="")
saved_distance = ttk.Label(lc_tab, text="")
saved_temperature = ttk.Label(lc_tab, text="")
saved_period = ttk.Label(lc_tab, text="")
saved_transit = ttk.Label(lc_tab, text="")
stellar_radius_lab = ttk.Label(lc_tab, text="")
star_temp_lab = ttk.Label(lc_tab, text="")

# Links
mast_link = ttk.Label(help_tab, text="MAST Archive", foreground="blue")
mast_link.grid(row=5, column=1, padx = 10, pady = 10)
mast_link.bind("<Button-1>", lambda e: callback("https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html"))
lk_link = ttk.Label(help_tab, text="Lightkurve", foreground="blue")
lk_link.grid(row=6, column=1, padx = 10, pady = 10)
lk_link.bind("<Button-1>", lambda e: callback("https://docs.lightkurve.org/index.html"))


def configure_resizing():
    # Set the weights for the rows and columns to allow resizing
    root.grid_rowconfigure(0, weight=1)
    root.grid_columnconfigure(0, weight=1)

    # Tab 1 - Light Curve Analysis
    lc_tab.grid_rowconfigure(1, weight=1)
    lc_tab.grid_rowconfigure(2, weight=1)
    lc_tab.grid_rowconfigure(3, weight=1)
    lc_tab.grid_rowconfigure(4, weight=1)
    lc_tab.grid_rowconfigure(5, weight=1)
    lc_tab.grid_rowconfigure(6, weight=1)
    lc_tab.grid_rowconfigure(7, weight=1)
    lc_tab.grid_rowconfigure(8, weight=1)
    lc_tab.grid_rowconfigure(9, weight=1)
    lc_tab.grid_rowconfigure(10, weight=1)
    lc_tab.grid_rowconfigure(11, weight=1)
    lc_tab.grid_rowconfigure(12, weight=1)
    lc_tab.grid_rowconfigure(13, weight=1)  # Added this line for the canvas

    lc_tab.grid_columnconfigure(0, weight=1)
    lc_tab.grid_columnconfigure(1, weight=1)
    lc_tab.grid_columnconfigure(2, weight=1)

    # Configure the canvas inside the first tab to expand with the window
    canvas.grid(row=1, column=2, rowspan=4, sticky="nsew", padx=10, pady=10)
    canvas2.grid(row=3, column=2, rowspan=4, sticky="nsew", padx=10, pady=10)
   
def on_resize(event):
    configure_resizing()

displayTabContents()

sv_ttk.set_theme("dark")

lc_tab.bind("<Configure>", on_resize)

root.mainloop()
