import tkinter as tk
import lightkurve as lk
from lightkurve import TessTargetPixelFile
from tkinter import filedialog, ttk
import tkinter.messagebox as msgbox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt
from astropy.timeseries import BoxLeastSquares
from urllib.request import urlopen
from io import BytesIO
from zipfile import ZipFile
import webbrowser
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

csv_filename = None


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
		star_temp = tpf.get_header()["TEFF"]
		lc = tpf.to_lightcurve(aperture_mask=tpf.pipeline_mask)
		saved_name.configure(text=obj)
	else:
		lc = lk.read(path_label["text"])
		obj = lc.meta["OBJECT"]
		stellar_radius = lc.meta["RADIUS"]
		star_temp = lc.meta["TEFF"]
		
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

	fig = Figure(figsize=(15, 4), dpi=100)
	ax = fig.add_subplot(131)
	ax2 = fig.add_subplot(132)
	ax3 = fig.add_subplot(133)
	

	# Plot all 3 figures
	lc.plot(ax=ax)
	bls.plot(ax=ax2)
	lc.fold(period=period, epoch_time=transit).scatter(ax= ax3)

	# Create a Tkinter canvas for embedding the figure
	# . Lightcurve
	canvas = FigureCanvasTkAgg(fig, master=root)
	canvas.get_tk_widget().delete("all")
	canvas.draw()

	# Place the canvas on the GUI
	canvas.get_tk_widget().grid(row=4, column=0, columnspan=3)
	lc_title.configure(text="Light curve")
	period_title.configure(text="Periodogram")
	phase_title.configure(text="Phase Folding")
	
	if saver.get() == 1:
		if not os.path.exists(obj):
			os.makedirs(obj)
		name = os.path.join(obj, f"{obj}_lightcurve.png" )
		fig.savefig(name)
		saved.configure(text=f"Saved as .png in {str(os.getcwd()) + '/' + str(obj)}")
	
	if within(phase_folded, moving_avg):
		if len(stellar_mass_entry.get()) > 0 :
			star_mass = float(stellar_mass_entry.get()) # * 1.89 * 10**30
		else:
			star_mass = 1 # 1.89 * 10**30
		period_years = period.value/365
		planet_distance = (star_mass * period_years**2)**(1/3)
		flux = np.mean(lc.flux)
		flux_label.configure(text=f"Received Flux: {flux:.2f}")
		
		# Calculate planet radius using the transit depth and stellar radius
		mean, minx = np.mean(moving_avg), np.min(moving_avg)
		transit_depth2 = abs(mean-minx)
		planet_radius = np.sqrt(transit_depth2/flux.value) * stellar_radius * 109.076 #109.076 is the conversion sun-rad to earth-rad

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
			c+=1
	if c > len(arr) * 0.50 and s < 200:
		return True

# -----------------------------------------------------------------------------------------------
# Functions used in the Temperature Simulator
# -----------------------------------------------------------------------------------------------

def calculate_surface_temperature(sigma, flux, co2, ch4, radius, pressure, density, heat):
	temp = ((flux * (1.0 - co2 - ch4)) / (4.0 * radius**2 * sigma)) ** 0.25
	temperature = temp * (1.0 + ((co2 + ch4) * pressure * 1000) / (density * heat))
	return temperature
       
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
	
def add_row():
	msgbox.showinfo("Work in Progress")

def download_csv():
	msgbox.showinfo("Work in Progress")


# -----------------------------------------------------------------------------------------------
# General handling of the second window (temperature)
# -----------------------------------------------------------------------------------------------

def temp_sim():
	
	def calculate_temperature():
		try:
			distance = float(distance_entry.get())
			radius = float(radius_entry.get())
			pressure = float(pressure_entry.get())
			density = float(density_entry.get())
			heat = float(specificHeat_entry.get())
			co2 = float(co2_entry.get())
			ch4 = float(ch4_entry.get())
			irr = float(irr_entry.get())

			# Calculate the surface temperature
			F = irr / (distance **2)
			CO2Parameter = co2 / 1000000 * 5.35 # Scaling factor for CO2 greenhouse effect (Earth-like climate) - https://www.ipcc.ch/report/ar5/wg1/
			CH4Parameter = ch4 / 1000000 * 0.036 # Scaling factor for CH4 greenhouse effect (Earth-like climate) - https://www.ipcc.ch/report/ar5/wg1/
			temperature = calculate_surface_temperature(sigma, F, CO2Parameter, CH4Parameter, radius, pressure, density, heat)
			
			# Esimtate ESI
			esi = (1-abs( (288 - temperature) / (288 + temperature)) ) * (1-abs( (1 - distance) / (1 + distance))**3 ) * (1-abs( (1 - radius) / (1 + radius))**3 )
			planet_ESI.configure(text=f"ESI: {esi:.2f}")

			# Update the result label
			result_label.configure(text="Surface Temperature: {:.2f} °C".format(temperature-273.15))
			saved_temperature.configure(text=temperature)
			# Update the circle size and color
			circle_size = min(int(radius * 10), 95)
			circle_color = rgbtohex(get_rgb(temperature))
			canvas.coords(circle, 100 - circle_size, 100 - circle_size, 100 + circle_size, 100 + circle_size)
			canvas.itemconfigure(circle, width=circle_size, fill=circle_color) 
			
			plot_comparison()
				  
		except ValueError:
			msgbox.showinfo(title="Bad Parameters", message="Invalid input parameters")
	
	def plot_comparison():
		parameters = ['D','R', 'P', 'Rho', 'Q', 'CO2', 'CH4']
		exoplanet_values = [float(distance_entry.get()),float(radius_entry.get()), float(pressure_entry.get()), float(density_entry.get()),
							float(specificHeat_entry.get()), float(co2_entry.get()), float(ch4_entry.get())]
		earth_values = [1.0 ,1.0, 100.0, 1.2, 700.0, 400.0, 1.8]

		fig, ax = plt.subplots(figsize=(3,4))
		bar_width = 0.3  # Adjust the bar width

		# Set the positions of the bars
		x_pos = np.arange(len(parameters))

		ax.bar(x_pos, exoplanet_values, width=bar_width,label='Exoplanet')
		ax.bar(x_pos + bar_width, earth_values, width=bar_width, label='Earth')

		ax.set_ylabel('Values')
		ax.set_title('Exoplanet vs Earth')
		ax.set_xticks(x_pos + bar_width / 2)
		ax.set_xticklabels(parameters, rotation=45, ha='right')  # Rotate x-axis labels by 45 degrees
		ax.legend()

		# Clear previous content in the bar_canvas
		bar_canvas.delete('all')

		# Create a FigureCanvasTkAgg instance and display it in the bar_canvas
		canvas2 = FigureCanvasTkAgg(fig, master=bar_canvas)
		canvas2.draw()
		canvas2.get_tk_widget().grid(row=0, column=0)

		# Add a toolbar
		toolbar = NavigationToolbar2Tk(canvas2, bar_canvas)
		toolbar.update()
		canvas2.get_tk_widget().grid(row=0, column=0)
	
	def compare_to_earth():
		
		if saved_radius["text"] == "":
			msgbox.showinfo(title="No planet", message="No planet found to compare")
			return 
		else:
			#Clear entries
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
			irr_entry.insert(0, 1361)
			
			calculate_temperature()
		
		
	# Create the GUI
	window = tk.Tk()
	window.title("Exoplanet Temperature Simulator")

	# Result label
	result_label = tk.Label(window, text=" ", font=("Arial", 14))
	result_label.grid(row=11, column=0, columnspan=2)

	# Planet ESI
	planet_ESI = tk.Label(window, text="ESI: ", font=("Arial", 14))
	planet_ESI.grid(row=12, column=0, columnspan=2)


	# Additional Notes Box
	T = tk.Text(window, height=5, width=30)
	T.grid(row=0, column=0, columnspan=2, sticky=tk.W+tk.E)
	quote = """ Earth parameters in given units for reference: 
	P = 100, Density = 1.2, Sp.Heat=700, CO2 = 400, CH4 = 1.8
	Solar Irradiance: 1361 W
	Marco Leonardi - 05/2023 - Unibo"""
	T.insert(tk.END, quote)
	T.config(state="disabled")


	# Distance input
	distance_label = tk.Label(window, text="Distance (AU):")
	distance_label.grid(row=1, column=0)
	distance_entry = tk.Entry(window)
	distance_entry.grid(row=1, column=1)

	# Radius input
	radius_label = tk.Label(window, text="Radius (Earth radii):")
	radius_label.grid(row=2, column=0)
	radius_entry = tk.Entry(window)
	radius_entry.grid(row=2, column=1)

	# Pressure input
	pressure_label = tk.Label(window, text="Pressure (kPa)")
	pressure_label.grid(row=3, column=0)
	pressure_entry = tk.Entry(window)
	pressure_entry.grid(row=3, column=1)

	# Density input
	density_label = tk.Label(window, text="Density (kg/m^-3)")
	density_label.grid(row=4, column=0)
	density_entry = tk.Entry(window)
	density_entry.grid(row=4, column=1)

	# Heat input
	heat_label = tk.Label(window, text="Specific Heat (J/Kg*K)")
	heat_label.grid(row=5, column=0)
	specificHeat_entry = tk.Entry(window)
	specificHeat_entry.grid(row=5, column=1)

	# CO2 input
	co2_label = tk.Label(window, text="Co2 Concentration (ppm \\ mg/L)")
	co2_label.grid(row=6, column=0)
	co2_entry = tk.Entry(window)
	co2_entry.grid(row=6, column=1)

	# CH4 input
	ch4_label = tk.Label(window, text="Methane Concentration (ppm \\ mg/L)")
	ch4_label.grid(row=7, column=0)
	ch4_entry = tk.Entry(window)
	ch4_entry.grid(row=7, column=1)

	# Irradiance input
	irr_label = tk.Label(window, text="Star Irradiance (Watts)")
	irr_label.grid(row=8, column=0)
	irr_entry = tk.Entry(window)
	irr_entry.grid(row=8, column=1)

	# Calculate button
	calculate_button = tk.Button(window, text="Calculate", command=calculate_temperature)
	calculate_button.grid(row=9, column=0, columnspan=2)
	
	# Other button
	earth_comp = tk.Button(window, text="Transfer Data", command=compare_to_earth)
	earth_comp.grid(row=10, column=0, columnspan=2)

	# Circle widget
	canvas = tk.Canvas(window, width=250, height=250)
	circle = canvas.create_oval(0, 0, 250, 250, outline="")
	canvas.grid(row=13, column=0, columnspan=1, sticky="se")

	# Bar graph canvas
	bar_canvas = tk.Canvas(window, width=250, height=250)
	bar_canvas.grid(row=13, column=1, columnspan = 1)


# -----------------------------------------------------------------------------------------------
# Help window with instructions:
# -----------------------------------------------------------------------------------------------

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

8. Star Irradiance (Watts):
   - Enter the irradiance of the host star in Watts.

9. Calculate:
   - Click 'Calculate' to estimate the surface temperature of the exoplanet.

10. Transfer Data:
   - Click 'Transfer Data' to automatically fill the parameters with the estimated values from the light curve analysis, 
   and fill the remaining atmospheric parameters with Earth's."""

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

def display_help_window():
    help_window = tk.Toplevel(root)
    help_window.title("Help")

    # General Help Button
    general_help_button = tk.Button(help_window, text="Finding and Downloading LCs", command=display_general_help)
    general_help_button.grid(row=0, column=0, padx=10, pady=5)

    # Light Curve Analysis Help Button
    lc_help_button = tk.Button(help_window, text="Light Curve Analysis Help", command=display_lc_help)
    lc_help_button.grid(row=1, column=0, padx=10, pady=5)

    # Temperature Simulator Help Button
    temp_sim_help_button = tk.Button(help_window, text="Temperature Simulator Help", command=display_temp_sim_help)
    temp_sim_help_button.grid(row=2, column=0, padx=10, pady=5)
    
    # Credits
    credits_button = tk.Button(help_window, text="Credits", command=cred)
    credits_button.grid(row=3,column=0, padx=10, pady=5)



# -----------------------------------------------------------------------------------------------
# GUI Section: buttons, labels and text boxes
# -----------------------------------------------------------------------------------------------


root = tk.Tk()
root.title("Lightcurve Analysis")


# Download name entry
pg_label = tk.Label(root, text="Product Group ID:")
pg_label.grid(row=0, column=0)
pg_entry = tk.Entry(root)
pg_entry.grid(row=1, column=0)
#pb = ttk.Progressbar(root,orient ="horizontal",length = 200, mode ="indeterminate")
varfile = tk.IntVar()
#pb.grid(row=1, column=4, columnspan=2)
checkbox = tk.Checkbutton(root, text="LC File", variable = varfile, onvalue=1, offvalue=0)
checkbox.grid(row=2, column=1)
# Download button
download_button = tk.Button(root, text="Download and Extract", command=download_and_unzip)
download_button.grid(row=2, column=0)

# Links
mast_link = tk.Label(root, text="MAST Archive", fg="blue")
mast_link.grid(row=9, column=3)
mast_link.bind("<Button-1>", lambda e: callback("https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html"))

# Path
path_label = tk.Label(root, text=" ")

# Load
load_button = tk.Button(root, text="Load File", command=load_file)
load_button.grid(row = 1, column = 1)

# Current workspace
current_label = tk.Label(root, text="No file selected.")
current_label.grid(row=0, column=1)
stellar_mass = tk.Label(root, text="Stellar Mass (in M_s)")
stellar_mass.grid(row=1, column=3)
stellar_mass_entry = tk.Entry(root)
stellar_mass_entry.grid(row=2, column=3)

lc_button = tk.Button(root, text="Analyse Light Curve", command=analyse_all)
lc_button.grid(row=0, column=2)

help_button = tk.Button(root, text="Help", command=display_help_window)
help_button.grid(row=0, column=3)

# Result graph canvas
canvas = tk.Canvas(root, bd=0, highlightthickness=0, width=250, height=350)
canvas.grid(row=4, column=0, rowspan = 2, sticky="nsew")
lc_title = tk.Label(root, text="")
lc_title.grid(row=5, column=0, sticky="n")
period_title = tk.Label(root, text="")
period_title.grid(row=5, column=1, sticky="n")
phase_title = tk.Label(root, text="")
phase_title.grid(row=5, column=2, sticky="n")

# Periodgram
period_label = tk.Label(root, text="Estimated Period: ")
period_label.grid(row=6, column=0, sticky="w")

# Planet info
planet_distance_lab = tk.Label(root, text="Estimated Semi-major axis: ")
planet_distance_lab.grid(row=7, column=0, sticky="w")
planet_radius_lab = tk.Label(root, text="Estimated Radius: ")
planet_radius_lab.grid(row=8, column=0, sticky="w")
flux_label = tk.Label(root, text="Received Flux: ")
flux_label.grid(row=9, column=0, sticky="w")

# Save figures
saver = tk.IntVar()
save_button = tk.Checkbutton(root, text="Auto Save figures", variable=saver, onvalue=1, offvalue=0)
save_button.select()
save_button.grid(row=2, column=2)
saved = tk.Label(root, text="")
saved.grid(row=3, column=2)

# Button to open the temp calculator
open_temp_sim = tk.Button(root, text="Temperature Simulation", command=temp_sim)
open_temp_sim.grid(row=1, column=2, sticky="n")

# CSV file -- This feature is a work in progress. Doesn't currently work - going to add 
# the possibility to make a csv file of all the light curves you analyse and export it
 
#csv_row = tk.Button(root, text="Add to Table", command = add_row)
#csv_row.grid(row=3, column=3)
#csv_download = tk.Button(root, text="Download Table", command=download_csv)
#csv_download.grid(row=4, column=3)

saved_name = tk.Label(root, text="")
saved_radius = tk.Label(root, text="")
saved_distance= tk.Label(root, text="")
saved_temperature= tk.Label(root, text="")
saved_period= tk.Label(root, text="")
saved_transit = tk.Label(root, text="")e

def configure_resizing():
    # Allow rows and columns to resize proportionally when the window is resized
    root.grid_rowconfigure(1, weight=1)
    root.grid_rowconfigure(2, weight=1)
    root.grid_rowconfigure(3, weight=1)
    root.grid_rowconfigure(4, weight=1)
    root.grid_rowconfigure(5, weight=1)
    root.grid_rowconfigure(6, weight=1)
    root.grid_rowconfigure(7, weight=1)
    root.grid_rowconfigure(8, weight=1)
    root.grid_rowconfigure(9, weight=1)
    root.grid_rowconfigure(10, weight=1)
    root.grid_rowconfigure(11, weight=1)
    root.grid_rowconfigure(12, weight=1)
    root.grid_columnconfigure(0, weight=1)
    root.grid_columnconfigure(1, weight=1)
    root.grid_columnconfigure(2, weight=1)
    root.grid_columnconfigure(3, weight=1)

# Call the configure_resizing function after creating the widgets
configure_resizing()

# Add a callback function to handle window resizing events
def on_resize(event):
    configure_resizing()

root.bind("<Configure>", on_resize)

# ---------------------------------------------------------------------------------------------------

root.mainloop()

