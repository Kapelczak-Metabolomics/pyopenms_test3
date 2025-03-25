# app.py  
import streamlit as st  
import plotly.graph_objects as go  
from pyopenms import MSExperiment, MzMLFile, MSSpectrum, SpectrumLookup, TraMLFile, FeatureXMLFile, FeatureMap, PeakMap  
import numpy as np  
import io  
  
# --- Helper Functions ---  
def load_mzml(file_bytes):  
    """  
    Load an mzML file from uploaded bytes and return an MSExperiment.  
    """  
    mzml_data = MSExperiment()  
    mzml_file = MzMLFile()  
    mzml_file.load(io.BytesIO(file_bytes), mzml_data)  
    return mzml_data  
  
def extract_tic(experiment):  
    """  
    Extract the Total Ion Chromatogram (TIC) from an MSExperiment.  
    Returns times and intensities.  
    """  
    times = []  
    intensities = []  
    # Iterate over individual spectra  
    for spectrum in experiment:  
        if spectrum.getMSLevel() == 1:  
            time = spectrum.getRT()  # retention time  
            intens = sum(spectrum.get_peaks()[1]) if len(spectrum.get_peaks()) > 1 else 0.0  
            times.append(time)  
            intensities.append(intens)  
    return np.array(times), np.array(intensities)  
  
def extract_xic(experiment, mass, tol=0.5):  
    """  
    Extract an extracted ion chromatogram (XIC) for a given mass.  
    The tolerance is in Da.  
    Returns times and intensities.  
    """  
    times = []  
    intensities = []  
    for spectrum in experiment:  
        # Only consider MS1 spectra for XIC extraction.  
        if spectrum.getMSLevel() == 1:  
            time = spectrum.getRT()  
            mz_array, intensity_array = spectrum.get_peaks()  
            mz_array = np.array(mz_array)  
            intensity_array = np.array(intensity_array)  
            # Select intensities within the tolerance window  
            selection = (mz_array >= mass - tol) & (mz_array <= mass + tol)  
            intens = intensity_array[selection].sum() if selection.any() else 0.0  
            times.append(time)  
            intensities.append(intens)  
    return np.array(times), np.array(intensities)  
  
def extract_ms2(experiment):  
    """  
    Extract MS2 spectra from the MSExperiment.  
    Returns a list of tuples (precursor_mass, fragment m/z, intensity) for each MS2 spectrum.  
    """  
    ms2_data = []  
    for spectrum in experiment:  
        if spectrum.getMSLevel() == 2:  
            # Get precursor information from the spectrum's precursors list  
            precursors = spectrum.getPrecursors()  
            if precursors:  
                precursor_mz = precursors[0].getMZ()  
            else:  
                precursor_mz = None  
            # Get fragment ion peaks  
            mz_array, intensity_array = spectrum.get_peaks()  
            ms2_data.append({  
                'precursor': precursor_mz,  
                'mz': mz_array,  
                'intensity': intensity_array  
            })  
    return ms2_data  
  
def parse_mgf(file_bytes):  
    """  
    Parse an MGF file and return a list of dictionaries with keys:  
    'title', 'pepmass', 'charge', 'mz_array', 'intensity_array'.  
    """  
    mgf_data = []  
    content = io.StringIO(file_bytes.decode("utf-8"))  
    current_record = {}  
    mzs = []  
    intensities = []  
    for line in content:  
        line = line.strip()  
        if not line:   
            continue  
        if line.startswith("BEGIN IONS"):  
            current_record = {}  
            mzs = []  
            intensities = []  
        elif line.startswith("END IONS"):  
            current_record['mz_array'] = np.array(mzs)  
            current_record['intensity_array'] = np.array(intensities)  
            mgf_data.append(current_record)  
        elif "=" in line:  
            key, value = line.split("=")  
            key = key.strip().lower()  
            value = value.strip()  
            # Convert appropriate types  
            if key == "pepmass":  
                parts = value.split()  
                current_record['pepmass'] = float(parts[0])  
            elif key == "charge":  
                # remove any trailing symbols e.g., "2+"  
                current_record['charge'] = int(''.join(filter(str.isdigit, value)))  
            else:  
                current_record[key] = value  
        else:  
            # This should be m/z and intensity values  
            parts = line.split()  
            if len(parts) == 2:  
                mzs.append(float(parts[0]))  
                intensities.append(float(parts[1]))  
    return mgf_data  
  
# --- Streamlit App Layout ---  
st.title("pyOpenMS Viewer and MS Analysis App")  
st.markdown("Upload an mzML file to view chromatograms, extract XICs, and inspect MS2 data. Optionally, upload an MGF file for fragmentation matching.")  
  
# File uploader for mzML file  
uploaded_mzml = st.file_uploader("Choose an mzML file", type=["mzML"])  
if uploaded_mzml is not None:  
    st.info("Loading mzML file...")  
    experiment = load_mzml(uploaded_mzml.getvalue())  
    st.success("mzML file successfully loaded!")  
      
    # Display Total Ion Chromatogram (TIC)  
    tic_times, tic_intensities = extract_tic(experiment)  
    st.subheader("Total Ion Chromatogram (TIC)")  
    fig_tic = go.Figure()  
    fig_tic.add_trace(go.Scatter(x=tic_times, y=tic_intensities, mode='lines', line=dict(color="#2563EB")))  
    fig_tic.update_layout(  
        title={"text": "TIC", "pad": {"t":15}},  
        xaxis_title="Retention Time",  
        yaxis_title="Total Intensity",  
        xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
        yaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
        template="simple_white"  
    )  
    st.plotly_chart(fig_tic, use_container_width=True)  
      
    # XIC Extraction Interface  
    st.subheader("Extracted Ion Chromatogram (XIC)")  
    mass_val = st.number_input("Enter the mass to extract XIC", min_value=0.0, value=500.0, step=0.1)  
    tol_val = st.number_input("Tolerance (Da)", min_value=0.0, value=0.5, step=0.1)  
    if st.button("Extract XIC"):  
        xic_times, xic_intens = extract_xic(experiment, mass_val, tol_val)  
        fig_xic = go.Figure()  
        fig_xic.add_trace(go.Scatter(x=xic_times, y=xic_intens, mode='lines', line=dict(color="#24EB84")))  
        fig_xic.update_layout(  
            title={"text": "XIC for mass " + str(mass_val), "pad": {"t":15}},  
            xaxis_title="Retention Time",  
            yaxis_title="Intensity",  
            xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
            yaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
            template="simple_white"  
        )  
        st.plotly_chart(fig_xic, use_container_width=True)  
      
    # MS2 Data Display  
    st.subheader("MS2 Data")  
    ms2_data = extract_ms2(experiment)  
    if ms2_data:  
        ms2_idx = st.selectbox("Select an MS2 spectrum", list(range(len(ms2_data))))  
        selected_ms2 = ms2_data[ms2_idx]  
        st.write("Precursor m/z: " + str(selected_ms2['precursor']))  
        fig_ms2 = go.Figure()  
        fig_ms2.add_trace(go.Bar(x=selected_ms2['mz'], y=selected_ms2['intensity'], marker_color="#B2EB24"))  
        fig_ms2.update_layout(  
            title={"text": "MS2 Fragmentation Spectrum", "pad": {"t":15}},  
            xaxis_title="m/z",  
            yaxis_title="Intensity",  
            xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
            yaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
            template="simple_white"  
        )  
        st.plotly_chart(fig_ms2, use_container_width=True)  
    else:  
        st.info("No MS2 data found in the uploaded mzML file.")  
  
# File uploader for MGF file  
uploaded_mgf = st.file_uploader("Choose an MGF file for fragmentation matching", type=["mgf"])  
if uploaded_mgf is not None:  
    st.info("Parsing MGF file for fragmentation matching...")  
    mgf_records = parse_mgf(uploaded_mgf.getvalue())  
    st.success("MGF file parsed successfully!")  
      
    st.subheader("MGF Fragmentation Data")  
    mgf_idx = st.selectbox("Select an MGF record", list(range(len(mgf_records))))  
    record = mgf_records[mgf_idx]  
    st.write("Peptide Mass (precursor): " + str(record.get('pepmass', 'N/A')))  
    st.write("Charge: " + str(record.get('charge', 'N/A')))  
      
    # Display fragmentation pattern for the selected record  
    fig_mgf = go.Figure()  
    fig_mgf.add_trace(go.Bar(x=record['mz_array'], y=record['intensity_array'], marker_color="#EB3424"))  
    fig_mgf.update_layout(  
        title={"text": "MGF Fragmentation Pattern", "pad": {"t":15}},  
        xaxis_title="m/z",  
        yaxis_title="Intensity",  
        xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
        yaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
        template="simple_white"  
    )  
    st.plotly_chart(fig_mgf, use_container_width=True)  
      
    # Simple Fragmentation Matching  
    st.subheader("Fragmentation Matching")  
    # For demonstration, define a simple matching: compare MS2 fragments with MGF fragments.  
    if 'ms2_data' in locals() and ms2_data:  
        match_idx = st.selectbox("Select a MS2 spectrum for matching", list(range(len(ms2_data))), key="match")  
        ms2_frag = ms
