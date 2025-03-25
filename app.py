# app.py  
import streamlit as st  
import plotly.graph_objects as go  
import numpy as np  
import pandas as pd  
import tempfile  
import os  
import io  
import re  
  
# Import pyopenms and install if missing  
try:  
    from pyopenms import MSExperiment, MzMLFile  
except ImportError:  
    import sys  
    !{sys.executable} -m pip install pyopenms  
    from pyopenms import MSExperiment, MzMLFile  
  
# -------------------------------------------------------------------  
# Streamlit page configuration  
st.set_page_config(page_title="PyOpenMS Mass Spectrometry Analyzer", layout="wide")  
st.title("Mass Spectrometry Data Analysis")  
st.markdown(  
    "Upload an mzML file to analyze your mass spectrometry data: view the Total Ion Chromatogram (TIC), extract an XIC, "  
    "explore MS2 spectra, and match fragmentation against a user-uploaded MGF file."  
)  
  
# -------------------------------------------------------------------  
# Helper Functions  
  
def load_mzml(uploaded_file):  
    """  
    Writes the uploaded mzML file to a temporary file and loads it  
    into an MSExperiment using pyopenms MzMLFile.load.  
    """  
    experiment = MSExperiment()  
    try:  
        # Write file to temporary location  
        with tempfile.NamedTemporaryFile(delete=False, suffix=".mzML") as tmp:  
            tmp.write(uploaded_file.read())  
            tmp.flush()  
            tmp_path = tmp.name  
        mzml_file = MzMLFile()  
        mzml_file.load(tmp_path, experiment)  
    except Exception as e:  
        st.error("Error loading mzML file: " + str(e))  
        experiment = None  
    finally:  
        if os.path.exists(tmp_path):  
            os.remove(tmp_path)  
    return experiment  
  
def extract_tic(experiment):  
    """  
    Extracts the Total Ion Chromatogram (TIC) from an MSExperiment.  
    Returns numpy arrays of retention times and summed intensities.  
    """  
    times, intensities = [], []  
    for spectrum in experiment:  
        if spectrum.getMSLevel() == 1:  
            time = spectrum.getRT()  
            peaks = spectrum.get_peaks()  
            intensity = sum(peaks[1]) if peaks[1].size > 0 else 0.0  
            times.append(time)  
            intensities.append(intensity)  
    return np.array(times), np.array(intensities)  
  
def extract_xic(experiment, mass, tol=0.5):  
    """  
    Extracts the Extracted Ion Chromatogram (XIC) for the given mass (m/z) and tolerance.  
    Returns numpy arrays of retention times and summed intensities within the m/z window.  
    """  
    times, intensities = [], []  
    for spectrum in experiment:  
        if spectrum.getMSLevel() == 1:  
            time = spectrum.getRT()  
            mz_array, intensity_array = spectrum.get_peaks()  
            mask = (mz_array >= mass - tol) & (mz_array <= mass + tol)  
            intensity = np.sum(intensity_array[mask]) if np.any(mask) else 0.0  
            times.append(time)  
            intensities.append(intensity)  
    return np.array(times), np.array(intensities)  
  
def extract_ms2_data(experiment):  
    """  
    Extracts MS2 spectra from the experiment.  
    Returns a list of dictionaries containing precursor, retention time, m/z array, and intensity array.  
    """  
    ms2_data = []  
    for spectrum in experiment:  
        if spectrum.getMSLevel() == 2:  
            precursors = spectrum.getPrecursors()  
            if precursors:  
                precursor_mz = precursors[0].getMZ()  
                mz_array, intensity_array = spectrum.get_peaks()  
                ms2_data.append({  
                    'rt': spectrum.getRT(),  
                    'precursor': precursor_mz,  
                    'mz': mz_array,  
                    'intensity': intensity_array  
                })  
    return ms2_data  
  
def parse_mgf(mgf_bytes):  
    """  
    Parses an MGF file (bytes) and returns a list of records.  
    Each record is a dictionary with keys like 'pepmass', 'charge', 'mz_array', and 'intensity_array'.  
    This parser assumes a simplified MGF format.  
    """  
    records = []  
    text = mgf_bytes.decode('utf-8')  
    blocks = re.split(r'BEGIN IONS', text, flags=re.IGNORECASE)  
    for block in blocks:  
        if "END IONS" in block:  
            record = {}  
            lines = block.strip().splitlines()  
            mzs = []  
            intensities = []  
            for line in lines:  
                if "=" in line:  
                    key, value = line.split("=", 1)  
                    key = key.strip().lower()  
                    value = value.strip()  
                    # We assume pepmass and charge are provided like 'pepmass=1234.56' and 'charge=2+'  
                    if key == "pepmass":  
                        try:  
                            record["pepmass"] = float(value.split()[0])  
                        except:  
                            record["pepmass"] = None  
                    elif key == "charge":  
                        record["charge"] = value  
                else:  
                    parts = line.split()  
                    if len(parts) == 2:  
                        try:  
                            mz_val = float(parts[0])  
                            inten_val = float(parts[1])  
                            mzs.append(mz_val)  
                            intensities.append(inten_val)  
                        except:  
                            continue  
            record["mz_array"] = mzs  
            record["intensity_array"] = intensities  
            records.append(record)  
    return records  
  
# -------------------------------------------------------------------  
# Main Application Logic  
  
uploaded_mzml = st.file_uploader("Upload an mzML file", type=["mzML"])  
if uploaded_mzml is not None:  
    st.info("Loading mzML file...")  
    experiment = load_mzml(uploaded_mzml)  
    if experiment is None:  
        st.error("Failed to load the mzML file.")  
    else:  
        st.success("mzML file loaded successfully!")  
          
        # Total Ion Chromatogram (TIC)  
        st.header("Total Ion Chromatogram (TIC)")  
        tic_times, tic_intensities = extract_tic(experiment)  
        fig_tic = go.Figure()  
        fig_tic.add_trace(go.Scatter(x=tic_times, y=tic_intensities, mode="lines", line=dict(color="#2563EB")))  
        fig_tic.update_layout(  
            title={"text": "TIC", "pad": {"t":15}},  
            xaxis_title="Retention Time",  
            yaxis_title="Intensity",  
            template="simple_white",  
            xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
            yaxis=dict(showgrid=True, gridcolor="#F3F4F6")  
        )  
        st.plotly_chart(fig_tic, use_container_width=True)  
          
        # Extracted Ion Chromatogram (XIC)  
        st.header("Extracted Ion Chromatogram (XIC)")  
        target_mass = st.number_input("Enter target m/z value", value=500.0, step=0.1)  
        tolerance = st.number_input("Enter tolerance (Da)", value=0.5, step=0.1)  
        xic_times, xic_intensities = extract_xic(experiment, target_mass, tol=tolerance)  
        fig_xic = go.Figure()  
        fig_xic.add_trace(go.Scatter(x=xic_times, y=xic_intensities, mode="lines", line=dict(color="#24EB84")))  
        fig_xic.update_layout(  
            title={"text": "XIC", "pad": {"t":15}},  
            xaxis_title="Retention Time",  
            yaxis_title="Intensity",  
            template="simple_white",  
            xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
            yaxis=dict(showgrid=True, gridcolor="#F3F4F6")  
        )  
        st.plotly_chart(fig_xic, use_container_width=True)  
          
        # MS2 Spectra Visualization  
        st.header("MS2 Spectra")  
        ms2_data = extract_ms2_data(experiment)  
        if ms2_data:  
            ms2_idx = st.selectbox("Select an MS2 spectrum", list(range(len(ms2_data))))  
            selected_ms2 = ms2_data[ms2_idx]  
            st.write("Selected MS2 Spectrum - Precursor m/z: " + str(selected_ms2['precursor']) + ", RT: " + str(selected_ms2['rt']))  
            fig_ms2 = go.Figure()  
            fig_ms2.add_trace(go.Bar(x=selected_ms2['mz'], y=selected_ms2['intensity'], marker_color="#B2EB24"))  
            fig_ms2.update_layout(  
                title={"text": "MS2 Fragmentation Spectrum", "pad": {"t":15}},  
                xaxis_title="m/z",  
                yaxis_title="Intensity",  
                template="simple_white",  
                xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                yaxis=dict(showgrid=True, gridcolor="#F3F4F6")  
            )  
            st.plotly_chart(fig_ms2, use_container_width=True)  
        else:  
            st.info("No MS2 spectra found in the uploaded mzML file.")  
          
        # MGF Fragmentation Matching Section  
        st.header("MGF Fragmentation Matching")  
        uploaded_mgf = st.file_uploader("Upload an MGF file", type=["mgf"], key="mgfUploader")  
        if uploaded_mgf is not None:  
            try:  
                st.info("Parsing MGF file...")  
                mgf_bytes = uploaded_mgf.getvalue()  
                mgf_records = parse_mgf(mgf_bytes)  
                if mgf_records:  
                    st.success("MGF file parsed successfully!")  
                      
                    mgf_idx = st.selectbox("Select an MGF record", list(range(len(mgf_records))), key="mgfSelect")  
                    record = mgf_records[mgf_idx]  
                    st.write("Peptide Mass (precursor): " + str(record.get("pepmass", "N/A")))  
                    st.write("Charge: " + str(record.get("charge", "N/A")))  
                      
                    fig_mgf = go.Figure()  
                    fig_mgf.add_trace(go.Bar(  
                        x=record["mz_array"],  
                        y=record["intensity_array"],  
                       
