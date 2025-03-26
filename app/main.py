# app.py  
import streamlit as st  
import plotly.graph_objects as go  
import numpy as np  
import pandas as pd  
import tempfile  
import os  
import io  
import re  
import subprocess  
import sys  
  
# Import pyopenms and install if missing  
try:  
    from pyopenms import MSExperiment, MzMLFile  
except ImportError:  
    st.info("Installing pyopenms...")  
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyopenms"])  
    from pyopenms import MSExperiment, MzMLFile  
  
# -------------------------------------------------------------------  
# Streamlit page configuration  
st.set_page_config(page_title="PyOpenMS Mass Spectrometry Analyzer", layout="wide")  
st.title("Mass Spectrometry Data Analysis")  
  
st.markdown(  
    "Upload an mzML file to analyze your mass spectrometry data: view the Total Ion Chromatogram (TIC), "  
    "extract an XIC, explore MS2 spectra, and match fragmentation against a user-uploaded MGF file."  
)  
  
# -------------------------------------------------------------------  
# Helper Functions  
  
def load_mzml(uploaded_file):  
    """  
    Writes the uploaded mzML file to a temporary file and loads it   
    into an MSExperiment using pyopenms MzMLFile.load.  
    """  
    experiment = MSExperiment()  
    tmp_path = None  
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
        if tmp_path and os.path.exists(tmp_path):  
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
            intensity = sum(peaks[1]) if len(peaks[1]) > 0 else 0  
            times.append(time)  
            intensities.append(intensity)  
    return np.array(times), np.array(intensities)  
  
def extract_xic(experiment, target_mz, tolerance):  
    """  
    Extracts an Extracted Ion Chromatogram (XIC) for a specific m/z value.  
    Returns numpy arrays of retention times and intensities for the target m/z.  
    """  
    times, intensities = [], []  
    for spectrum in experiment:  
        if spectrum.getMSLevel() == 1:  
            time = spectrum.getRT()  
            mz_array, intensity_array = spectrum.get_peaks()  
              
            # Find peaks within the tolerance range  
            extracted_intensity = 0  
            for i, mz in enumerate(mz_array):  
                if abs(mz - target_mz) <= tolerance:  
                    extracted_intensity += intensity_array[i]  
              
            times.append(time)  
            intensities.append(extracted_intensity)  
      
    return np.array(times), np.array(intensities)  
  
def extract_ms2_data(experiment):  
    """  
    Extracts MS2 spectra data from an MSExperiment.  
    Returns a list of dictionaries containing MS2 spectra information.  
    """  
    ms2_data = []  
    for spectrum in experiment:  
        if spectrum.getMSLevel() == 2:  
            precursors = spectrum.getPrecursors()  
            precursor_mz = precursors[0].getMZ() if precursors else None  
              
            mz_array, intensity_array = spectrum.get_peaks()  
              
            ms2_data.append({  
                'rt': spectrum.getRT(),  
                'precursor': precursor_mz,  
                'mz': mz_array,  
                'intensity': intensity_array  
            })  
      
    return ms2_data  
  
def find_matches(ms2_mz_array, ms2_intensity_array, mgf_record, tolerance=0.5):  
    """  
    Finds matches between MS2 spectrum and MGF record.  
    Returns indices of matched peaks.  
    """  
    ms2_indices = []  
    mgf_indices = []  
      
    for i, ms2_mz in enumerate(ms2_mz_array):  
        for j, mgf_mz in enumerate(mgf_record['mz_array']):  
            if abs(ms2_mz - mgf_mz) <= tolerance:  
                ms2_indices.append(i)  
                mgf_indices.append(j)  
                break  
      
    return {  
        'ms2_indices': ms2_indices,  
        'mgf_indices': mgf_indices  
    }  
  
def create_stem_data(x, y):  
    """  
    Prepares data to simulate a stem (stick) plot.  
    For each point, create two points: one at y=0 and one at the actual intensity,  
    with None added to separate the lines.  
    """  
    x_stem = []  
    y_stem = []  
    for xi, yi in zip(x, y):  
        x_stem += [xi, xi, None]  
        y_stem += [0, yi, None]  
    return x_stem, y_stem  
  
# -------------------------------------------------------------------  
# Main App  
  
# File uploader for mzML files  
uploaded_mzml = st.file_uploader("Upload an mzML file", type=["mzML"])  
  
if uploaded_mzml is not None:  
    st.info("Loading mzML file...")  
    experiment = load_mzml(uploaded_mzml)  
      
    if experiment is not None:  
        st.success("mzML file loaded successfully!")  
          
        # Display basic information  
        st.header("File Information")  
        num_spectra = experiment.size()  
        ms1_count = sum(1 for spec in experiment if spec.getMSLevel() == 1)  
        ms2_count = sum(1 for spec in experiment if spec.getMSLevel() == 2)  
        col1, col2, col3 = st.columns(3)  
        col1.metric("Total Spectra", num_spectra)  
        col2.metric("MS1 Spectra", ms1_count)  
        col3.metric("MS2 Spectra", ms2_count)  
          
        # TIC Visualization  
        st.header("Total Ion Chromatogram (TIC)")  
        tic_times, tic_intensities = extract_tic(experiment)  
        fig_tic = go.Figure()  
        fig_tic.add_trace(go.Scatter(  
            x=tic_times,  
            y=tic_intensities,  
            mode='lines',  
            line=dict(color="#2563EB", width=2),  
            name="TIC"  
        ))  
        fig_tic.update_layout(  
            title={"text": "Total Ion Chromatogram", "pad": {"t":15}},  
            xaxis_title="Retention Time",  
            yaxis_title="Intensity",  
            template="simple_white",  
            xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
            yaxis=dict(showgrid=True, gridcolor="#F3F4F6")  
        )  
        st.plotly_chart(fig_tic, use_container_width=True)  
          
        # XIC Visualization with option to show/hide points  
        st.header("Extracted Ion Chromatogram (XIC)")  
        col1, col2, col3 = st.columns([2, 2, 1])  
        with col1:  
            mass_value = st.number_input("Enter target m/z value:", value=500.0, step=1.0)  
        with col2:  
            tol = st.number_input("Enter tolerance (Da):", value=0.5, step=0.1)  
        with col3:  
            show_points = st.checkbox("Show points", value=True)  
          
        xic_times, xic_intensities = extract_xic(experiment, mass_value, tol)  
        fig_xic = go.Figure()  
          
        # Determine plot mode based on user preference  
        plot_mode = 'lines+markers' if show_points else 'lines'  
          
        fig_xic.add_trace(go.Scatter(  
            x=xic_times,  
            y=xic_intensities,  
            mode=plot_mode,  
            line=dict(color="#24EB84", width=2),  
            marker=dict(size=6),  
            name="XIC"  
        ))  
        fig_xic.update_layout(  
            title={"text": f"Extracted Ion Chromatogram (m/z: {mass_value} ± {tol})", "pad": {"t":15}},  
            xaxis_title="Retention Time",  
            yaxis_title="Intensity",  
            template="simple_white",  
            xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
            yaxis=dict(showgrid=True, gridcolor="#F3F4F6")  
        )  
        st.plotly_chart(fig_xic, use_container_width=True)  
          
        # MS2 Spectrum Display as stem plot  
        st.header("MS2 Fragmentation Analysis")  
        ms2_data = extract_ms2_data(experiment)  
        if ms2_data:  
            # Let user pick a spectrum based on retention time and precursor m/z  
            options = [f"RT: {spec['rt']:.2f}, Prec: {spec['precursor']:.4f}" if spec['precursor'] is not None else f"RT: {spec['rt']:.2f}, Prec: N/A"   
                       for spec in ms2_data]  
            selection = st.selectbox("Select an MS2 spectrum to view:", options)  
            selected_index = options.index(selection)  
            selected_ms2 = ms2_data[selected_index]  
              
            # Validate and plot the selected MS2 spectrum as a stem plot  
            if len(selected_ms2['mz']) > 0 and len(selected_ms2['intensity']) > 0:  
                fig_ms2 = go.Figure()  
                  
                # Create stem plot data  
                x_stem, y_stem = create_stem_data(selected_ms2['mz'], selected_ms2['intensity'])  
                  
                # Plot MS2 spectrum as stem plot  
                fig_ms2.add_trace(go.Scatter(  
                    x=x_stem,  
                    y=y_stem,  
                    mode='lines',  
                    line=dict(color="#B2EB24", width=1.5),  
                    name="MS2 Peaks"  
                ))  
                  
                fig_ms2.update_layout(  
                    title={"text": "MS2 Mass Spectrum", "pad": {"t":15}},  
                    xaxis_title="m/z",  
                    yaxis_title="Intensity",  
                    template="simple_white",  
                    xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                    yaxis=dict(showgrid=True, gridcolor="#F3F4F6")  
                )  
                st.plotly_chart(fig_ms2, use_container_width=True)  
            else:  
                st.warning("The selected MS2 spectrum does not contain any peaks.")  
                  
            # Optional: MGF Matching section (if provided by user)  
            st.subheader("MGF Matching")  
            mgf_file = st.file_uploader("Upload an MGF file for matching (optional)", type=["mgf"])  
            if mgf_file is not None:  
                try:  
                    mgf_content = mgf_file.read().decode("utf-8").splitlines()  
                    # Parse MGF file (basic parsing)  
                    mgf_records = []  
                    current = {}  
                    for line in mgf_content:  
                        line = line.strip()  
                        if line.startswith("BEGIN IONS"):  
                            current = {"mz_array": [], "intensity_array": []}  
                        elif line.startswith("END IONS"):  
                            if current and len(current["mz_array"]) > 0:  
                                mgf_records.append(current)  
                            current = {}  
                        elif "=" in line:  
                            key, value = line.split("=", 1)  
                            current[key] = value  
                        elif line and not line.startswith("#"):  
                            # Parse peak data (m/z intensity pairs)  
                            parts = line.split()  
                            if len(parts) >= 2:  
                                try:  
                                    mz = float(parts[0])  
                                    intensity = float(parts[1])  
                                    current["mz_array"].append(mz)  
                                    current["intensity_array"].append(intensity)  
                                except ValueError:  
                                    pass  
                      
                    if mgf_records:  
                        # Let user select an MGF record to match  
                        mgf_options = [f"Record {i+1}" for i in range(len(mgf_records))]  
                        mgf_selection = st.selectbox("Select an MGF record to match:", mgf_options)  
                        mgf_index = mgf_options.index(mgf_selection)  
                        selected_mgf = mgf_records[mgf_index]  
                          
                        # Set tolerance for matching  
                        match_tolerance = st.number_input("Matching tolerance (Da):", value=0.5, step=0.1)  
                          
                        # Find matches  
                        matches = find_matches(  
                            selected_ms2['mz'],   
                            selected_ms2['intensity'],   
                            selected_mgf,   
                            tolerance=match_tolerance  
                        )  
                          
                        # Extract matched peaks  
                        matched_mz = [selected_ms2['mz'][i] for i in matches['ms2_indices']]  
                        matched_intensity = [selected_ms2['intensity'][i] for i in matches['ms2_indices']]  
                          
                        # Display matching results  
                        st.write(f"Found {len(matched_mz)} matching peaks.")  
                          
                        # Plot matched spectrum  
                        if matched_mz:  
                            fig_match = go.Figure()  
                              
                            # Create stem plot for all MS2 peaks  
                            x_stem, y_stem = create_stem_data(selected_ms2['mz'], selected_ms2['intensity'])  
                            fig_match.add_trace(go.Scatter(  
                                x=x_stem,  
                                y=y_stem,  
                                mode='lines',  
                                line=dict(color="#B2EB24", width=1.5),  
                                name="MS2 Peaks"  
                            ))  
                              
                            # Create stem plot for matched peaks  
                            x_matched_stem, y_matched_stem = create_stem_data(matched_mz, matched_intensity)  
                            fig_match.add_trace(go.Scatter(  
                                x=x_matched_stem,  
                                y=y_matched_stem,  
                                mode='lines',  
                                line=dict(color="#EB3424", width=2.5),  
                                name="Matched Fragments"  
                            ))  
                              
                            fig_match.update_layout(  
                                title={"text": "MS2 Spectrum with Matched Fragments", "pad": {"t":15}},  
                                xaxis_title="m/z",  
                                yaxis_title="Intensity",  
                                template="simple_white",  
                                xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                                yaxis=dict(showgrid=True, gridcolor="#F3F4F6")  
                            )  
                            st.plotly_chart(fig_match, use_container_width=True)  
                              
                            # Display matched peaks in a table  
                            match_data = []  
                            for i, idx in enumerate(matches['ms2_indices']):  
                                ms2_mz = selected_ms2['mz'][idx]  
                                ms2_intensity = selected_ms2['intensity'][idx]  
                                mgf_mz = selected_mgf['mz_array'][matches['mgf_indices'][i]]  
                                match_data.append({  
                                    "MS2 m/z": round(ms2_mz, 4),  
                                    "MS2 Intensity": int(ms2_intensity),  
                                    "MGF m/z": round(mgf_mz, 4),  
                                    "Difference (Da)": round(ms2_mz - mgf_mz, 4)  
                                })  
                              
                            if match_data:  
                                st.dataframe(pd.DataFrame(match_data))  
                        else:  
                            st.warning("No matches found with the current tolerance.")  
                    elif len(mgf_records) == 0:  
                        st.error("No valid m/z values found in the MGF file.")  
                    else:  
                        st.error("No valid records found in the MGF file.")  
                except Exception as e:  
                    st.error(f"Error processing MGF file: {str(e)}")  
        else:  
            st.warning("No MS2 spectra found in the file.")  
    else:  
        st.error("Failed to load the mzML file. Please check the file format.")  
else:  
    # Show features when no file is uploaded  
    st.info("Please upload an mzML file to begin analysis.")  
      
    st.subheader("Features")  
    st.markdown("""  
    - **Total Ion Chromatogram (TIC)**: View the total ion current across retention time  
    - **Extracted Ion Chromatogram (XIC)**: Extract and visualize specific m/z values  
    - **MS2 Fragmentation Analysis**: View MS2 spectra and analyze fragmentation patterns  
    - **MGF Matching**: Match experimental MS2 data with theoretical fragments from MGF files  
    """)  
