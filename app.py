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
            intensity = sum(peaks[1]) if peaks[1].size > 0 else 0.0  
            times.append(time)  
            intensities.append(intensity)  
    return np.array(times), np.array(intensities)  
  
def extract_xic(experiment, mass, tol=0.5):  
    """  
    Extracts an Extracted Ion Chromatogram (XIC) from an MSExperiment.  
    Returns numpy arrays of retention times and intensities for the specified m/z.  
    """  
    times, intensities = [], []  
    for spectrum in experiment:  
        if spectrum.getMSLevel() == 1:  
            time = spectrum.getRT()  
            peaks = spectrum.get_peaks()  
            mz_array, intensity_array = peaks  
              
            # Find peaks within the tolerance range  
            matches = np.where((mz_array >= mass - tol) & (mz_array <= mass + tol))[0]  
            if matches.size > 0:  
                # Sum intensities of all matching peaks  
                intensity = np.sum(intensity_array[matches])  
            else:  
                intensity = 0.0  
                  
            times.append(time)  
            intensities.append(intensity)  
    return np.array(times), np.array(intensities)  
  
def extract_ms2_data(experiment):  
    """  
    Extracts MS2 spectra from an MSExperiment.  
    Returns a list of dictionaries containing MS2 data.  
    """  
    ms2_data = []  
    for i, spectrum in enumerate(experiment):  
        if spectrum.getMSLevel() == 2:  
            precursors = spectrum.getPrecursors()  
            if precursors:  
                precursor_mz = precursors[0].getMZ()  
                mz_array, intensity_array = spectrum.get_peaks()  
                  
                # Only add if we have valid data  
                if len(mz_array) > 0 and len(intensity_array) > 0:  
                    ms2_data.append({  
                        'index': i,  
                        'rt': spectrum.getRT(),  
                        'precursor': precursor_mz,  
                        'mz': mz_array,  
                        'intensity': intensity_array  
                    })  
    return ms2_data  
  
def parse_mgf(mgf_content):  
    """  
    Parses an MGF file content and returns a list of records.  
    Each record is a dictionary with mz_array, intensity_array, and other metadata.  
    """  
    records = []  
    current_record = None  
      
    for line in mgf_content.splitlines():  
        line = line.strip()  
        if not line:  
            continue  
              
        if line == "BEGIN IONS":  
            current_record = {  
                'mz_array': [],  
                'intensity_array': [],  
                'title': 'Unknown'  
            }  
        elif line == "END IONS":  
            if current_record:  
                # Convert lists to numpy arrays  
                current_record['mz_array'] = np.array(current_record['mz_array'])  
                current_record['intensity_array'] = np.array(current_record['intensity_array'])  
                records.append(current_record)  
                current_record = None  
        elif current_record:  
            if line.startswith("TITLE="):  
                current_record['title'] = line[6:]  
            elif line.startswith("PEPMASS="):  
                current_record['pepmass'] = float(line.split('=')[1].split()[0])  
            elif line.startswith("CHARGE="):  
                current_record['charge'] = line[7:]  
            elif "=" not in line:  # This is a peak  
                try:  
                    mz, intensity = map(float, line.split())  
                    current_record['mz_array'].append(mz)  
                    current_record['intensity_array'].append(intensity)  
                except ValueError:  
                    pass  # Skip invalid lines  
      
    return records  
  
def match_fragments(ms2_spec, mgf_record, tolerance=0.5):  
    """  
    Matches fragments between an MS2 spectrum and an MGF record.  
    Returns indices of matching fragments.  
    """  
    ms2_indices = []  
    mgf_indices = []  
      
    for i, ms2_mz in enumerate(ms2_spec['mz']):  
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
          
        # XIC Visualization  
        st.header("Extracted Ion Chromatogram (XIC)")  
        mass_value = st.number_input("Enter target m/z value:", value=500.0, step=1.0)  
        tol = st.number_input("Enter tolerance (Da):", value=0.5, step=0.1)  
        xic_times, xic_intensities = extract_xic(experiment, mass_value, tol)  
          
        fig_xic = go.Figure()  
        fig_xic.add_trace(go.Scatter(  
            x=xic_times,  
            y=xic_intensities,  
            mode='lines+markers',  
            line=dict(color="#24EB84", width=2),  
            name="XIC"  
        ))  
        fig_xic.update_layout(  
            title={"text": "Extracted Ion Chromatogram", "pad": {"t":15}},  
            xaxis_title="Retention Time",  
            yaxis_title="Intensity",  
            template="simple_white",  
            xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
            yaxis=dict(showgrid=True, gridcolor="#F3F4F6")  
        )  
        st.plotly_chart(fig_xic, use_container_width=True)  
          
        # MS2 Spectrum Display with Stem Plot  
        st.header("MS2 Fragmentation Analysis")  
        ms2_data = extract_ms2_data(experiment)  
        if ms2_data:  
            # Let user pick a spectrum based on retention time and precursor m/z  
            options = [f"RT: {spec['rt']:.2f}, Prec: {spec['precursor']:.4f}" for spec in ms2_data]  
            selection = st.selectbox("Select an MS2 spectrum to view:", options)  
            selected_index = options.index(selection)  
            selected_ms2 = ms2_data[selected_index]  
              
            # Validate and plot the selected MS2 spectrum  
            if len(selected_ms2['mz']) > 0 and len(selected_ms2['intensity']) > 0:  
                # Create stem plot data  
                x_stem, y_stem = create_stem_data(selected_ms2['mz'], selected_ms2['intensity'])  
                  
                # Plot MS2 spectrum as stem plot  
                fig_ms2 = go.Figure()  
                fig_ms2.add_trace(go.Scatter(  
                    x=x_stem,  
                    y=y_stem,  
                    mode='lines',  
                    line=dict(color="#B2EB24", width=1.5),  
                    name="MS2 Peaks"  
                ))  
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
                st.warning("The selected MS2 spectrum does not contain any peaks.")  
              
            # MGF Matching Section  
            st.header("MGF Matching")  
            uploaded_mgf = st.file_uploader("Upload an MGF file", type=["mgf"])  
            if uploaded_mgf is not None:  
                try:  
                    mgf_content = uploaded_mgf.read().decode('utf-8')  
                    # Parse MGF file  
                    records = parse_mgf(mgf_content)  
                    if records:  
                        st.success(f"MGF file loaded successfully! Found {len(records)} records.")  
                          
                        # Let user select a record to match against  
                        record_options = [f"Record {i+1}: {record.get('title', 'No Title')}" for i, record in enumerate(records)]  
                        record_selection = st.selectbox("Select an MGF record to match against:", record_options)  
                        record_index = record_options.index(record_selection)  
                        record = records[record_index]  
                          
                        if len(record.get("mz_array", [])) > 0:  
                            # Define tolerance in Da for matching  
                            tolerance = st.number_input("Set mass tolerance (Da):", value=0.5, step=0.1)  
                              
                            # Match fragments  
                            match_results = match_fragments(selected_ms2, record, tolerance)  
                              
                            if match_results["ms2_indices"]:  
                                st.success(f"Found {len(match_results['ms2_indices'])} matching fragments!")  
                                  
                                # Plot matches on the MS2 spectrum  
                                fig_match = go.Figure()  
                                  
                                # Plot all MS2 peaks as stem plot  
                                x_stem, y_stem = create_stem_data(selected_ms2['mz'], selected_ms2['intensity'])  
                                fig_match.add_trace(go.Scatter(  
                                    x=x_stem,  
                                    y=y_stem,  
                                    mode='lines',  
                                    line=dict(color="#B2EB24", width=1.5),  
                                    name="MS2 Peaks"  
                                ))  
                                  
                                # Highlight matched peaks with a different color  
                                matched_mz = [selected_ms2['mz'][i] for i in match_results['ms2_indices']]  
                                matched_intensity = [selected_ms2['intensity'][i] for i in match_results['ms2_indices']]  
                                  
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
                                    title={"text": "MS2 Spectrum with Matching Fragments", "pad": {"t":15}},  
                                    xaxis_title="m/z",  
                                    yaxis_title="Intensity",  
                                    template="simple_white",  
                                    xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                                    yaxis=dict(showgrid=True, gridcolor="#F3F4F6")  
                                )  
                                st.plotly_chart(fig_match, use_container_width=True)  
                                  
                                # Display matching details in a table  
                                match_data = []  
                                for ms2_idx, mgf_idx in zip(match_results['ms2_indices'], match_results['mgf_indices']):  
                                    ms2_mz = selected_ms2['mz'][ms2_idx]  
                                    ms2_intensity = selected_ms2['intensity'][ms2_idx]  
                                    mgf_mz = record['mz_array'][mgf_idx]  
                                    mgf_intensity = record['intensity_array'][mgf_idx]  
                                    delta = ms2_mz - mgf_mz  
                                      
                                    match_data.append({  
                                        "MS2 m/z": round(ms2_mz, 4),  
                                        "MS2 Intensity": round(ms2_intensity, 0),  
                                        "MGF m/z": round(mgf_mz, 4),  
                                        "MGF Intensity": round(mgf_intensity, 0),  
                                        "Delta (Da)": round(delta, 4)  
                                    })  
                                  
                                match_df = pd.DataFrame(match_data)  
                                st.dataframe(match_df)  
                            else:  
                                st.warning("No matching fragments found with the current tolerance.")  
                        else:  
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
