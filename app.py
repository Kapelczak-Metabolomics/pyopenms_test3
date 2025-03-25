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
    print("Installing pyopenms...")  
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyopenms"])  
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
    Extracts an Extracted Ion Chromatogram (XIC) for a given mass with tolerance.  
    Returns numpy arrays of retention times and summed intensities within the mass window.  
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
                ms2_data.append({  
                    'index': i,  
                    'rt': spectrum.getRT(),  
                    'precursor': precursor_mz,  
                    'mz': mz_array,  
                    'intensity': intensity_array  
                })  
    return ms2_data  
  
def parse_mgf(mgf_bytes):  
    """  
    Parses an MGF file from bytes.  
    Returns a list of dictionaries containing MGF records.  
    """  
    mgf_records = []  
    try:  
        mgf_text = mgf_bytes.decode('utf-8')  
        # Split by BEGIN IONS/END IONS blocks  
        blocks = re.split(r'BEGIN IONS|END IONS', mgf_text)  
        for i in range(1, len(blocks), 2):  
            if i < len(blocks):  
                block = blocks[i].strip()  
                record = {}  
                mz_values = []  
                intensity_values = []  
                  
                lines = block.split('\n')  
                for line in lines:  
                    line = line.strip()  
                    if not line:  
                        continue  
                      
                    # Parse metadata lines (e.g., PEPMASS, CHARGE)  
                    if '=' in line:  
                        key, value = line.split('=', 1)  
                        record[key.lower()] = value  
                    # Parse peak lines (m/z intensity pairs)  
                    elif ' ' in line:  
                        parts = line.split()  
                        if len(parts) >= 2:  
                            try:  
                                mz = float(parts[0])  
                                intensity = float(parts[1])  
                                mz_values.append(mz)  
                                intensity_values.append(intensity)  
                            except ValueError:  
                                pass  
                  
                if mz_values and intensity_values:  
                    record['mz_array'] = np.array(mz_values)  
                    record['intensity_array'] = np.array(intensity_values)  
                    mgf_records.append(record)  
    except Exception as e:  
        st.error(f"Error parsing MGF file: {str(e)}")  
      
    return mgf_records  
  
def match_fragments(ms2_spec, mgf_record, tolerance=0.5):  
    """  
    Matches fragments between an MS2 spectrum and an MGF record.  
    Returns a dictionary with matching indices.  
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
          
        # XIC Extraction  
        st.header("Extracted Ion Chromatogram (XIC)")  
        col1, col2 = st.columns(2)  
        with col1:  
            xic_mass = st.number_input("m/z value", min_value=50.0, max_value=2000.0, value=500.0, step=0.1)  
        with col2:  
            xic_tol = st.number_input("Tolerance (Da)", min_value=0.01, max_value=10.0, value=0.5, step=0.1)  
          
        xic_times, xic_intensities = extract_xic(experiment, xic_mass, xic_tol)  
          
        fig_xic = go.Figure()  
        fig_xic.add_trace(go.Scatter(  
            x=xic_times,  
            y=xic_intensities,  
            mode='lines',  
            line=dict(color="#24EB84", width=2),  
            name=f"XIC (m/z {xic_mass}±{xic_tol})"  
        ))  
        fig_xic.update_layout(  
            title={"text": f"Extracted Ion Chromatogram (m/z {xic_mass}±{xic_tol})", "pad": {"t":15}},  
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
            st.write(f"Selected MS2 Spectrum - Precursor m/z: {selected_ms2['precursor']}, RT: {selected_ms2['rt']}")  
              
            fig_ms2 = go.Figure()  
            fig_ms2.add_trace(go.Bar(  
                x=selected_ms2['mz'],  
                y=selected_ms2['intensity'],  
                marker_color="#B2EB24"  
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
                    st.success(f"MGF file parsed successfully! Found {len(mgf_records)} records.")  
                      
                    mgf_idx = st.selectbox("Select an MGF record", list(range(len(mgf_records))), key="mgfSelect")  
                    record = mgf_records[mgf_idx]  
                    st.write(f"Peptide Mass (precursor): {record.get('pepmass', 'N/A')}")  
                    st.write(f"Charge: {record.get('charge', 'N/A')}")  
                      
                    fig_mgf = go.Figure()  
                    fig_mgf.add_trace(go.Bar(  
                        x=record["mz_array"],  
                        y=record["intensity_array"],  
                        marker_color="#EB3424"  
                    ))  
                    fig_mgf.update_layout(  
                        title={"text": "MGF Fragmentation Pattern", "pad": {"t":15}},  
                        xaxis_title="m/z",  
                        yaxis_title="Intensity",  
                        template="simple_white",  
                        xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                        yaxis=dict(showgrid=True, gridcolor="#F3F4F6")  
                    )  
                    st.plotly_chart(fig_mgf, use_container_width=True)  
                      
                    # Fragmentation matching  
                    if ms2_data:  
                        st.subheader("Fragment Matching")  
                        match_tol = st.number_input("Matching tolerance (Da)", min_value=0.01, max_value=2.0, value=0.5, step=0.1, key="matchTol")  
                          
                        ms2_match_idx = st.selectbox("Select MS2 spectrum for matching", list(range(len(ms2_data))), key="ms2MatchSelect")  
                        ms2_spec = ms2_data[ms2_match_idx]  
                          
                        match_results = match_fragments(ms2_spec, record, match_tol)  
                          
                        if match_results['ms2_indices']:  
                            st.success(f"Found {len(match_results['ms2_indices'])} matching fragments!")  
                              
                            # Create a visualization of the matches  
                            fig_match = go.Figure()  
                              
                            # MS2 spectrum (positive y-axis)  
                            fig_match.add_trace(go.Bar(  
                                x=ms2_spec['mz'],  
                                y=ms2_spec['intensity'],  
                                name="MS2 Spectrum",  
                                marker_color="#B2EB24"  
                            ))  
                              
                            # MGF spectrum (negative y-axis for comparison)  
                            fig_match.add_trace(go.Bar(  
                                x=record['mz_array'],  
                                y=-record['intensity_array'],  # Negative for display below x-axis  
                                name="MGF Spectrum",  
                                marker_color="#EB3424"  
                            ))  
                              
                            # Add lines connecting matching fragments  
                            for i, j in zip(match_results['ms2_indices'], match_results['mgf_indices']):  
                                ms2_mz = ms2_spec['mz'][i]  
                                ms2_intensity = ms2_spec['intensity'][i]  
                                mgf_mz = record['mz_array'][j]  
                                mgf_intensity = -record['intensity_array'][j]  # Negative for display  
                                  
                                fig_match.add_shape(  
                                    type="line",  
                                    x0=ms2_mz, y0=ms2_intensity,  
                                    x1=mgf_mz, y1=mgf_intensity,  
                                    line=dict(color="#2563EB", width=1, dash="dot")  
                                )  
                              
                            fig_match.update_layout(  
                                title={"text": "Fragment Matching Visualization", "pad": {"t":15}},  
                                xaxis_title="m/z",  
                                yaxis_title="Intensity",  
                                template="simple_white",  
                                xaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                                yaxis=dict(showgrid=True, gridcolor="#F3F4F6"),  
                                legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)  
                            )  
                            st.plotly_chart(fig_match, use_container_width=True)  
                              
                            # Show matching table  
                            match_data = {  
                                "MS2 m/z": [ms2_spec['mz'][i] for i in match_results['ms2_indices']],  
                                "MGF m/z": [record['mz_array'][j] for j in match_results['mgf_indices']],  
                                "Difference (Da)": [ms2_spec['mz'][i] - record['mz_array'][j] for i, j in zip(match_results['ms2_indices'], match_results['mgf_indices'])]  
                            }  
                            match_df = pd.DataFrame(match_data)  
                            st.dataframe(match_df)  
                        else:  
                            st.warning("No matching fragments found with the current tolerance.")  
                else:  
                    st.error("No valid records found in the MGF file.")  
            except Exception as e:  
                st.error(f"Error processing MGF file: {str(e)}")  
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
