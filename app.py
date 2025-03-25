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
            mz_array, intensity_array = spectrum.get_peaks()  
              
            # Find peaks within the tolerance range  
            matches = np.where(np.abs(mz_array - mass) <= tol)[0]  
            intensity = sum(intensity_array[matches]) if matches.size > 0 else 0.0  
              
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
    Each record contains the peptide sequence, charge, and m/z array.  
    """  
    records = []  
    current_record = None  
      
    for line in mgf_content.split('\n'):  
        line = line.strip()  
        if not line:  
            continue  
              
        if line == "BEGIN IONS":  
            current_record = {  
                'title': '',  
                'peptide': '',  
                'charge': 0,  
                'mz_array': [],  
                'intensity_array': []  
            }  
        elif line == "END IONS":  
            if current_record:  
                records.append(current_record)  
                current_record = None  
        elif current_record:  
            if line.startswith("TITLE="):  
                current_record['title'] = line[6:]  
            elif line.startswith("PEPMASS="):  
                parts = line[8:].split()  
                if parts:  
                    current_record['pepmass'] = float(parts[0])  
            elif line.startswith("CHARGE="):  
                charge_str = line[7:].replace('+', '')  
                current_record['charge'] = int(charge_str)  
            elif line.startswith("SEQ="):  
                current_record['peptide'] = line[4:]  
            elif re.match(r'^\d', line):  
                parts = line.split()  
                if len(parts) >= 2:  
                    mz = float(parts[0])  
                    intensity = float(parts[1])  
                    current_record['mz_array'].append(mz)  
                    current_record['intensity_array'].append(intensity)  
      
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
          
        # MS2 Spectrum Display  
        st.header("MS2 Fragmentation Analysis")  
        ms2_data = extract_ms2_data(experiment)  
          
        if ms2_data:  
            # Let user pick a spectrum based on retention time  
            options = [f"RT: {spec['rt']:.2f}, Precursor m/z: {spec['precursor']:.4f}" for spec in ms2_data]  
            choice = st.selectbox("Select an MS2 spectrum to view:", options)  
              
            # Get the selected spectrum  
            selected_index = options.index(choice)  
            selected_ms2 = ms2_data[selected_index]  
              
            # Check if the selected MS2 spectrum has valid data  
            if len(selected_ms2['mz']) > 0 and len(selected_ms2['intensity']) > 0:  
                # Plot MS2 spectrum  
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
                  
                # Display data table  
                df = pd.DataFrame({  
                    'm/z': selected_ms2['mz'],  
                    'Intensity': selected_ms2['intensity']  
                })  
                st.dataframe(df)  
            else:  
                st.warning("The selected MS2 spectrum does not contain any peaks.")  
        else:  
            st.warning("No MS2 spectra found in the file.")  
          
        # MGF Matching  
        st.header("MGF Matching")  
        uploaded_mgf = st.file_uploader("Upload an MGF file for fragment matching", type=["mgf"])  
          
        if uploaded_mgf is not None and ms2_data:  
            try:  
                mgf_content = uploaded_mgf.read().decode('utf-8')  
                records = parse_mgf(mgf_content)  
                  
                if records:  
                    # Use the first MS2 spectrum for matching  
                    ms2_spec = ms2_data[0]  
                    record = records[0]  
                      
                    if len(record.get("mz_array", [])) > 0:  
                        # Define tolerance in Da for matching  
                        tolerance = st.number_input("Set mass tolerance (Da):", value=0.5, step=0.1)  
                          
                        match_results = match_fragments(ms2_spec, record, tolerance)  
                          
                        if match_results['ms2_indices']:  
                            st.success("Matching fragments found!")  
                              
                            # Plot the matching on the MS2 spectrum plot  
                            match_marker = dict(color="#EB3424", size=10, symbol="circle-open")  
                            fig_match = go.Figure()  
                            fig_match.add_trace(go.Scatter(  
                                x=ms2_spec['mz'],  
                                y=ms2_spec['intensity'],  
                                mode="markers",  
                                marker=dict(color="#2563EB", size=8),  
                                name="MS2 peaks"  
                            ))  
                            # Add markers for matched fragments  
                            fig_match.add_trace(go.Scatter(  
                                x=[ms2_spec['mz'][i] for i in match_results['ms2_indices']],  
                                y=[ms2_spec['intensity'][i] for i in match_results['ms2_indices']],  
                                mode="markers",  
                                marker=match_marker,  
                                name="Matched fragments"  
                            ))  
                            fig_match.update_layout(  
                                title={"text":"MS2 Spectrum with Matching Fragments", "pad": {"t":15}},  
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
                        st.error("No valid m/z values found in the MGF file.")  
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
